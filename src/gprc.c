/*
  libgpr - a library for genetic programming
  Copyright (C) 2013  Bob Mottram <bob@sluggish.dyndns.org>

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  3. Neither the name of the University nor the names of its contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.
  .
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE HOLDERS OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "gprc.h"


/* returns the value of an actuator */
static float gprc_get_ADF_module_actuator(gprc_function * f, int index,
										  int module,
										  int rows, int columns,
										  int sensors)
{
	float * state = f->genome[module].state;
	return state[sensors + (rows*columns) + index];
}

/* returns the number of sensors for the given ADF_module */
int gprc_get_sensors(int ADF_module, int sensors)
{
	if (ADF_module == 0) {
		return sensors;
	}
	return GPRC_MAX_ADF_MODULE_SENSORS;
}

/* get the number of actuators for the given ADF_module */
int gprc_get_actuators(int ADF_module, int actuators)
{
	if (ADF_module == 0) {
		return actuators;
	}
	return 1;
}

/* returns the number of arguments that a function has */
static int gprc_function_args(int function_type,
							  float value,
							  int connections_per_gene,
							  int argc)
{
	/* no arguments */
	if (function_type==GPR_FUNCTION_VALUE) {
		return 0;
	}

	/* single argument functions */
	if ((function_type==GPR_FUNCTION_NEGATE) ||
		(function_type==GPR_FUNCTION_WEIGHT) ||
		(function_type==GPR_FUNCTION_FLOOR) ||
		(function_type==GPR_FUNCTION_EXP) ||
		(function_type==GPR_FUNCTION_SQUARE_ROOT) ||
		(function_type==GPR_FUNCTION_ABS) ||
		(function_type==GPR_FUNCTION_SINE) ||
		(function_type==GPR_FUNCTION_ARCSINE) ||
		(function_type==GPR_FUNCTION_COSINE) ||
		(function_type==GPR_FUNCTION_NOOP1) ||
		(function_type==GPR_FUNCTION_NOOP2) ||
		(function_type==GPR_FUNCTION_NOOP3) ||
		(function_type==GPR_FUNCTION_NOOP4) ||
		(function_type==GPR_FUNCTION_ARCCOSINE)) {
		return 1;
	}
	/* two argument functions */
	if ((function_type==GPR_FUNCTION_CUSTOM) ||
		(function_type==GPR_FUNCTION_MODULUS) ||
		(function_type==GPR_FUNCTION_DIVIDE) ||
		(function_type==GPR_FUNCTION_GREATER_THAN) ||
		(function_type==GPR_FUNCTION_LESS_THAN) ||
		(function_type==GPR_FUNCTION_EQUALS) ||
		(function_type==GPR_FUNCTION_AND) ||
		(function_type==GPR_FUNCTION_OR) ||
		(function_type==GPR_FUNCTION_XOR) ||
		(function_type==GPR_FUNCTION_NOT) ||
		(function_type==GPR_FUNCTION_POW) ||	
		(function_type==GPR_FUNCTION_COPY_FUNCTION) ||	
		(function_type==GPR_FUNCTION_COPY_CONSTANT) ||	
		(function_type==GPR_FUNCTION_COPY_CONNECTION1) ||	
		(function_type==GPR_FUNCTION_COPY_CONNECTION2) ||	
		(function_type==GPR_FUNCTION_COPY_CONNECTION3) ||	
		(function_type==GPR_FUNCTION_COPY_CONNECTION4) ||	
		(function_type==GPR_FUNCTION_GET) ||
		(function_type==GPR_FUNCTION_SET) ||
		(function_type==GPR_FUNCTION_COPY_BLOCK)) {
		return 2;
	}

	if ((function_type==GPR_FUNCTION_ADD) ||
		(function_type==GPR_FUNCTION_MULTIPLY) ||
		(function_type==GPR_FUNCTION_SIGMOID) ||
		(function_type==GPR_FUNCTION_AVERAGE) ||
		(function_type==GPR_FUNCTION_MIN) ||
		(function_type==GPR_FUNCTION_MAX) ||
		(function_type==GPR_FUNCTION_HEBBIAN) ||
		(function_type==GPR_FUNCTION_SUBTRACT)) {
		return 1 + (abs((int)value)%(connections_per_gene-1));
	}

	if (function_type==GPR_FUNCTION_ADF) {
		argc = 1 + abs(((int)argc)%GPRC_MAX_ADF_MODULE_SENSORS);
		if (argc >= connections_per_gene) {
			argc = connections_per_gene-1;
		}
		return argc;
	}

	return connections_per_gene;
}

/* clears the state of an individual */
void gprc_clear_state(gprc_function * f,
					  int rows, int columns,
					  int sensors, int actuators)
{
	for (int m = 0; m < f->ADF_modules+1; m++) {
		memset((void*)f->genome[m].state, '\0',
			   ((rows*columns)+
				gprc_get_sensors(m, sensors)+
				gprc_get_actuators(m, actuators))*
			   sizeof(float));	
	}
}

/* clears the used flags within a module */
void gprc_clear_used_module(gprc_function * f,
							int module,
							int rows, int columns,
							int sensors, int actuators)
{
	memset((void*)f->genome[module].used,'\0',
		   ((rows*columns)+
			gprc_get_sensors(module, sensors)+
			gprc_get_actuators(module, actuators))*
		   sizeof(unsigned char));	
}

/* clear the usage arrays */
void gprc_clear_used(gprc_function * f,
					 int rows, int columns,
					 int sensors, int actuators)
{
	for (int m = 0; m < f->ADF_modules+1; m++) {
		gprc_clear_used_module(f, m,
							   rows, columns,
							   sensors, actuators);
	}
}

/* initialise sensor sources for the given system */
void gprc_init_sensor_sources(gprc_system * system,
							  int no_of_sensor_sources,
							  unsigned int * random_seed)
{
	int i,j,k;
	gprc_function * f;
	gprc_population * population;

	if (no_of_sensor_sources <= 0) return;

	for (i = 0; i < system->size; i++) {
		/* for every island in the system */
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			f = &population->individual[j];
			f->no_of_sensor_sources = no_of_sensor_sources;
			/* create the array */
			f->sensor_source =
				(int*)malloc(population->sensors*sizeof(int));
			/* set random values */
			for (k = 0; k < population->sensors; k++) {
				f->sensor_source[k] =
					rand_num(random_seed)%no_of_sensor_sources;
			}
		}
	}
}							 

/* initialise actuator_detinations for the given system */
void gprc_init_actuator_destinations(gprc_system * system,
									 int no_of_actuator_destinations,
									 unsigned int * random_seed)
{
	int i,j,k;
	gprc_function * f;
	gprc_population * population;

	if (no_of_actuator_destinations <= 0) return;

	for (i = 0; i < system->size; i++) {
		/* for every island in the system */
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			f = &population->individual[j];
			f->no_of_actuator_destinations =
				no_of_actuator_destinations;
			/* create the array */
			f->actuator_destination =
				(int*)malloc(population->actuators*sizeof(int));
			/* initial random values */
			for (k = 0; k < population->actuators; k++) {
				f->actuator_destination[k] =
					rand_num(random_seed)%no_of_actuator_destinations;
			}
		}
	}
}							 

/* Clears all gene values */
static void gprc_clear_genome(gprc_function * f,
							  int rows, int columns,
							  int actuators,
							  int connections_per_gene)
{
	int m;
	float * gene;

	for (m = 0; m < f->ADF_modules+1; m++) {
		gene = f->genome[m].gene;
		memset((void*)gene,'\0',
			   ((rows*columns*GPRC_GENE_SIZE(connections_per_gene)) +
				gprc_get_actuators(m, actuators))*sizeof(float));
	}
}

/* initialize an individual */
void gprc_init(gprc_function * f,
			   int rows, int columns, int sensors, int actuators,
			   int connections_per_gene,
			   int ADF_modules,
			   unsigned int * random_seed)
{
	int m, sens, act;

	/* allocate arrays */
	f->ADF_modules = ADF_modules;
	for (m = 0; m < ADF_modules+1; m++) {
		sens = gprc_get_sensors(m, sensors);
		act = gprc_get_actuators(m, actuators);
		f->genome[m].gene =
			(float*)malloc(((rows*columns*
							 GPRC_GENE_SIZE(connections_per_gene)) +
							act)*sizeof(float));
		f->genome[m].state =
			(float*)malloc(((rows*columns) + sens + act)*
						   sizeof(float));
		f->genome[m].used =
			(unsigned char*)malloc(((rows*columns) + sens + act)*
								   sizeof(unsigned char));
	}

	/* clear the state */
	gprc_clear_state(f, rows, columns, sensors, actuators);

	/* clear the used genes */
	gprc_clear_used(f, rows, columns, sensors, actuators);

	/* clear the genome */
	gprc_clear_genome(f, rows, columns,
					  actuators, connections_per_gene);

	/* don't do anything with sensor sources
	   or actuator destinations.  We can enable this later if it
	   is needed */
	f->no_of_sensor_sources = 0;
	f->sensor_source = 0;
	f->no_of_actuator_destinations = 0;
	f->actuator_destination = 0;

	/* intial random seed */
	f->random_seed = rand_num(random_seed);

	f->age = 0;

	f->temp_genes = (int*)malloc(GPRC_MAX_ADF_GENES*3*sizeof(int));
}

/* deallocate memory for an individual */
void gprc_free(gprc_function * f)
{
	for (int m = 0; m < f->ADF_modules+1; m++) {
		free(f->genome[m].gene);
		free(f->genome[m].state);
		free(f->genome[m].used);
	}

	if (f->no_of_sensor_sources>0) {
		free(f->sensor_source);
	}

	if (f->no_of_actuator_destinations>0) {
		free(f->actuator_destination);
	}
	free(f->temp_genes);
}

/* does the given ADF_module contain an ADF call to
   the given ADF_module number */
int gprc_contains_ADFs(gprc_function * f,
					   int ADF_module, int call_ADF_module,
					   int rows, int columns,
					   int connections_per_gene,
					   int sensors)
{
	int n=0;
	int function_type, call_ADF_module2;
	int index,sens = gprc_get_sensors(ADF_module,sensors);

	for (index = sens; index < sens + (rows*columns); index++,
			 n += GPRC_GENE_SIZE(connections_per_gene)) {
		/* has been traced */
		if (f->genome[ADF_module].used[index] == 1) {
			function_type =
				(int)f->genome[ADF_module].gene[n];
			if (function_type == GPR_FUNCTION_ADF) {
				call_ADF_module2 =
					1 + (abs((int)f->genome[ADF_module].gene[n+GPRC_GENE_CONSTANT])%
						 f->ADF_modules);
				if (call_ADF_module == call_ADF_module2) {
					return n;
				}
			}
		}
	}
	return -1;
}

/* detects a subgraph up to a certain maximum depth */
int gprc_get_subgraph(gprc_function * f,
					  int ADF_module,
					  int index, int parent_index, int connection,
					  int rows, int columns,
					  int connections_per_gene,
					  int sensors,
					  int depth,
					  int max_depth, int max_genes,
					  int * genes, int * no_of_genes,
					  int * no_of_inputs,
					  int random_termination)
{
	int n=index,function_type,min,max,connection_index;
	float value;
	int sens = gprc_get_sensors(ADF_module,sensors);

	if ((*no_of_inputs >= GPRC_MAX_ADF_MODULE_SENSORS) ||
		(*no_of_genes >= max_genes)) {
		if (*no_of_genes > 0) {
			genes[((*no_of_genes-1) * 3)+1] =
				-(genes[((*no_of_genes-1) * 3)+1]+1);
		}
		return 0;
	}

	function_type = GPR_FUNCTION_VALUE;
	if (index >= sens) {
		n = (index - sens)*GPRC_GENE_SIZE(connections_per_gene);
		function_type = (int)f->genome[ADF_module].gene[n];
	}

	if ((function_type == GPR_FUNCTION_ADF) ||
		(depth >= max_depth) || (index < sens) ||
		((random_termination > 0) &&
		 ((rand_num(&f->random_seed)%10000) < 1000))) {
		*no_of_inputs = *no_of_inputs + 1;
		/* index of the input */
		genes[*no_of_genes * 3] = index;
		/* negative value indicates that this is an input */
		genes[(*no_of_genes * 3)+1] = -(parent_index+1);
		/* the connection index of the parent gene */
		genes[(*no_of_genes * 3)+2] = connection;
		*no_of_genes = *no_of_genes + 1;
		return 0;
	}

	genes[*no_of_genes * 3] = index;
	genes[(*no_of_genes * 3)+1] = parent_index;
	genes[(*no_of_genes * 3)+2] = connection;
	*no_of_genes = *no_of_genes + 1;

	value = f->genome[ADF_module].gene[n+GPRC_GENE_CONSTANT];

	min = 0;
	max = gprc_function_args(function_type,
							 value, connections_per_gene,
							 (int)f->genome[ADF_module].gene[n+GPRC_INITIAL]);

	if (function_type==GPR_FUNCTION_ADF) {
		min = 1;
		max++;
	}

	for (connection = min; connection < max; connection++) {
		/* get the prior connection */
		connection_index =
			(int)f->genome[ADF_module].gene[n + GPRC_INITIAL +
											connection];
		/* recursively traverse the graph */
		gprc_get_subgraph(f, ADF_module,
						  connection_index, index, connection,
						  rows, columns,
						  connections_per_gene,
						  sensors, depth+1,
						  max_depth, max_genes,
						  genes, no_of_genes,
						  no_of_inputs,
						  random_termination);
	}

	return 1;
}

/* returns the index of an unused ADF,
   or -1 if there are no unused ADFs */
static int gprc_get_unused_ADF(gprc_function * f,
							   int rows, int columns,
							   int connections_per_gene,
							   int sensors)
{
	int m;

	for (m = 0; m < f->ADF_modules; m++) {
		if (gprc_contains_ADFs(f,0,m+1,
							   rows, columns,
							   connections_per_gene,
							   sensors)==-1) {
			return m+1;
		}
	}
	return -1;
}

/* returns the number of arguments of the given ADF */
int get_ADF_args(gprc_function * f, int ADF_module)
{
	int i,ctr=0;
	unsigned char * used = f->genome[ADF_module].used;

	for (i=0;i<GPRC_MAX_ADF_MODULE_SENSORS;i++) {
		if (used[i]==1) ctr++;
	}
	return ctr;
}

/* move a subgraph from the main program into one of the ADF modules */
static void gprc_move_code_to_ADF(gprc_function * f,
								  int ADF_module,
								  int call_ADF_module,
								  int rows, int columns,
								  int connections_per_gene,
								  int sensors,
								  int no_of_genes,
								  int no_of_inputs)
{
	int i,index,parent_index,n,c,argc=0,parent_connection;
	int function_type;
	int sens = gprc_get_sensors(ADF_module,sensors);
	int call_sens = gprc_get_sensors(call_ADF_module,sensors);
	float value;

	/* clear the used flags */
	gprc_clear_used_module(f, call_ADF_module,
						   rows, columns,
						   sensors, 1);

	/* set which inputs are used */
	memset((void*)f->genome[call_ADF_module].used,'\0',call_sens);		 
	for (i = 0; i < no_of_inputs; i++) {
		f->genome[call_ADF_module].used[i] = 1;
	}

	for (i = 0; i < no_of_genes; i++) {
		/* index of the parent gene */
		parent_index = f->temp_genes[(i*3) + 1];
		/* is this an input? */
		if (parent_index < 0) {
			parent_index = (-parent_index)-1;

			n = (parent_index-sens) *
				GPRC_GENE_SIZE(connections_per_gene);

			parent_connection = f->temp_genes[(i*3) + 2];

			f->genome[call_ADF_module].gene[n + GPRC_INITIAL +
											parent_connection] =
				argc++;

			continue;
		}
		index = f->temp_genes[i*3];
		n = (index-sens) * GPRC_GENE_SIZE(connections_per_gene);
		if (i==0) {
			/* set the actuator source */
			f->genome[call_ADF_module].gene[rows*columns*
											GPRC_GENE_SIZE(connections_per_gene)] =
				index-sens+call_sens;
			f->genome[call_ADF_module].used[(rows*columns) +
											call_sens] = 1;
		}

		f->genome[call_ADF_module].used[index-sens+call_sens] = 1;

		/* copy function */
		function_type = f->genome[ADF_module].gene[n];
		f->genome[call_ADF_module].gene[n] = function_type;
		/* copy value */
		value = f->genome[ADF_module].gene[n+GPRC_GENE_CONSTANT];
		f->genome[call_ADF_module].gene[n+GPRC_GENE_CONSTANT] = value;
		/* copy connections */
		for (c = 0; c < connections_per_gene; c++) {
			f->genome[call_ADF_module].gene[n+GPRC_INITIAL+c] =
				f->genome[ADF_module].gene[n+GPRC_INITIAL+c]-
				sens+call_sens;
		}
	}
#ifdef DEBUG
	assert(argc == no_of_inputs);
#endif
}

/* update te number of ADF arguments */
static void gprc_update_ADF_arguments(gprc_function * f,
									  int ADF_module,
									  int call_ADF_module,
									  int rows, int columns,
									  int connections_per_gene,
									  int sensors)
{
	int argc, n=0, function_type, call_ADF_module2;
	int previous_values, row, col;
	float * gene = f->genome[ADF_module].gene;

	/* get the number of arguments for the ADF */
	argc = get_ADF_args(f, call_ADF_module);
	if (argc >= connections_per_gene) {
		argc = connections_per_gene-1;
		printf("Number of ADF arguments greater than " \
			   "connections per gene\n%d %d\n",
			   argc, connections_per_gene);
	}

	for (col = 0; col < columns; col++) {
		for (row = 0; row < rows; row++,
				 n += GPRC_GENE_SIZE(connections_per_gene)) {
			function_type = (int)gene[n];
			if (function_type == GPR_FUNCTION_ADF) {
				/* index of the ADF_module being called */
				call_ADF_module2 =
					1 + (abs((int)gene[n+GPRC_GENE_CONSTANT])%
						 f->ADF_modules);
				if (call_ADF_module == call_ADF_module2) {
					if (argc <= 0) {
						gene[n] = GPR_FUNCTION_VALUE;
						previous_values =
							gprc_get_sensors(ADF_module,
											 sensors) + (col*rows);
						gene[n+GPRC_INITIAL] =
							rand_num(&f->random_seed)%previous_values;
					}
					else {
						/* set the number of arguments
						   within the ADF gene */
						gene[n+GPRC_INITIAL] = argc-1;
					}
				}
			}
		}
	}
}

/* Update the connections for each ADF module */
static void gprc_update_ADF_modules(gprc_function * f,
									int rows, int columns,
									int connections_per_gene,
									int sensors)
{
#pragma omp parallel for
	for (int m = 1; m < f->ADF_modules+1; m++) {
		gprc_update_ADF_arguments(f, 0, m,
								  rows, columns,
								  connections_per_gene,
								  sensors);
	}
}

/* remove code and replace it with an ADF */
static void gprc_remove_code(gprc_function * f,
							 int ADF_module,
							 int call_ADF_module,
							 int rows, int columns,
							 int connections_per_gene,
							 int sensors,
							 int no_of_genes,
							 int no_of_inputs)
{
	int i,n=0,index,parent_index;
	int sens = gprc_get_sensors(ADF_module,sensors);
	int argc = 0;

	for (i = 0; i < no_of_genes; i++) {
		index = f->temp_genes[i*3];
		parent_index = f->temp_genes[(i*3) + 1];

		if (i == 0) {
			n = (index - sens) * GPRC_GENE_SIZE(connections_per_gene);
			f->genome[ADF_module].gene[n+GPRC_GENE_FUNCTION_TYPE] =
				GPR_FUNCTION_ADF;
			f->genome[ADF_module].gene[n+GPRC_GENE_CONSTANT] =
				call_ADF_module-1;
			f->genome[ADF_module].gene[n+GPRC_INITIAL] = 0;
		}
		else {
			if (parent_index < 0) {
				if (argc < connections_per_gene-1) {
					f->genome[ADF_module].gene[n+1+GPRC_INITIAL+argc]=
						index;
					argc++;
					f->genome[ADF_module].gene[n+GPRC_INITIAL] = argc;
				}
			}
			else {
				f->genome[ADF_module].used[index] = 0;
			}
		}
	}
	if (f->genome[ADF_module].gene[n+GPRC_INITIAL] > 0) {
		f->genome[ADF_module].gene[n+GPRC_INITIAL] -= 1;
	}
}

/* returns the index of a random used gene */
static int gprc_get_used_gene(gprc_function * f,
							  int ADF_module,
							  int rows, int columns,
							  int sensors)
{
	int index,idx,ctr=0;
	int no_of_used_genes=0;
	int sens = gprc_get_sensors(ADF_module,sensors);

	for (index = sens; index < sens + (rows*columns); index++) {
		if (f->genome[ADF_module].used[index] == 1) {
			no_of_used_genes++;
		}
	}

	/* no used genes were found */
	if (no_of_used_genes==0) return -1;
	
	/* pick a gene at random */
	idx = rand_num(&f->random_seed)%no_of_used_genes;

	for (index = sens; index < sens + (rows*columns); index++) {
		if (f->genome[ADF_module].used[index] == 1) {
			if (ctr==idx) return index;
			ctr++;
		}
	}

	return -1;
}

/* the purpose of this is to discover which functions within
   the grid are actually used as part of the input -> output
   transformation. */
static void gprc_used_genes(gprc_ADF_module * f,
							int rows, int columns,
							int connections_per_gene,
							int sensors, int actuators)
{
    int row,col,index=0,index2,c,connection_index,ctr,prev_ctr,n;
	int min,max,function_type;
	int array_bytes =
		(sensors+actuators+(rows*columns))*sizeof(unsigned char);

	/* clear the array */
	memset((void*)f->used,'\0',array_bytes);

	/* mark actuators as traced */
	for (index=0; index < actuators; index++) {
		f->used[sensors+(rows*columns)+index]=1;
	}

	/* propagate backwards */
	ctr = 0;
	prev_ctr = -1;
	while (ctr > prev_ctr) {
		prev_ctr = ctr;
		index = 0;
		n = 0;
		for (col = 0; col < columns; col++) {
			for (row = 0; row < rows; row++, index++,
					 n += GPRC_GENE_SIZE(connections_per_gene)) {
				/* has been traced */
				if (f->used[index + sensors] == 1) {

					function_type = (int)f->gene[n];
					min=0;
					max =
						gprc_function_args(function_type,
										   f->gene[n+GPRC_GENE_CONSTANT],
										   connections_per_gene,
										   (int)f->gene[n+GPRC_INITIAL]);

					if (function_type==GPR_FUNCTION_ADF) {
						min=1;
						max++;
					}

					for (c = min; c < max; c++) {
						/* get the prior connection */
						connection_index =
							(int)f->gene[n + GPRC_INITIAL + c];
						/* prior has not been traced */
						if (f->used[connection_index] == 0) {
							/* mark the prior as traced */
							f->used[connection_index] = 1;
							ctr++;
						}
					}
				}
			}
		}
		for (index2 = 0; index2 < actuators; index2++, n++) {
			connection_index = (int)f->gene[n];
			/* prior has not been traced */
			if (f->used[connection_index] == 0) {
				/* mark the prior as traced */
				f->used[connection_index]=1;
				ctr++;
			}
		}
	}
}

/* the purpose of this is to discover which functions within
   the grid are actually used as part of the input -> output
   transformation. */
void gprc_used_functions(gprc_function * f,
						 int rows, int columns,
						 int connections_per_gene,
						 int sensors, int actuators)
{
	/*#pragma omp parallel for*/
	for (int m = 0; m < f->ADF_modules+1; m++) {
		gprc_used_genes(&f->genome[m],
						rows, columns,
						connections_per_gene,
						gprc_get_sensors(m, sensors),
						gprc_get_actuators(m,actuators));
	}
}							  

/* Tries to convert code within the given module into
   an automatically defined function */
int gprc_compress_ADF(gprc_function * f,
					  int ADF_module,
					  int start_index,
					  int rows, int columns,
					  int connections_per_gene,
					  int sensors, int actuators,
					  float min_value, float max_value,
					  int max_depth,
					  int random_termination)
{
	int n, call_ADF_module, index, function_type;
	int no_of_genes = 0;
	int no_of_inputs = 0;

	if (f->ADF_modules == 0) return -1;

	/* does an unused ADF exist? */
	call_ADF_module =
		gprc_get_unused_ADF(f, rows, columns,
							connections_per_gene,
							sensors);

	if (start_index > -1) {
		index = start_index;
	}
	else {
		/* get a used gene at random */
		index = gprc_get_used_gene(f, ADF_module,
								   rows, columns,
								   sensors);
	}

	/* was a gene found? */
	if (index == -1) {
		return -2;
	}

	n = (index - gprc_get_sensors(ADF_module,sensors)) *
		GPRC_GENE_SIZE(connections_per_gene);
	function_type = f->genome[ADF_module].gene[n];
	if ((function_type == GPR_FUNCTION_ADF) ||
		(function_type == GPR_FUNCTION_VALUE) ||
		(function_type == GPR_FUNCTION_NONE)) {
		/* don't try to replace existing ADFs or values */
		return 0;
	}

	if (call_ADF_module == -1) {
		/* use one of the existing ADFs */		
		f->genome[ADF_module].gene[n+GPRC_GENE_FUNCTION_TYPE] =
			GPR_FUNCTION_ADF;
		f->genome[ADF_module].gene[n+GPRC_GENE_CONSTANT] =
			rand_num(&f->random_seed)%(f->ADF_modules);

		f->genome[ADF_module].gene[n+GPRC_INITIAL] =
			rand_num(&f->random_seed)%GPRC_MAX_ADF_MODULE_SENSORS;
		if (f->genome[ADF_module].gene[n+GPRC_INITIAL] >=
			connections_per_gene-1) {
			f->genome[ADF_module].gene[n+GPRC_INITIAL] =
				connections_per_gene-1;
		}
		return 0;
	}

	/* extract a subgraph beginning with this gene */
	gprc_get_subgraph(f, ADF_module,
					  index, 0, 0,
					  rows, columns,
					  connections_per_gene,
					  sensors,
					  0, max_depth, GPRC_MAX_ADF_GENES,
					  f->temp_genes, &no_of_genes,
					  &no_of_inputs,random_termination);

	if (no_of_genes < GPRC_MIN_ADF_GENES) {
		return -3;
	}
	if (no_of_genes >= GPRC_MAX_ADF_GENES) {
		return -4;
	}
	if (no_of_inputs >= GPRC_MAX_ADF_MODULE_SENSORS) {
		return -5;
	}
	if (no_of_inputs < 1) {
		return -6;
	}

	/* copy subgraph code into the new ADF module */
	gprc_move_code_to_ADF(f, ADF_module, call_ADF_module,
						  rows, columns,
						  connections_per_gene,
						  sensors, no_of_genes, no_of_inputs);

	/* remove subgraph code from the original module */	
	gprc_remove_code(f, ADF_module, call_ADF_module,
					 rows, columns, connections_per_gene,
					 sensors, no_of_genes, no_of_inputs);

	return 0;
}

/* do connections point to the same location?
   If yes then return the array index */
static int gprc_same_connections(float * gene,
								 int connections_per_gene)
{
	int i,j;

	/* one connection */
	if (connections_per_gene<=1) return 0;

	/* multiple connections */
	for (i = 0; i < connections_per_gene; i++) {
		for (j = i+1; j < connections_per_gene; j++) {
			if ((int)gene[GPRC_INITIAL+i] ==
				(int)gene[GPRC_INITIAL+j]) return i;
		}
	}
	return -1;
}

/* removes all ADF functions from the main program */
void gprc_remove_ADFs(gprc_function * f,
					  int rows, int columns,
					  int connections_per_gene)
{
	int index,n,function_type,m;
	float * gene;

	for (m = 0; m < f->ADF_modules+1; m++) {
		gene = f->genome[m].gene;
		n = 0;
		for (index = 0; index < rows*columns; index++,
				 n += GPRC_GENE_SIZE(connections_per_gene)) {
			function_type = (int)gene[n];
			if (function_type == GPR_FUNCTION_ADF) {
				gene[n] = GPR_FUNCTION_VALUE;
			}
		}
	}
}

void gprc_valid_ADFs(gprc_function * f,
					 int rows, int columns,
					 int connections_per_gene,
					 int sensors,
					 float min_value, float max_value)
{
	int m, n, function_type, ADF_module_index;
	int new_connection,row,col,previous_values;

	/*for (m = 0; m < f->ADF_modules+1; m++) {*/
	for (m = 0; m < 1; m++) {
		n=0;
		for (col = 0; col < columns; col++) {
			previous_values =
				gprc_get_sensors(m,sensors) + (col*rows);
			for (row = 0; row < rows; row++,
					 n += GPRC_GENE_SIZE(connections_per_gene)) {
				/* for each gene */				
				function_type = (int)f->genome[m].gene[n];
				if (function_type==0) {
					function_type = GPR_FUNCTION_VALUE;
				}
				
				if (function_type == GPR_FUNCTION_ADF) {
					if ((m > 0) || (f->ADF_modules == 0)) {
						f->genome[m].gene[n] = GPR_FUNCTION_VALUE;

						f->genome[m].gene[n+GPRC_GENE_CONSTANT] =
							gpr_random_value(min_value, max_value,
											 &f->random_seed);

						new_connection =
							rand_num(&f->random_seed)%
							previous_values;
						f->genome[m].gene[n+GPRC_INITIAL] =
							new_connection;
					}
					else {
						ADF_module_index =
							1 + (abs((int)f->genome[m].gene[n+GPRC_GENE_CONSTANT])%
								 f->ADF_modules);
						f->genome[m].gene[n+GPRC_GENE_CONSTANT] =
							ADF_module_index-1;
						/*
						  argc = get_ADF_args(f, ADF_module_index);
						  f->genome[m].gene[n+GPRC_INITIAL] = argc-1;
						  if (f->genome[m].gene[n+GPRC_INITIAL] < 0) {
						  f->genome[m].gene[n] = GPR_FUNCTION_VALUE;
						  }
						*/
					}
				}
			}
		}
	}
}								 

/* for functions which have two inputs make sure that the
   inputs are from different sources in some cases */
static void gprc_ADF_valid_logical_operators(gprc_ADF_module * f,
											 int rows, int columns,
											 int connections_per_gene,
											 int sensors,
											 unsigned int * random_seed)
{
	int row,col,n=0,previous_values,function_type,index;
	int attempts,max;

	for (col = 0; col < columns; col++) {
		previous_values = (col*rows) + sensors;
		for (row = 0; row < rows; row++,
				 n += GPRC_GENE_SIZE(connections_per_gene)) {
			/* for each gene */
			function_type = (int)f->gene[n+GPRC_GENE_FUNCTION_TYPE];

			/* ignore ADFs */
			if (function_type==GPR_FUNCTION_ADF) continue;

			/* maximum number of arguments for this function */
			max = gprc_function_args(function_type,
									 f->gene[n+GPRC_GENE_CONSTANT],
									 connections_per_gene,
									 (int)f->gene[n+GPRC_INITIAL]);

			if (max <= 1) continue;

			/* look for connections which are the same */
			index = gprc_same_connections(&f->gene[n],max);
			attempts=0;
			while ((index>-1) && (attempts<5)) {
				/* change the connection */
				f->gene[n+GPRC_INITIAL+index] =
					rand_num(random_seed)%previous_values;
				index = gprc_same_connections(&f->gene[n],max);
				attempts++;
			}
			if ((attempts==5) && (max>2)) {
				f->gene[n+GPRC_GENE_CONSTANT] = 0;
			}
		}
	}
}

/* for functions which have two inputs make sure that the
   inputs are from different sources in some cases */
void gprc_valid_logical_operators(gprc_function * f,
								  int rows, int columns,
								  int connections_per_gene,
								  int sensors,
								  unsigned int * random_seed)
{
	int m;

	for (m = 0; m < f->ADF_modules+1; m++) {
		gprc_ADF_valid_logical_operators(&f->genome[m],
										 rows, columns,
										 connections_per_gene,
										 gprc_get_sensors(m,sensors),
										 random_seed);
	}
}

/* ensures that the output sources are unique */
void gprc_unique_outputs(gprc_function * f,
						 int rows, int columns,
						 int connections_per_gene,
						 int sensors, int actuators,
						 unsigned int * random_seed)
{
	int i,j,changes,attempts,m,act,sens;
	int n = rows*columns*GPRC_GENE_SIZE(connections_per_gene);
	float * gene;
	const int max_attempts = 20;

	for (m = 0; m < f->ADF_modules+1; m++) {
		gene = f->genome[m].gene;
		sens = gprc_get_sensors(m,sensors);
		act = gprc_get_actuators(m,actuators);
		if (act <= 1) continue;

		/* actuators should point to different genes */
		changes=1;
		attempts=0;
		while ((changes > 0) && (attempts<max_attempts)) {
			changes = 0;
			for (i = 0; i < act; i++) {
				for (j = i+1; j < act; j++) {
					if ((int)gene[n+i] == (int)gene[n+j]) {
						gene[n+i] =
							sens + (int)rand_num(random_seed)%
							(rows*columns);
						changes++;
						break;
					}
				}
			}
			attempts++;
		}
	}
}

/* forces the given individual to be valid */
void gprc_tidy(gprc_function * f,
			   int rows, int columns,
			   int sensors, int actuators,
			   int connections_per_gene,
			   float min_value, float max_value)
{
	/* check that output connections are unique */
	gprc_unique_outputs(f, rows, columns,
						connections_per_gene,
						sensors, actuators,
						&f->random_seed);

	/* ensure that any logical operators have valid inputs */
	gprc_valid_logical_operators(f, rows, columns,
								 connections_per_gene,
								 sensors,
								 &f->random_seed);

	/* ensure that ADF calls are valid */
	gprc_valid_ADFs(f, rows, columns,
					connections_per_gene,
					sensors, min_value, max_value);

	/* update the used functions */
	gprc_used_functions(f, rows, columns,
						connections_per_gene,
						sensors, actuators);
}

/* creates an initial random state for an individual */
void gprc_random(gprc_function * f,
				 int rows, int columns,
				 int sensors, int actuators,
				 int connections_per_gene,
				 float min_value, float max_value,
				 int integers_only, unsigned int * random_seed,
				 int * instruction_set, int no_of_instructions)
{
	int col,row,n,i,j,w,previous_values,m,act,function_type;
	float * gene;

	/* for each ADF_module */
	for (m = 0; m < f->ADF_modules+1; m++) {
		/* get the genome for this ADF_module */
		gene = f->genome[m].gene;
		n=0;
		for (col = 0; col < columns; col++) {
			previous_values = (col*rows) + gprc_get_sensors(m,sensors);
			for (row = 0; row < rows; row++) {
				/* instruction type */
				function_type =
					gpr_random_function(instruction_set,
										no_of_instructions,
										random_seed);
				if (function_type == GPR_FUNCTION_ADF) {
					function_type = GPR_FUNCTION_VALUE;
				}
				gene[n++] = function_type;

				/* value */
				if (integers_only <= 0) {
					gene[n++] = gpr_random_value(min_value, max_value,
												 random_seed);
				}
				else {
					gene[n++] = gpr_random_value(min_value, max_value,
												 random_seed);
				}
				/* other values */
				for (j = 2; j < GPRC_INITIAL; j++) {
					gene[n++] = 0;
				}
				/* input indexes */
				for(i = 0; i < connections_per_gene; i++) {
					gene[n++] = rand_num(random_seed)%previous_values;
				}
				/* connection weights */
				for (w = 1; w < GPRC_WEIGHTS_PER_CONNECTION; w++) {
					for(i = 0; i < connections_per_gene; i++) {
						if (integers_only<=0) {
							gene[n++] =
								gpr_random_value(GPRC_MIN_WEIGHT,
												 GPRC_MAX_WEIGHT,
												 random_seed);
						}
						else {
							gene[n++] =
								gpr_random_value(GPRC_MIN_WEIGHT,
												 GPRC_MAX_WEIGHT,
												 random_seed);
						}
					}
				}
			}
		}

		/* random outputs */
		gene = f->genome[m].gene;
		act = gprc_get_actuators(m,actuators);
		for (i = 0; i < act; i++) {
			gene[n++] = gprc_get_sensors(m,sensors) +
				(int)rand_num(random_seed)%(rows*columns);
		}
	}

	gprc_tidy(f, rows, columns,
			  sensors, actuators,
			  connections_per_gene,
			  min_value, max_value);
}

/* prints the state to the console */
void print_gprc(gprc_function * f,
				int ADF_module,
				int rows, int columns,
				int sensors, int actuators,
				int connections_per_gene,
				int integers_only)
{
	int row,col,i,index=0,index2=0;
	float v;
	float * gene = f->genome[ADF_module].gene;
	float * state = f->genome[ADF_module].state;

	actuators = gprc_get_actuators(ADF_module, actuators);

	printf("\n\nSensors\n\n");

	for (i = 0; i <
			 gprc_get_sensors(ADF_module,sensors);
		 i++, index++) {
		v = state[index];
		if (v >= 0) {
			printf("  ");
		}
		else {
			printf(" -");
			v = -v;
		}		
		if (integers_only > 0) {
			printf("%05d", (int)v);
		}
		else {
			printf("%05d.%01d", (int)v,
				   (int)((v - (int)v)*10));
		}		
	}

	printf("\n\nSystem state\n\n");

	for (col = 0; col < columns; col++) {
		for (row = 0; row < rows; row++,index++) {
			v = state[index];
			if (v >= 0) {
				printf("  ");
			}
			else {
				printf(" -");
				v = -v;
			}		
			if (integers_only > 0) {
				printf("%05d", (int)v);
			}
			else {
				printf("%05d.%01d", (int)v,
					   (int)((v - (int)v)*10));
			}		
		}
		printf("\n");
	}

	printf("\nFunctions\n\n");

	index2 = 0;
	for (col = 0; col < columns; col++) {
		for (row = 0; row < rows; row++,index2++) {
			printf(" %02d",
				   (int)gene[index2*
							 GPRC_GENE_SIZE(connections_per_gene)]);
		}
		printf("\n");
	}

	printf("\nConnections\n\n");

	for (col = 0; col < columns; col++) {
		for (i = 0; i < connections_per_gene; i++) {
			for (row = 0; row < rows; row++) {
				index2 = (col*rows) + row;
				printf(" %05d",
					   (int)gene[index2*
								 GPRC_GENE_SIZE(connections_per_gene)+
								 GPRC_INITIAL+i]);
			}
			printf("\n");
		}
		printf("\n");
	}

	printf("\nActuators\n\n");
	for (i = 0; i < actuators; i++, index++) {
		v = state[index];
		if (v >= 0) {
			printf("  ");
		}
		else {
			printf(" -");
			v = -v;
		}		
		if (integers_only > 0) {
			printf("%05d", (int)v);
		}
		else {
			printf("%05d.%01d", (int)v,
				   (int)((v - (int)v)*10));
		}		
	}

	printf("\n\n");
}

void gprc_dot_label(gprc_function * f,
					int ADF_module,
					int rows, int columns,
					int connections_per_gene,
					int sensors, int actuators,
					char * sensor_names[],
					char * actuator_names[],
					FILE * fp)
{
	char name[256];
	int function_type,col,row,index=0,index2,used_ctr;
	int connection_index,n=0;
	float constant_value;
	float * gene;
	unsigned char * used;
	int sens = gprc_get_sensors(ADF_module,sensors);
	int offset =
		ADF_module*(((rows*columns) + sensors + actuators)+100);

	gene = f->genome[ADF_module].gene;
	used = f->genome[ADF_module].used;

	if (ADF_module > 0) {
		actuators = 1;
	}

	/* inputs */
	used_ctr=0;
	for (index = 0; index < sens; index++) {
		if (used[index] == 1) {
			if (ADF_module == 0) {
				fprintf(fp,"  a%d [label=\"%s\", shape=box];\n",
						index, sensor_names[index]);
			}
			else {
				fprintf(fp,
						"  a%d [label=\"ADF %d Sensor %d\", ",
						offset + index, ADF_module, used_ctr);
				fprintf(fp,"%s","shape=box];\n");
				used_ctr++;
			}
		}
	}

	/* grid functions */
	index = 0;
	for (col = 0; col < columns; col++) {
		for (row = 0; row < rows;
			 row++, index++, n+=GPRC_GENE_SIZE(connections_per_gene)) {
			if (used[sens + index] == 1) {
				function_type = (int)gene[n];
				constant_value = gene[n+GPRC_GENE_CONSTANT];

				sprintf(name,"Unknown %d", function_type);
				gpr_get_function_name(function_type,
									  constant_value,
									  1, name);
				fprintf(fp,"  a%d [label=\"%s\"];\n",
						offset + sens+index, name);
			}
		}
	}

	/* actuator names */
	for (index2 = 0; index2 < actuators; index2++, index++) {
		/* include the previous gene */
		connection_index = (int)gene[n + index2];
		if (used[connection_index] == 0) {
			if (connection_index >= sensors) {
				connection_index -= sensors;
				function_type =
					(int)gene[connection_index*
							  GPRC_GENE_SIZE(connections_per_gene)];
				constant_value =
					gene[(connection_index*
						  GPRC_GENE_SIZE(connections_per_gene)) + 1];

				sprintf(name,"Unknown %d", function_type);
				gpr_get_function_name(function_type,
									  constant_value,
									  1, name);

				fprintf(fp,"  a%d [label=\"%s\"];\n",
						offset + connection_index,name);
			}
			else {
				if (ADF_module == 0) {
					fprintf(fp,"  a%d [label=\"%s\"];\n",
							offset + connection_index,
							sensor_names[connection_index]);
				}
				else {
					fprintf(fp,"  a%d [label=\"Arg %d\"];\n",
							offset + connection_index,
							connection_index);
				}
			}
		}
		if (ADF_module == 0) {
			fprintf(fp,"  a%d [label=\"%s\", shape=box];\n",
					offset + sens + index, actuator_names[index2]);
		}
		else {
			fprintf(fp,"  a%d [label=\"Return\", shape=box];\n",
					offset + sens + index);
		}
	}
}

void gprc_dot_links(gprc_function * f,
					int ADF_module,
					int rows, int columns,
					int connections_per_gene,
					int sensors, int actuators,
					FILE * fp)
{
	int col,row,index=0,index2,c,connection_index,n=0;
	int min,max,function_type;
	float * gene;
	unsigned char * used;
	int sens = gprc_get_sensors(ADF_module,sensors);
	int offset =
		ADF_module*(((rows*columns) + sensors + actuators)+100);

	gene = f->genome[ADF_module].gene;
	used = f->genome[ADF_module].used;

	if (ADF_module > 0) {
		actuators = 1;
	}

	/* grid links */
	for (col = 0; col < columns; col++) {
		for (row = 0; row < rows; row++,
				 index++, n += GPRC_GENE_SIZE(connections_per_gene)) {
			if (used[sens + index] == 1) {

				function_type = (int)gene[n];
				min=0;
				max = gprc_function_args(function_type,
										 gene[n+GPRC_GENE_CONSTANT],
										 connections_per_gene,
										 (int)gene[n+GPRC_INITIAL]);

				if (function_type==GPR_FUNCTION_ADF) {
					min=1;
					max++;
				}

				for (c = min; c < max; c++) {
					/* prior connection */
					connection_index =
						(int)gene[n + GPRC_INITIAL + c];

					if (used[connection_index] == 1) {
						fprintf(fp,"  a%d -> a%d;\n",
								offset + connection_index,
								offset + sens + index);
					}
				}
			}
		}
	}

	/* actuator links */
	for (index2 = 0; index2 < actuators; index2++,index++,n++) {
		connection_index = (int)gene[n];
		fprintf(fp,"  a%d -> a%d;\n",
				offset + connection_index, offset + sens + index);
	}
}

void gprc_dot(gprc_function * f,
			  gprc_population * population,
			  char * sensor_names[],
			  char * actuator_names[],
			  FILE * fp)
{
	int m;

	/* ensure that all ADF values are valid */
	gprc_valid_ADFs(f,
					population->rows, population->columns,
					population->connections_per_gene,
					population->sensors,
					population->min_value,
					population->max_value);

	fprintf(fp,"%s","digraph graphname {\n");
	for (m = 0; m < f->ADF_modules+1; m++) {

		if ((m == 0) ||
			(gprc_contains_ADFs(f, 0, m,
								population->rows, population->columns,
								population->connections_per_gene,
								population->sensors)>-1)) {
			fprintf(fp,"%s","subgraph {\n");
			gprc_dot_label(f, m,
						   population->rows,
						   population->columns,
						   population->connections_per_gene,
						   population->sensors,
						   population->actuators,
						   sensor_names,
						   actuator_names,
						   fp);
			gprc_dot_links(f, m,
						   population->rows,
						   population->columns,
						   population->connections_per_gene,
						   population->sensors,
						   population->actuators,
						   fp);
			fprintf(fp,"%s","}\n");
		}
	}
	fprintf(fp,"%s","}\n");
}

/* sets teh source index for an output */
static void gprc_set_output_source(gprc_function * f,
								   int ADF_module,
								   int rows, int columns,
								   int connections_per_gene,
								   int output_index, int source_index)
{
	float * gene = f->genome[ADF_module].gene;
	int index =
		(rows*columns*GPRC_GENE_SIZE(connections_per_gene)) +
		output_index;

	gene[index] = source_index;
}

/* swaps the positions of two chromosomes */
static void gprc_swap_chromosomes(gprc_function * f,
								  int ADF_module,
								  int sensors,
								  int rows, int columns,
								  int connections_per_gene,
								  int chromosome_index1,
								  int chromosome_index2,
								  int chromosomes,
								  int shift)
{
	int col,col2,row,n,n2,i,i2,j,c;
	int start_row1, end_row1, start_row2, previous_values;
	float temp_value;
	float * gene = f->genome[ADF_module].gene;
	float * state = f->genome[ADF_module].state;

	if ((chromosome_index1 == chromosome_index2) &&
		(shift==0)) {
		return;
	}

	start_row1 = chromosome_index1 * rows / chromosomes;
	end_row1 = (chromosome_index1+1) * rows / chromosomes;
	start_row2 = chromosome_index2 * rows / chromosomes;

	for (col = 0; col < columns; col++) {
		col2 = col + shift;
		if (col2>=columns) col2 -= columns;
		if (col2<0) col2 += columns;

		if ((col2<col) || (chromosome_index1!=chromosome_index2)) {
			i = (col*rows) + start_row1;
			i2 = (col2*rows) + start_row2;
		}
		else {
			i = (col2*rows) + start_row1;
			i2 = (col*rows) + start_row2;
		}
		n = i * GPRC_GENE_SIZE(connections_per_gene);
		n2 = i2 * GPRC_GENE_SIZE(connections_per_gene);
		for (row = start_row1; row < end_row1; row++, i++, i2++,
				 n += GPRC_GENE_SIZE(connections_per_gene),
				 n2 += GPRC_GENE_SIZE(connections_per_gene)) {
			/* swap the gene */
			for (j = 0;
				 j < GPRC_GENE_SIZE(connections_per_gene); j++) {
				temp_value = gene[n+j];
				gene[n+j] = gene[n2+j];
				gene[n2+j] = temp_value;
			}
		
			/* swap the state */
			temp_value = state[i];
			state[i] = state[i2];
			state[i2] = temp_value;
		}
	}

	/* valid connections */
	n=0;
	for (col = 0; col < columns; col++) {
		previous_values =
			(col*rows) + gprc_get_sensors(ADF_module,sensors);
		for (row = 0; row < rows; row++,
				 n += GPRC_GENE_SIZE(connections_per_gene)) {
			for (c = 0; c < connections_per_gene; c++) {
				i = (int)gene[n+GPRC_INITIAL+c];
				if (i>=previous_values) {
					gene[n+GPRC_INITIAL+c] =
						rand_num(&f->random_seed)%previous_values;
				}
			}
		}
	}
}

/* randomly permutes connections within genes */
static void	gprc_mutation_permute(gprc_function * f,
								  int rows, int columns,
								  int connections_per_gene,
								  unsigned int * random_seed,
								  float prob)
{
	int no_of_mutations, i, index, con1, con2, col;
	int previous_values, m, n, function_type, min, max;
	float temp;
	float * gene;

	no_of_mutations =
		(int)((prob/(float)(f->ADF_modules+1))*rows*columns);
	for (m = 0; m < f->ADF_modules+1; m++) {
		gene = f->genome[m].gene;
		for (i = 0; i < no_of_mutations; i++) {
			/* pick a gene at random */
			index = rand_num(random_seed)%(rows*columns);
			col = (index / columns);
			if (col<=1) continue;
			previous_values = (col-1)*rows;
			n = index * GPRC_GENE_SIZE(connections_per_gene);
			function_type = (int)gene[n];
			/* pick connections */
			if (function_type != GPR_FUNCTION_ADF) {
				min = 0;
				max = connections_per_gene;
			}
			else {
				min = 1;
				max = get_ADF_args(f,m)+1;
			}
			if (max-min <= 0) continue;
			con1 = min + (rand_num(random_seed)%(max-min));
			con2 = min + (rand_num(random_seed)%(max-min));
			if (con1==con2) continue;


			if ((gene[n+GPRC_INITIAL+con1] < 0) ||
				(gene[n+GPRC_INITIAL+con1] >= previous_values)) {
				gene[n+GPRC_INITIAL+con1] =
					rand_num(random_seed)%previous_values;
			}
			if ((gene[n+GPRC_INITIAL+con2] < 0) ||
				(gene[n+GPRC_INITIAL+con2] >= previous_values)) {
				gene[n+GPRC_INITIAL+con2] =
					rand_num(random_seed)%previous_values;
			}

			/* swap */
			temp = gene[n+GPRC_INITIAL+con1];
			gene[n+GPRC_INITIAL+con1] = gene[n+GPRC_INITIAL+con2];
			gene[n+GPRC_INITIAL+con2] = temp;
		}
	}
}

/* mutates an individual */
void gprc_mutate(gprc_function * f,
				 int rows, int columns,
				 int sensors, int actuators,
				 int connections_per_gene,
				 int chromosomes,
				 float prob,
				 float chromosomes_prob,
				 float min_value, float max_value,
				 int integers_only,
				 int * instruction_set, int no_of_instructions)
{
	int no_of_mutations,function_type,call_ADF_module;
	int i,m,index,locn,col,new_connection;
	int sens, act, gene_index, conn_index;
	float * gene;
	int step = GPRC_GENE_SIZE(connections_per_gene);

	/* mutate sensor sources */
	if (f->no_of_sensor_sources > 0) {
		for (i = 0; i < sensors; i++) {
			if (rand_num(&f->random_seed)%10000 < prob*10000) {
				/* point mutation */
				f->sensor_source[i] =
					rand_num(&f->random_seed)%
					f->no_of_sensor_sources;
			}
		}
	}

	/* mutate actuator destinations */
	if (f->no_of_actuator_destinations > 0) {
		for (i = 0; i < actuators; i++) {
			if (rand_num(&f->random_seed)%10000 < prob*10000) {
				/* point mutation */
				f->actuator_destination[i] =
					rand_num(&f->random_seed)%
					f->no_of_actuator_destinations;
			}
		}
	}

	no_of_mutations =
		(int)(prob*0.75f*rows*columns)/(f->ADF_modules+1);
	for (m = 0; m < f->ADF_modules+1; m++) {

		/* swap chromosome positions */
		if (m == 0) {
			if ((rand_num(&f->random_seed)%10000) <
				chromosomes_prob*10000) {
				gprc_swap_chromosomes(f, m,
									  sensors, rows, columns,
									  connections_per_gene,
									  rand_num(&f->random_seed)%
									  chromosomes,
									  rand_num(&f->random_seed)%
									  chromosomes,
									  chromosomes,
									  (rand_num(&f->random_seed)%5)-2);
			}
		}

		/* mutate the genome */
		gene = f->genome[m].gene;
		act = gprc_get_actuators(m,actuators);
		sens = gprc_get_sensors(m,sensors);
		for (i = 0; i < no_of_mutations; i++) {
			/* pick a gene at random */
			index =
				rand_num(&f->random_seed)%(rows*columns*step + act);
			if (index < act) {
				/* mutate actuators */
				gprc_set_output_source(f, m, rows, columns,
									   connections_per_gene,
									   index,
									   sens + rand_num(&f->random_seed)%
									   (rows*columns));
			}
			else {
				index -= act;
				locn = index % step;
				/* first value */
				if (locn == GPRC_GENE_FUNCTION_TYPE) {
					/* function type */
					function_type =
						gpr_random_function(instruction_set,
											no_of_instructions,
											&f->random_seed);
					if (function_type != GPR_FUNCTION_ADF) {
						if (gene[index] == GPR_FUNCTION_ADF) {
							/* if this was previously an ADF
							   then change the first connection
							   back into the expected range */
							col =
								(index/
								 GPRC_GENE_SIZE(connections_per_gene)) /
								rows;
							new_connection =
								rand_num(&f->random_seed)%
								((col*rows) +
								 gprc_get_sensors(m,sensors));
							gene[index+GPRC_INITIAL] = new_connection;
						}
						gene[index] = function_type;
					}
					else {
						if ((m == 0) && (f->ADF_modules > 0)) {
							/* ADF module index */
							call_ADF_module =
								1 + (abs((int)gene[index+GPRC_GENE_CONSTANT])%
									 f->ADF_modules);
							if (gprc_contains_ADFs(f, 0,
												   call_ADF_module,
												   rows, columns,
												   connections_per_gene,
												   sensors) != -1) {
								gene[index] = function_type;
								gene[index+GPRC_INITIAL] =
									get_ADF_args(f, call_ADF_module)-1;
							}
						}
					}
				}
				/* second value */
				else if (locn < GPRC_INITIAL) {
					/* value */
					if (integers_only<=0) {
						if (rand_num(&f->random_seed)%10000>5000) {
							/* incremental */
							gene[index] =
								gpr_mutate_value(gene[index]+0.01f,
												 GPR_MUTATE_VALUE_PERCENT,
												 &f->random_seed);
						}
						else {
							/* random */
							gene[index] =
								gpr_random_value(min_value, max_value,
												 &f->random_seed);
						}
					}
					else {
						if (rand_num(&f->random_seed)%10000>5000) {
							/* incremental */
							if (rand_num(&f->random_seed)%2==0) {
								gene[index] = (int)gene[index] + 1;
							}
							else {
								gene[index] = (int)gene[index] - 1;
							}
						}
						else {
							/* random */
							gene[index] =
								(int)gpr_random_value(min_value,
													  max_value,
													  &f->random_seed);
						}
					}
				}
				/* connections */
				else {
					function_type = (int)gene[index-GPRC_INITIAL];
					if (!((function_type == GPR_FUNCTION_ADF) &&
						  (index%step == GPRC_INITIAL) && (m == 0))) {
						gene_index = index / step;
						conn_index = index % step;
						if (conn_index <
							GPRC_INITIAL+connections_per_gene) {
							/* connections */
							col = gene_index / rows;
							new_connection =
								rand_num(&f->random_seed)%
								((col*rows) +
								 gprc_get_sensors(m,sensors));
							gene[index] = new_connection;
						}
						else {
							/* weights */
							if (rand_num(&f->random_seed)%10000>5000) {
								/* incremental */
								gene[index] =
									gpr_mutate_value(gene[index]+0.01f,
													 GPR_MUTATE_VALUE_PERCENT,
													 &f->random_seed);
							}
							else {
								gene[index] =
									gpr_random_value(GPRC_MIN_WEIGHT,
													 GPRC_MIN_WEIGHT,
													 &f->random_seed);
							}
						}
					}
					else {
						/* ADF number of connections */
						call_ADF_module =
							1 + (abs((int)gene[index-1])%f->ADF_modules);
						gene[index] =
							get_ADF_args(f,call_ADF_module) - 1;
						/* if there are no arguments then make
						   this into a value function */
						if (gene[index] < 0) {
							gene[index-GPRC_INITIAL] =
								GPR_FUNCTION_VALUE;
							new_connection =
								rand_num(&f->random_seed)%
								((col*rows) +
								 gprc_get_sensors(m,sensors));
							gene[index] = new_connection;
						}
					}
				}
			}
		}
	}

	/* connection permutations */
	gprc_mutation_permute(f, rows, columns, connections_per_gene,
						  &f->random_seed,prob*0.25f);

	/* make sure thet output sources are unique */
	gprc_unique_outputs(f, rows, columns, connections_per_gene,
						sensors, actuators, &f->random_seed);

	/* ensure that any logical operators have valid inputs */
	gprc_valid_logical_operators(f, rows, columns,
								 connections_per_gene,
								 sensors, &f->random_seed);

	/* ensure that ADF calls are valid */
	gprc_valid_ADFs(f, rows, columns,
					connections_per_gene,
					sensors,
					min_value, max_value);
}

/* validate the genome */
int gprc_validate(gprc_function * f,
				  int rows, int columns,
				  int sensors, int actuators,
				  int connections_per_gene,
				  int integers_only,
				  int * instruction_set,
				  int no_of_instructions)
{
	int col,row,n,con,i,index,m,sens,act,ADF_args;
	int gene_args,call_ADF_module,min,max;
	int min_function_type=9999,max_function_type=0;
	int function_type;
	int previous_values;
	float * gene, * state;
	char name[256];

	/* get the minimum and maximum function type */
	for (i = 0; i < no_of_instructions; i++) {
		if (instruction_set[i]<min_function_type) {
			min_function_type = instruction_set[i];
		}
		if (instruction_set[i]>max_function_type) {
			max_function_type = instruction_set[i];
		}
	}

	/*for (m = 0; m < f->ADF_modules+1; m++) {*/
	for (m = 0; m < 1; m++) {
		sens = gprc_get_sensors(m,sensors);
		gene = f->genome[m].gene;
		state = f->genome[m].state;
		act = gprc_get_actuators(m,actuators);

		for (i = sens; i < sens + (rows*columns); i++) {
			if ((state[i] > GPR_MAX_CONSTANT) ||
				(state[i] < -GPR_MAX_CONSTANT)) {
				return GPR_VALIDATE_STATE_VALUE_OUT_OF_RANGE;
			}
		}

		/* check each gene */
		n=0;
		index=0;
		for (col = 0; col < columns; col++) {
			for (row = 0; row < rows; row++, index++,
					 n+=GPRC_GENE_SIZE(connections_per_gene)) {
				function_type = (int)gene[n];
				if (function_type==GPR_FUNCTION_ADF) {
					/* which ADF is being called */
					call_ADF_module =
						1 + (abs((int)gene[n+GPRC_GENE_CONSTANT])%
							 f->ADF_modules);
					/* how many arguments does the called ADF have? */
					ADF_args = get_ADF_args(f, call_ADF_module);
					/* how many arguments does the gene
					   say that it has? */
					gene_args =
						1 + (abs((int)gene[n+GPRC_INITIAL])%
							 GPRC_MAX_ADF_MODULE_SENSORS);
					if (gene_args != ADF_args) {
						printf("\nModule: %d  Gene: %d  " \
							   "Max_connections: %d\n", ADF_args,
							   gene_args, connections_per_gene);
						return GPR_VALIDATE_ADF_NO_OF_ARGS;
					}

					if (m > 0) {
						return GPR_VALIDATE_ADF_WITHIN_ADF;
					}
				}
				/* does the function belong in the set ? */
				if (gpr_function_in_set(function_type,
										instruction_set,
										no_of_instructions)==0) {
					printf("Function type: %d does not belong in set," \
						   " module %d\n",(int)gene[n],m);
					return GPR_VALIDATE_FUNCTION_TYPE_NOT_IN_SET;
				}
				/* check the function type */
				if ((int)gene[n] < min_function_type) {
					return GPR_VALIDATE_FUNCTION_LESS_THAN_MIN;
				}
				if (function_type > max_function_type) {
					printf("Function type: %d/%d\n",
						   (int)gene[n],max_function_type);
					return GPR_VALIDATE_FUNCTION_TYPE_TOO_LARGE;
				}

				/* get the number of connections */
				min = 0;
				max = gprc_function_args(function_type,
										 gene[1], connections_per_gene,
										 (int)gene[n+GPRC_INITIAL]);

				if (function_type==GPR_FUNCTION_ADF) {
					min = 1;
					max++;
				}

				/* check the connections */
				previous_values = (col*rows) + sens;
				for (con = min; con < max; con++) {
					if ((int)gene[n+GPRC_INITIAL+con] < 0) {
						printf("\nConnection: %d\n",con);
						return GPR_VALIDATE_CONNECTION_LESS_THAN_ZERO;
					}
					if ((int)gene[n+GPRC_INITIAL+con] >=
						previous_values) {
						gpr_get_function_name((int)gene[n],
											  1.0f, 0, name);

						printf("\nFunction: %s\n",name);
						printf("Connection index %d: %d/%d  %d %d\n",
							   con, (int)gene[n+GPRC_INITIAL+con],
							   previous_values,
							   max,connections_per_gene);
						return GPR_VALIDATE_CONNECTION_TOO_LARGE;
					}
				}
			}
		}

		/* check actuators */
		for (i = 0; i < act; i++, n++) {
			if ((int)gene[n] < 0) {
				return GPR_VALIDATE_ACTUATOR_CONNECTION_LESS_THAN_ZERO;
			}
			if ((int)gene[n] >= (col*rows) + sens) {
				return GPR_VALIDATE_ACTUATOR_CONNECTION_OUT_OF_RANGE;
			}
		}
	}

	return GPR_VALIDATE_OK;
}

/* runs an ADF */
static void gprc_c_run_ADF(gprc_function * f,
						   int ADF_module, int i,
						   float * gp,
						   int rows, int columns,
						   int connections_per_gene,
						   int sensors, int actuators,
						   float dropout_prob, int dynamic,
						   float (*custom_function)(float,float,float),
						   int integers_only)
{
	int call_ADF_module,itt,s,sens,argc=1,used_ctr;
	float * ADF_state, * state;
	unsigned char * ADF_used;

	if ((ADF_module != 0) || (f->ADF_modules == 0)) return;

	sens = gprc_get_sensors(ADF_module,sensors);
	state = f->genome[ADF_module].state;

	/* index of the ADF_module */
	call_ADF_module =
		1 + (abs((int)gp[GPRC_GENE_CONSTANT])%f->ADF_modules);

	/* get the number of arguments for the ADF */
	argc = 1 + ((abs((int)gp[GPRC_INITIAL]))%
				GPRC_MAX_ADF_MODULE_SENSORS);
	if (argc >= connections_per_gene) {
		argc = connections_per_gene-1;
	}

	/* clear the values of all ADF sensors to avoid
	   any residue from previous calls */
	memset((void*)f->genome[call_ADF_module].state,'\0',
		   GPRC_MAX_ADF_MODULE_SENSORS*sizeof(float));

	/* set the inputs to the ADF_module */
	ADF_used = f->genome[call_ADF_module].used;
	ADF_state = f->genome[call_ADF_module].state;
	used_ctr = 0;
	for (s = 0; s < argc; s++) {
		/* proceed to the next used ADF argument */
		while (ADF_used[used_ctr]==0) {
			used_ctr++;
		}
		if (integers_only < 1) {
			ADF_state[used_ctr] =
				state[(int)gp[1+GPRC_INITIAL+s]];
		}
		else {
			ADF_state[used_ctr] =
				(int)state[(int)gp[1+GPRC_INITIAL+s]];
		}
		used_ctr++;
	}

	/* run the ADF_module */
	for (itt = 0; itt < 2; itt++) {
		if (integers_only < 1) {
			gprc_run_float(f, call_ADF_module,
						   rows, columns,
						   connections_per_gene,
						   sensors, actuators,
						   dropout_prob, dynamic,
						   (*custom_function));
		}
		else {
			gprc_run_int(f, call_ADF_module,
						 rows, columns,
						 connections_per_gene,
						 sensors, actuators,
						 dropout_prob, dynamic,
						 (*custom_function));
		}
	}

	/* get the actuator value */
	state[sens+i] =
		gprc_get_ADF_module_actuator(f, 0, ADF_module,
									 rows, columns, sens);
	if (integers_only > 0) {
		state[sens+i] = (int)state[sens+i];
	}
}

/* run an individual */
void gprc_run_float(gprc_function * f,
					int ADF_module,
					int rows, int columns,
					int connections_per_gene,
					int sensors, int actuators,
					float dropout_prob,
					int dynamic,
					float (*custom_function)(float,float,float))
{
	int row,col,n=0,i=0,j,k,g,ctr,src,dest,no_of_args;
	float * gp;
	int gene_size = GPRC_GENE_SIZE(connections_per_gene);
	int block_from, block_to, act;
	int dropout = (int)(dropout_prob*10000);
	int sens = gprc_get_sensors(ADF_module,sensors);
	float * gene = f->genome[ADF_module].gene;
	unsigned char * used = f->genome[ADF_module].used;
	float * state = f->genome[ADF_module].state;

	act = gprc_get_actuators(ADF_module,actuators);

	for (col = 0; col < columns; col++) {
		for (row = 0; row < rows; row++, i++,n+=gene_size) {

			/* if this function is not on the path
			   between sensors and actuators then skip it
			   since it has no effect upon the program behavior */
			if ((dynamic <= 0) && (used[i+sens] == 0)) {
				continue;
			}

			/* occasional dropout helps to avoid overfitting*/
			if (rand_num(&f->random_seed)%10000<dropout) continue;

			gp = &gene[n];
			switch((int)gp[GPRC_GENE_FUNCTION_TYPE]) {
			case GPR_FUNCTION_GET: {				
				j = abs((int)state[(int)gp[GPRC_INITIAL]] +
						(int)state[(int)gp[1+GPRC_INITIAL]])
					%(rows*columns);
				state[sens+i] = state[sens+j];
				break;
			}
			case GPR_FUNCTION_SET: {
				j = abs((int)state[(int)gp[1+GPRC_INITIAL]])
					%(rows*columns);
				state[sens+i] = gp[GPRC_GENE_CONSTANT]*
					state[(int)gp[GPRC_INITIAL]];
				state[sens+j] = state[sens+i];
				if (state[sens+j] > GPR_MAX_CONSTANT) {
					state[sens+j] = GPR_MAX_CONSTANT;
				}
				if (state[sens+j] < -GPR_MAX_CONSTANT) {
					state[sens+j] = -GPR_MAX_CONSTANT;
				}
				break;
			}
			case GPR_FUNCTION_ADF: {
				gprc_c_run_ADF(f, ADF_module, i,
							   gp, rows, columns,
							   connections_per_gene,
							   sensors, actuators,
							   dropout_prob, dynamic,
							   (*custom_function),0);
				break;
			}
			case GPR_FUNCTION_CUSTOM: {
				if (*custom_function) {
					state[sens+i] =
						(*custom_function)(gp[GPRC_GENE_CONSTANT],
										   gp[GPRC_INITIAL],
										   gp[GPRC_GENE_CONSTANT]);
				}
				break;
			}
			case GPR_FUNCTION_VALUE: {
				state[sens+i] = gp[GPRC_GENE_CONSTANT];
				break;
			}
			case GPR_FUNCTION_SIGMOID: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = 0;
				for (j = 0; j < no_of_args; j++) {
					state[sens+i] +=
						state[(int)gp[GPRC_INITIAL+j]]*
						gp[GPRC_INITIAL+j+connections_per_gene];
				}

				state[sens+i] =
					1.0f / (1.0f + exp(-state[sens+i]));
				break;
			}
			case GPR_FUNCTION_ADD: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = 0;
				for (j = 0; j < no_of_args; j++) {
					state[sens+i] += state[(int)gp[GPRC_INITIAL+j]];
				}
				break;
			}
			case GPR_FUNCTION_SUBTRACT: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					state[sens+i] -= state[(int)gp[GPRC_INITIAL+j]];
				}
				break;
			}
			case GPR_FUNCTION_NEGATE: {
				state[sens+i] = -state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_MULTIPLY: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					state[sens+i] *= state[(int)gp[GPRC_INITIAL+j]];
				}
				break;
			}
			case GPR_FUNCTION_WEIGHT: {
				state[sens+i] = state[(int)gp[GPRC_INITIAL]] *
					gp[GPRC_GENE_CONSTANT];
				break;
			}
			case GPR_FUNCTION_DIVIDE: {
				if((state[(int)gp[1+GPRC_INITIAL]] <= 1e-1) &&
				   (state[(int)gp[1+GPRC_INITIAL]] >= -1e-1)) {
					state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				}
				else {
					state[sens+i] =
						state[(int)gp[GPRC_INITIAL]] /
						state[(int)gp[1+GPRC_INITIAL]];
				}
				break;
			}
			case GPR_FUNCTION_MODULUS: {
				if((int)state[(int)gp[1+GPRC_INITIAL]] == 0) {
					state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				}
				else {
					state[sens+i] =
						(int)state[(int)gp[GPRC_INITIAL]]%
						(int)state[(int)gp[1+GPRC_INITIAL]];
				}
				break;
			}
			case GPR_FUNCTION_FLOOR: {
				state[sens+i] = floor(state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_AVERAGE: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					state[sens+i] += state[(int)gp[GPRC_INITIAL+j]];
				}
				state[sens+i] /= no_of_args;
				break;
			}
			case GPR_FUNCTION_NOOP1: {
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_NOOP2: {
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_NOOP3: {
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_NOOP4: {
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_GREATER_THAN: {
				if (state[(int)gp[GPRC_INITIAL]] >
					state[(int)gp[1+GPRC_INITIAL]]) {
					state[sens+i] = gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_LESS_THAN: {
				if (state[(int)gp[GPRC_INITIAL]] <
					state[(int)gp[1+GPRC_INITIAL]]) {
					state[sens+i] = gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_EQUALS: {
				if ((int)state[(int)gp[GPRC_INITIAL]] ==
					(int)state[(int)gp[1+GPRC_INITIAL]]) {
					state[sens+i] = gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_AND: {
				if ((state[(int)gp[GPRC_INITIAL]]>0) &&
					(state[(int)gp[1+GPRC_INITIAL]]>0)) {
					state[sens+i] = gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_OR: {
				if ((state[(int)gp[GPRC_INITIAL]]>0) ||
					(state[(int)gp[1+GPRC_INITIAL]]>0)) {
					state[sens+i] = gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_XOR: {
				if ((state[(int)gp[GPRC_INITIAL]]>0) !=
					(state[(int)gp[1+GPRC_INITIAL]]>0)) {
					state[sens+i] = gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_NOT: {
				if (((int)state[(int)gp[GPRC_INITIAL]]) !=
					((int)state[(int)gp[1+GPRC_INITIAL]])) {
					state[sens+i] = gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_HEBBIAN: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				/* update the output */
				state[sens+i] = 0;
				for (j = 0; j < no_of_args; j++) {
					state[sens+i] +=
						state[(int)gp[GPRC_INITIAL+j]] *
						gp[GPRC_INITIAL+j+connections_per_gene];
				}
				/* adjust weights */
				for (j = 0; j < no_of_args; j++) {
					gp[GPRC_INITIAL+j+connections_per_gene] +=
						state[sens+i] * state[(int)gp[GPRC_INITIAL+j]] *
						GPR_HEBBIAN_LEARNING_RATE;
				}
				break;
			}
			case GPR_FUNCTION_EXP: {
				state[sens+i] = (float)exp(state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_SQUARE_ROOT: {
				state[sens+i] =
					(float)sqrt(fabs(state[(int)gp[GPRC_INITIAL]]));
				break;
			}
			case GPR_FUNCTION_ABS: {
				state[sens+i] =
					(float)fabs(state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_SINE: {
				state[sens+i] =
					(float)sin(state[(int)gp[GPRC_INITIAL]])*256;
				break;
			}
			case GPR_FUNCTION_ARCSINE: {
				state[sens+i] =
					(float)asin(state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_COSINE: {
				state[sens+i] =
					(float)cos(state[(int)gp[GPRC_INITIAL]])*256;
				break;
			}
			case GPR_FUNCTION_ARCCOSINE: {
				state[sens+i] =
					(float)acos(state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_POW: {
				state[sens+i] =
					(float)pow(state[(int)gp[GPRC_INITIAL]],
							   state[(int)gp[1+GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_MIN: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					if (state[(int)gp[GPRC_INITIAL+j]] < state[sens+i]) {
						state[sens+i] = state[(int)gp[GPRC_INITIAL+j]];
					}
				}
				break;
			}
			case GPR_FUNCTION_MAX: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					if (state[(int)gp[GPRC_INITIAL+j]] > state[sens+i]) {
						state[sens+i] = state[(int)gp[GPRC_INITIAL+j]];
					}
				}
				break;
			}
			case GPR_FUNCTION_COPY_FUNCTION: {
				if (((int)gp[GPRC_INITIAL] > sens) &&
					((int)gp[1+GPRC_INITIAL] > sens)) {
					src = ((int)gp[GPRC_INITIAL]-sens) * gene_size;
					dest = ((int)gp[1+GPRC_INITIAL]-sens) * gene_size;
					gene[dest] = gene[src];					
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONSTANT: {
				if (((int)gp[GPRC_INITIAL] > sens) &&
					((int)gp[1+GPRC_INITIAL] > sens)) {
					src = ((int)gp[GPRC_INITIAL]-sens) * gene_size;
					dest = ((int)gp[1+GPRC_INITIAL]-sens) * gene_size;
					gene[dest+GPRC_GENE_CONSTANT] =
						gene[src+GPRC_GENE_CONSTANT];
				}
				break;
			}
			case GPR_FUNCTION_COPY_STATE: {
				state[(int)gp[1+GPRC_INITIAL]] =
					state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_COPY_BLOCK: {
				block_from = (int)gp[GPRC_INITIAL];
				block_to = (int)gp[1+GPRC_INITIAL];
				if (block_from < block_to) {
					block_from = (int)gp[1+GPRC_INITIAL];
					block_to = (int)gp[GPRC_INITIAL];
				}
				k = block_to - GPR_BLOCK_WIDTH;
				for (j = block_from - GPR_BLOCK_WIDTH;
					 j <= block_from + GPR_BLOCK_WIDTH; j++,k++) {
					if ((j>sens) &&
						(k>sens) &&
						(j<i) && (k<i)) {
						for (g = 0; g < gene_size; g++) {
							gene[(j-sens)*gene_size + g] =
								gene[(k-sens)*gene_size + g];
						}
					}
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION1: {
				if (gp[GPRC_INITIAL] > sens) {
					src = ((int)gp[GPRC_INITIAL] - sens) * gene_size;
					gp[1+GPRC_INITIAL] = gene[src+GPRC_INITIAL];
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION2: {
				if (gp[1+GPRC_INITIAL] > sens) {
					src = ((int)gp[1+GPRC_INITIAL] - sens) * gene_size;
					gp[GPRC_INITIAL] = gene[src+GPRC_INITIAL];
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION3: {
				if (gp[1+GPRC_INITIAL] > sens) {
					src = ((int)gp[1+GPRC_INITIAL] - sens) * gene_size;
					gp[GPRC_INITIAL] = gene[src+1+GPRC_INITIAL];
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION4: {
				if (gp[GPRC_INITIAL] > sens) {
					src = ((int)gp[GPRC_INITIAL] - sens) * gene_size;
					gp[1+GPRC_INITIAL] = gene[src+1+GPRC_INITIAL];
				}
				break;
			}
			}
			/* prevent values from going out of range */
			if (is_nan(state[sens+i])) {
				state[sens+i] = 0;
			}
			if (state[sens+i] > GPR_MAX_CONSTANT) {
				state[sens+i] = GPR_MAX_CONSTANT;
			}
			if (state[sens+i] < -GPR_MAX_CONSTANT) {
				state[sens+i] = -GPR_MAX_CONSTANT;
			}
		}
	}

	/* set the actuator values */
	ctr = sens + i;
	for (i = 0; i < act; i++, ctr++, n++) {
		state[ctr] = state[(int)gene[n]];
	}
}

/* an integer version of the run function */
void gprc_run_int(gprc_function * f,
				  int ADF_module,
				  int rows, int columns,
				  int connections_per_gene,
				  int sensors, int actuators,
				  float dropout_prob, int dynamic,
				  float (*custom_function)(float,float,float))
{
	int row,col,n=0,i=0,j,k,g,ctr,src,dest,no_of_args;
	float * gp;
	int gene_size = GPRC_GENE_SIZE(connections_per_gene);
	int block_from,block_to;
	int dropout = (int)(dropout_prob*10000);
	int sens = gprc_get_sensors(ADF_module,sensors);
	float * gene = f->genome[ADF_module].gene;
	unsigned char * used = f->genome[ADF_module].used;
	float * state = f->genome[ADF_module].state;

	actuators = gprc_get_actuators(ADF_module,actuators);

	for (col = 0; col < columns; col++) {
		for (row = 0; row < rows; row++, i++,n+=gene_size) {

			/* if this function is not on the path
			   between sensors and actuators then skip it
			   since it has no effect upon the program behavior */
			if ((dynamic <= 0) && (used[i+sens] == 0)) {
				continue;
			}

			/* occasional dropout helps to avoid overfitting*/
			if (rand_num(&f->random_seed)%10000<dropout) continue;

			gp = &gene[n];
			switch((int)gp[GPRC_GENE_FUNCTION_TYPE]) {
			case GPR_FUNCTION_GET: {				
				j = abs((int)state[(int)gp[GPRC_INITIAL]] +
						(int)state[(int)gp[1+GPRC_INITIAL]])
					%(rows*columns);
				state[sens+i] = (int)state[sens+j];
				break;
			}
			case GPR_FUNCTION_SET: {
				j = abs((int)state[(int)gp[1+GPRC_INITIAL]])
					%(rows*columns);
				state[sens+i] =
					(int)gp[GPRC_GENE_CONSTANT]*
					(int)state[(int)gp[GPRC_INITIAL]];
				state[sens+j] = (int)state[sens+i];
				if (state[sens+j] > GPR_MAX_CONSTANT) {
					state[sens+j] = GPR_MAX_CONSTANT;
				}
				if (state[sens+j] < -GPR_MAX_CONSTANT) {
					state[sens+j] = -GPR_MAX_CONSTANT;
				}
				break;
			}
			case GPR_FUNCTION_ADF: {
				gprc_c_run_ADF(f, ADF_module, i,
							   gp, rows, columns,
							   connections_per_gene,
							   sensors, actuators,
							   dropout_prob, dynamic,
							   (*custom_function),1);
				break;
			}
			case GPR_FUNCTION_CUSTOM: {
				if (*custom_function) {
					state[sens+i] =
						(*custom_function)((int)gp[GPRC_GENE_CONSTANT],
										   (int)gp[GPRC_INITIAL],
										   (int)gp[GPRC_GENE_CONSTANT]);
				}
				break;
			}
			case GPR_FUNCTION_VALUE: {
				state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				break;
			}
			case GPR_FUNCTION_SIGMOID: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = 0;
				for (j = 0; j < no_of_args; j++) {
					state[sens+i] +=
						state[(int)gp[GPRC_INITIAL+j]]*
						gp[GPRC_INITIAL+j+connections_per_gene];
				}

				state[sens+i] =
					1.0f / (1.0f + exp(-state[sens+i]));
				break;
			}
			case GPR_FUNCTION_ADD: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = 0;
				for (j = 0; j < no_of_args; j++) {
					state[sens+i] += (int)state[(int)gp[GPRC_INITIAL+j]];
				}
				break;
			}
			case GPR_FUNCTION_SUBTRACT: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					state[sens+i] -= (int)state[(int)gp[GPRC_INITIAL+j]];
				}
				break;
			}
			case GPR_FUNCTION_NEGATE: {
				state[sens+i] = -(int)state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_MULTIPLY: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					state[sens+i] *= (int)state[(int)gp[GPRC_INITIAL+j]];
				}
				break;
			}
			case GPR_FUNCTION_WEIGHT: {
				state[sens+i] = state[(int)gp[GPRC_INITIAL]] *
					(int)gp[GPRC_GENE_CONSTANT];
				break;
			}
			case GPR_FUNCTION_DIVIDE: {
				if ((int)state[(int)gp[1+GPRC_INITIAL]] == 0) {
					state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				}
				else {
					state[sens+i] =
						(int)state[(int)gp[GPRC_INITIAL]] /
						(int)state[(int)gp[1+GPRC_INITIAL]];
				}
				break;
			}
			case GPR_FUNCTION_MODULUS: {
				if ((int)state[(int)gp[1+GPRC_INITIAL]] == 0) {
					state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				}
				else {
					state[sens+i] =
						(int)state[(int)gp[GPRC_INITIAL]] %
						(int)state[(int)gp[1+GPRC_INITIAL]];
				}
				break;
			}
			case GPR_FUNCTION_FLOOR: {
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_AVERAGE: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = 0;
				for (j = 0; j < no_of_args; j++) {
					state[sens+i] += (int)state[(int)gp[GPRC_INITIAL+j]];
				}
				state[sens+i] /= no_of_args;
				break;
			}
			case GPR_FUNCTION_NOOP1: {
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_NOOP2: {
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_NOOP3: {
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_NOOP4: {
				state[sensors+i] = (int)state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_GREATER_THAN: {
				if ((int)state[(int)gp[GPRC_INITIAL]] >
					(int)state[(int)gp[1+GPRC_INITIAL]]) {
					state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_LESS_THAN: {
				if ((int)state[(int)gp[GPRC_INITIAL]] <
					(int)state[(int)gp[1+GPRC_INITIAL]]) {
					state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_EQUALS: {
				if ((int)state[(int)gp[GPRC_INITIAL]] ==
					(int)state[(int)gp[1+GPRC_INITIAL]]) {
					state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_AND: {
				if (((int)state[(int)gp[GPRC_INITIAL]]>0) &&
					((int)state[(int)gp[1+GPRC_INITIAL]]>0)) {
					state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_OR: {
				if (((int)state[(int)gp[GPRC_INITIAL]]>0) ||
					((int)state[(int)gp[1+GPRC_INITIAL]]>0)) {
					state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_XOR: {
				if (((int)state[(int)gp[GPRC_INITIAL]]>0) !=
					((int)state[(int)gp[1+GPRC_INITIAL]]>0)) {
					state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_NOT: {
				if (((int)state[(int)gp[GPRC_INITIAL]]) !=
					((int)state[(int)gp[1+GPRC_INITIAL]])) {
					state[sens+i] = (int)gp[GPRC_GENE_CONSTANT];
				}
				else {
					state[sens+i] = 0;
				}
				break;
			}
			case GPR_FUNCTION_HEBBIAN: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				/* update the output */
				state[sens+i] = 0;
				for (j = 0; j < no_of_args; j++) {
					state[sens+i] +=
						state[(int)gp[GPRC_INITIAL+j]] *
						gp[GPRC_INITIAL+j+connections_per_gene];
				}
				/* adjust weights */
				for (j = 0; j < no_of_args; j++) {
					gp[GPRC_INITIAL+j+connections_per_gene] +=
						state[sens+i] * state[(int)gp[GPRC_INITIAL+j]] *
						GPR_HEBBIAN_LEARNING_RATE;
				}
				break;
			}
			case GPR_FUNCTION_EXP: {
				state[sens+i] =
					(int)exp((int)state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_SQUARE_ROOT: {
				state[sens+i] =
					(int)sqrt((int)fabs(state[(int)gp[GPRC_INITIAL]]));
				break;
			}
			case GPR_FUNCTION_ABS: {
				state[sens+i] =
					(int)abs((int)state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_SINE: {
				state[sens+i] =
					(int)(sin((int)state[(int)gp[GPRC_INITIAL]])*256);
				break;
			}
			case GPR_FUNCTION_ARCSINE: {
				state[sens+i] =
					(int)asin((int)state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_COSINE: {
				state[sens+i] =
					(int)(cos((int)state[(int)gp[GPRC_INITIAL]])*256);
				break;
			}
			case GPR_FUNCTION_ARCCOSINE: {
				state[sens+i] =
					(int)acos((int)state[(int)gp[GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_POW: {
				state[sens+i] =
					(int)pow((int)state[(int)gp[GPRC_INITIAL]],
							 (int)state[(int)gp[1+GPRC_INITIAL]]);
				break;
			}
			case GPR_FUNCTION_MIN: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					if ((int)state[(int)gp[GPRC_INITIAL+j]] <
						state[sens+i]) {
						state[sens+i] =
							(int)state[(int)gp[GPRC_INITIAL+j]];
					}
				}
				break;
			}
			case GPR_FUNCTION_MAX: {
				no_of_args =
					1 + (abs((int)gp[GPRC_GENE_CONSTANT])%
						 (connections_per_gene-1));
				state[sens+i] = (int)state[(int)gp[GPRC_INITIAL]];
				for (j = 1; j < no_of_args; j++) {
					if ((int)state[(int)gp[GPRC_INITIAL+j]] >
						state[sens+i]) {
						state[sens+i] =
							(int)state[(int)gp[GPRC_INITIAL+j]];
					}
				}
				break;
			}
			case GPR_FUNCTION_COPY_FUNCTION: {
				if (((int)gp[GPRC_INITIAL] > sens) &&
					((int)gp[1+GPRC_INITIAL] > sens)) {
					src = ((int)gp[GPRC_INITIAL]-sens) * gene_size;
					dest = ((int)gp[1+GPRC_INITIAL]-sens) * gene_size;
					gene[dest] = gene[src];
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONSTANT: {
				if (((int)gp[GPRC_INITIAL] > sens) &&
					((int)gp[1+GPRC_INITIAL] > sens)) {
					src = ((int)gp[GPRC_INITIAL]-sens) * gene_size;
					dest = ((int)gp[1+GPRC_INITIAL]-sens) * gene_size;
					gene[dest+GPRC_GENE_CONSTANT] =
						gene[src+GPRC_GENE_CONSTANT];
				}
				break;
			}
			case GPR_FUNCTION_COPY_STATE: {
				state[(int)gp[1+GPRC_INITIAL]] =
					state[(int)gp[GPRC_INITIAL]];
				break;
			}
			case GPR_FUNCTION_COPY_BLOCK: {
				block_from = (int)gp[GPRC_INITIAL];
				block_to = (int)gp[1+GPRC_INITIAL];
				if (block_from<block_to) {
					block_from = (int)gp[1+GPRC_INITIAL];
					block_to = (int)gp[GPRC_INITIAL];
				}
				k = block_to - GPR_BLOCK_WIDTH;
				for (j = block_from - GPR_BLOCK_WIDTH;
					 j <= block_from + GPR_BLOCK_WIDTH; j++,k++) {
					if ((j>sens) &&
						(k>sens) &&
						(j<i) && (k<i)) {

						for (g = 0; g < gene_size; g++) {
							gene[(j-sens)*gene_size + g] =
								gene[(k-sens)*gene_size + g];
						}
					}
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION1: {
				if (gp[GPRC_INITIAL] > sens) {
					src = ((int)gp[GPRC_INITIAL] - sens) * gene_size;
					gp[1+GPRC_INITIAL] = gene[src+GPRC_INITIAL];
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION2: {
				if (gp[1+GPRC_INITIAL] > sens) {
					src = ((int)gp[1+GPRC_INITIAL] - sens) * gene_size;
					gp[GPRC_INITIAL] = gene[src+GPRC_INITIAL];
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION3: {
				if (gp[1+GPRC_INITIAL] > sens) {
					src = ((int)gp[1+GPRC_INITIAL] - sens) * gene_size;
					gp[GPRC_INITIAL] = gene[src+1+GPRC_INITIAL];
				}
				break;
			}
			case GPR_FUNCTION_COPY_CONNECTION4: {
				if (gp[GPRC_INITIAL] > sens) {
					src = ((int)gp[GPRC_INITIAL] - sens) * gene_size;
					gp[1+GPRC_INITIAL] = gene[src+1+GPRC_INITIAL];
				}
				break;
			}
			}
			/* prevent values from going out of range */
			if (is_nan(state[sens+i])) {
				state[sens+i] = 0;
			}
			if (state[sens+i] > GPR_MAX_CONSTANT) {
				state[sens+i] = GPR_MAX_CONSTANT;
			}
			if (state[sens+i] < -GPR_MAX_CONSTANT) {
				state[sens+i] = -GPR_MAX_CONSTANT;
			}
		}
	}

	/* set the actuator values */
	ctr = sens + i;
	for (i = 0; i < actuators; i++, ctr++, n++) {
		state[ctr] = (int)state[(int)gene[n]];
	}
}

void gprc_run(gprc_function * f, gprc_population * population,
			  float dropout_prob, int dynamic,
			  float (*custom_function)(float,float,float))
{
	if (population->integers_only<=0) {
		gprc_run_float(f, 0,
					   population->rows, population->columns,
					   population->connections_per_gene,
					   population->sensors, population->actuators,
					   dropout_prob, dynamic, (*custom_function));
	}
	else {
		gprc_run_int(f, 0,
					 population->rows, population->columns,
					 population->connections_per_gene,
					 population->sensors, population->actuators,
					 dropout_prob, dynamic, (*custom_function));
	}
}

void gprc_run_environment(gprc_function * f, gprc_environment * population,
						  float dropout_prob, int dynamic,
						  float (*custom_function)(float,float,float))
{
	if (population->integers_only<=0) {
		gprc_run_float(f, 0,
					   population->rows, population->columns,
					   population->connections_per_gene,
					   population->sensors, population->actuators,
					   dropout_prob, dynamic, (*custom_function));
	}
	else {
		gprc_run_int(f, 0,
					 population->rows, population->columns,
					 population->connections_per_gene,
					 population->sensors, population->actuators,
					 dropout_prob, dynamic, (*custom_function));
	}
}

/* initialize the population */
void gprc_init_population(gprc_population * population,
						  int size,
						  int rows, int columns,
						  int sensors, int actuators,
						  int connections_per_gene,
						  int ADF_modules,
						  int chromosomes,
						  float min_value, float max_value,
						  int integers_only, unsigned int * random_seed,
						  int * instruction_set, int no_of_instructions)
{
	int i;

	population->individual =
		(gprc_function*)malloc(size*sizeof(gprc_function));
	population->size = size;
	population->rows = rows;
	population->columns = columns;
	population->sensors = sensors;
	population->actuators = actuators;
	population->connections_per_gene = connections_per_gene;
	population->chromosomes = chromosomes;
	population->ADF_modules = ADF_modules;
	population->min_value = min_value;
	population->max_value = max_value;
	population->integers_only = integers_only;
	population->fitness = (float*)malloc(size*sizeof(float));

	population->history.index = 0;
	population->history.interval = 1;
	population->history.tick = 0;

	for (i = 0; i < size; i++) {
		/* initialise the individual */
		gprc_init(&population->individual[i],
				  rows, columns, sensors, actuators,
				  connections_per_gene, ADF_modules, random_seed);

		/* initialise individuals randomly */
		gprc_random(&population->individual[i],
					rows, columns,
					sensors, actuators,
					connections_per_gene,
					min_value, max_value,
					integers_only, random_seed,
					instruction_set, no_of_instructions);

		/* remove any ADF functions */
		gprc_remove_ADFs(&population->individual[i],
						 rows, columns,
						 connections_per_gene);

		/* ensure that any ADFs are valid */
		gprc_valid_ADFs(&population->individual[i],
						rows, columns,
						connections_per_gene,
						sensors,
						min_value, max_value);

		/* clear the fitness value */
		population->fitness[i] = 0;

		/* update the used functions */
		gprc_used_functions(&population->individual[i],
							rows, columns,
							connections_per_gene,
							sensors, actuators);
	}
}

/* initialize the population of an environment */
void gprc_init_environment(gprc_environment * population,
						   int max_population_size,
						   int initial_population_size,
						   int rows, int columns,
						   int sensors, int actuators,
						   int connections_per_gene,
						   int ADF_modules,
						   int chromosomes,
						   float min_value, float max_value,
						   int integers_only,
						   unsigned int * random_seed,
						   int * instruction_set,
						   int no_of_instructions)
{
	int i;

	population->individual =
		(gprc_function*)malloc(max_population_size*
							   sizeof(gprc_function));
	population->max_population_size = max_population_size;
	population->population_size = initial_population_size;
	population->rows = rows;
	population->columns = columns;
	population->sensors = sensors;
	population->actuators = actuators;
	population->connections_per_gene = connections_per_gene;
	population->chromosomes = chromosomes;
	population->ADF_modules = ADF_modules;
	population->min_value = min_value;
	population->max_value = max_value;
	population->integers_only = integers_only;
	population->mating =
		(int*)malloc(max_population_size*3*sizeof(int));
	population->matings = 0;

	for (i = 0; i < max_population_size; i++) {
		/* initialise the individual */
		gprc_init(&population->individual[i],
				  rows, columns, sensors, actuators,
				  connections_per_gene, ADF_modules, random_seed);

		/* initialise individuals randomly */
		gprc_random(&population->individual[i],
					rows, columns,
					sensors, actuators,
					connections_per_gene,
					min_value, max_value,
					integers_only, random_seed,
					instruction_set, no_of_instructions);

		/* remove any ADF functions */
		gprc_remove_ADFs(&population->individual[i],
						 rows, columns,
						 connections_per_gene);

		/* ensure that any ADFs are valid */
		gprc_valid_ADFs(&population->individual[i],
						rows, columns,
						connections_per_gene,
						sensors,
						min_value, max_value);

		/* update the used functions */
		gprc_used_functions(&population->individual[i],
							rows, columns,
							connections_per_gene,
							sensors, actuators);
	}
}

/* initialise a system which contains multiple sub-populations */
void gprc_init_system(gprc_system * system,
					  int islands,
					  int population_per_island,
					  int rows, int columns,
					  int sensors, int actuators,
					  int connections_per_gene,
					  int ADF_modules,
					  int chromosomes,
					  float min_value, float max_value,
					  int integers_only, unsigned int * random_seed,
					  int * instruction_set, int no_of_instructions)
{
	int i;

	system->size = islands;
	system->migration_tick=0;
	system->island =
		(gprc_population*)malloc(islands*sizeof(gprc_population));
	system->fitness = (float*)malloc(islands*sizeof(float));

	/* clear the fitness values */
	for (i = 0; i < islands; i++) {

		/* create a population for the island */
		gprc_init_population(&system->island[i],
							 population_per_island,
							 rows, columns,
							 sensors, actuators,
							 connections_per_gene,
							 ADF_modules,
							 chromosomes,
							 min_value, max_value,
							 integers_only, random_seed,
							 instruction_set, no_of_instructions);

		/* clear the average fitness for the population */
		system->fitness[i] = 0;
	}
}

/* frees memory for a system */
void gprc_free_system(gprc_system * system)
{
	for (int i = 0; i < system->size; i++) {
		gprc_free_population(&system->island[i]);
	}
	free(system->island);
	free(system->fitness);
}

/* deallocates memory for the given population */
void gprc_free_population(gprc_population * population)
{
	for (int i = 0; i < population->size; i++) {
		gprc_free(&population->individual[i]);
	}
	free(population->individual);
	free(population->fitness);
}

/* deallocates memory for the given environment */
void gprc_free_environment(gprc_environment * population)
{
	for (int i = 0; i < population->max_population_size; i++) {
		gprc_free(&population->individual[i]);
	}
	free(population->individual);
	free(population->mating);
}

/* Evaluates the fitness of all individuals in the population.
   Here we use openmp to speed up the process, since each
   evaluation is independent */
void gprc_evaluate(gprc_population * population,
				   int time_steps, int reevaluate,
				   float (*evaluate_program)
				   (int,gprc_population*,int,int))
{
	int i;

#pragma omp parallel for
	for (i = 0; i < population->size; i++) {
		if ((population->fitness[i]==0) ||
			(reevaluate>0)) {
			int s;
			gprc_function * f = &population->individual[i];
			unsigned char * used = f->genome[0].used;			
			/* clear the retained state */
			gprc_clear_state(f,
							 population->rows, population->columns,
							 population->sensors,
							 population->actuators);

			/* is there a path which links sensors to actuators? */
			for (s = 0; s < population->sensors; s++) {
				if (used[s] != 0) break;
			}
			
			if (s < population->sensors) {
				/* run the evaluation function */
				population->fitness[i] =
					(*evaluate_program)(time_steps,population,i,0);
			}
			else {
				/* don't evaluate, since there is no path between
				   sensors and actuators */
				population->fitness[i] = 0;
			}
		}
		/* if individual gets too old */
		(&population->individual[i])->age++;
		if ((&population->individual[i])->age>GPR_MAX_AGE) {
			population->fitness[i] = 0;
		}
	}
}

/* evaluates a system containing multiple sub-populations */
void gprc_evaluate_system(gprc_system * system,
						  int time_steps, int reevaluate,
						  float (*evaluate_program)
						  (int,gprc_population*,int,int))
{
	int i;

#pragma omp parallel for
	for (i = 0; i < system->size; i++) {
		/* evaluate the island population */
		gprc_evaluate(&system->island[i],
					  time_steps, reevaluate,
					  (*evaluate_program));
		/* set the average fitness */
		system->fitness[i] = gprc_average_fitness(&system->island[i]);
	}
}

/* returns the highest fitness value */
float gprc_best_fitness(gprc_population * population)
{
	return population->fitness[0];
}

/* returns the median fitness value */
float gprc_median_fitness(gprc_population * population)
{
	return population->fitness[population->size/2];
}

/* returns the highest fitness value for the given system */
float gprc_best_fitness_system(gprc_system * system)
{
	return gprc_best_fitness(&system->island[0]);
}

/* returns the lowest fitness value */
float gprc_worst_fitness(gprc_population * population)
{
	return population->fitness[population->size-1];
}

/* returns the average fitness of the population */
float gprc_average_fitness(gprc_population * population)
{
	int i;
	float av = 0;

	for (i = 0; i < population->size; i++) {
		av += population->fitness[i];
	}
	return av / population->size;
}

/* returns the fittest individual in the population */
gprc_function * gprc_best_individual(gprc_population * population)
{
	return &population->individual[0];
}

/* returns the fittest individual in the given system */
gprc_function * gprc_best_individual_system(gprc_system * system)
{
	return gprc_best_individual(&system->island[0]);
}

/* set a sensor to the given value */
void gprc_set_sensor(gprc_function * f, int index, float value)
{
	float * state = f->genome[0].state;
	state[index] = value;
}

/* returns the value of a sensor */
float gprc_get_sensor(gprc_function * f, int index)
{
	float * state = f->genome[0].state;
	return state[index];
}

/* returns the sensor source identifier for the given sensor */
int gprc_get_sensor_source(gprc_function * f, int index)
{
	return f->sensor_source[index];
}

/* returns the value of an actuator */
float gprc_get_actuator(gprc_function * f, int index,
						int rows, int columns, int sensors)
{
	return gprc_get_ADF_module_actuator(f,index,0,rows,columns,sensors);
}

/* returns the actuator destination identifier
   for the given actuator */
int gprc_get_actuator_destination(gprc_function * f, int index)
{
	return f->actuator_destination[index];
}

/* sorts individuals in order of fitness */
void gprc_sort(gprc_population * population)
{
	int i, j, best;
	float max, temp_fitness;
	gprc_function temp_individual;

	for (i = 0; i < population->size-1; i++) {
		max = population->fitness[i];
		best = i;
		for (j = i+1; j < population->size; j++) {
			if (population->fitness[j] > max) {
				max = population->fitness[j];
				best = j;
			}
		}
		if (best != i) {
			/* swap fitness */
			temp_fitness = population->fitness[i];
			population->fitness[i] = population->fitness[best];
			population->fitness[best] = temp_fitness;

			/* swap individual */
			temp_individual = population->individual[i];
			population->individual[i] = population->individual[best];
			population->individual[best] = temp_individual;
		}
	}
}

/* sorts populations in order of average fitness */
void gprc_sort_system(gprc_system * system)
{
	int i, j, best;
	float max, temp_fitness;
	gprc_population temp_island;

	for (i = 0; i < system->size-1; i++) {
		max = system->fitness[i];
		best = i;
		for (j = i+1; j < system->size; j++) {
			if (system->fitness[j] > max) {
				max = system->fitness[j];
				best = j;
			}
		}
		if (best != i) {
			/* swap fitness */
			temp_fitness = system->fitness[i];
			system->fitness[i] = system->fitness[best];
			system->fitness[best] = temp_fitness;

			/* swap island */
			temp_island = system->island[i];
			system->island[i] = system->island[best];
			system->island[best] = temp_island;
		}
	}
}

/* copy a ADF_module from one individual to another */
void gprc_copy_ADF_module(gprc_function * source, gprc_function * dest,
						  int ADF_module,
						  int rows, int columns,
						  int connections_per_gene,
						  int sensors, int actuators)
{
	float * source_gene, * dest_gene;
	unsigned char * source_used, * dest_used;

	source_gene = source->genome[ADF_module].gene;
	dest_gene = dest->genome[ADF_module].gene;
	source_used = source->genome[ADF_module].used;
	dest_used = dest->genome[ADF_module].used;

	/* copy the genome */
	memcpy((void*)dest_gene,(void*)source_gene,
		   ((rows*columns*GPRC_GENE_SIZE(connections_per_gene)) +
			gprc_get_actuators(ADF_module,actuators))*sizeof(float));

	/* copy the used functions */
	memcpy((void*)dest_used,(void*)source_used,
		   ((rows*columns) +
			gprc_get_sensors(ADF_module,sensors) +
			gprc_get_actuators(ADF_module,actuators))*
		   sizeof(unsigned char));
}

/* copies the source genome to the destination genome */
void gprc_copy(gprc_function * source, gprc_function * dest,
			   int rows, int columns, int connections_per_gene,
			   int sensors, int actuators)
{
	int m, min_ADF_modules = source->ADF_modules;

	if (dest->ADF_modules < min_ADF_modules) {
		min_ADF_modules = dest->ADF_modules;
	}

	/* copy the sensor sources */
	if (source->no_of_sensor_sources>0) {
		if (dest->sensor_source==0) {
			dest->no_of_sensor_sources = source->no_of_sensor_sources;
			dest->sensor_source = (int*)malloc(sensors*sizeof(int));
		}
		memcpy((void*)dest->sensor_source,
			   (void*)source->sensor_source,
			   sensors*sizeof(int));
	}

	/* copy the actuator destinations */
	if (source->no_of_actuator_destinations>0) {
		if (dest->actuator_destination==0) {
			dest->no_of_actuator_destinations =
				source->no_of_actuator_destinations;
			dest->actuator_destination =
				(int*)malloc(actuators*sizeof(int));
		}
		memcpy((void*)dest->actuator_destination,
			   (void*)source->actuator_destination,
			   actuators*sizeof(int));
	}

	/* copy each ADF_module */
	for (m = 0; m < min_ADF_modules+1; m++) {
		gprc_copy_ADF_module(source, dest, m,
							 rows, columns, connections_per_gene,
							 sensors, actuators);
	}
}

/* copies a chromosome from the parent to the child */
static void gprc_copy_chromosome(gprc_function * parent,
								 gprc_function * child,
								 int sensors,
								 int rows, int columns,
								 int connections_per_gene,
								 int ADF_module,
								 int chromosome_index, int chromosomes)
{
	int col,row,n,i,j;
	int start_row = chromosome_index * rows / chromosomes;
	int end_row = (chromosome_index+1) * rows / chromosomes;
	float * parent_gene = parent->genome[ADF_module].gene;
	float * child_gene = child->genome[ADF_module].gene;
	float * parent_state = parent->genome[ADF_module].state;
	float * child_state = child->genome[ADF_module].state;

	for (col = 0; col < columns; col++) {
		i = (col*rows) + start_row;
		n = i * GPRC_GENE_SIZE(connections_per_gene);
		for (row = start_row; row < end_row; row++, i++,
				 n += GPRC_GENE_SIZE(connections_per_gene)) {
			/* copy the genome */
			for (j = 0;
				 j < GPRC_GENE_SIZE(connections_per_gene); j++) {
				child_gene[n+j] = parent_gene[n+j];
			}
		
			/* copy the state */
			child_state[i] = parent_state[i];
		}
	}
}

/* Produces a child by crossing the parents */
void gprc_crossover(gprc_function *parent1, gprc_function *parent2,
					int rows, int columns,
					int sensors, int actuators,
					int connections_per_gene,
					int chromosomes,
					int allocate_memory,
					gprc_function *child)
{
	int crossover_point,i,c,n,m,min_ADF_modules,p1,p2;
	float * parent1_gene, * parent2_gene, * child_gene;
	gprc_function * parent;

	min_ADF_modules = parent1->ADF_modules;
	if (parent2->ADF_modules < min_ADF_modules) {
		min_ADF_modules = parent2->ADF_modules;
	}

	/* initialise the child */
	if (allocate_memory > 0) {
		gprc_init(child,
				  rows, columns,
				  sensors, actuators,
				  connections_per_gene,
				  min_ADF_modules, &parent1->random_seed);
	}

	for (m = 0; m < min_ADF_modules+1; m++) {

		if (m == 0) {
			/* copy chromosomes from the parents */
			for (c = 0; c < chromosomes; c++) {
				if (rand_num(&parent1->random_seed)%10000 > 5000) {
					/* copy chromosome from the first parent */
					gprc_copy_chromosome(parent1, child, sensors,
										 rows, columns,
										 connections_per_gene,
										 m, c, chromosomes);
				}
				else {
					/* copy chromosome from the second parent */
					gprc_copy_chromosome(parent2, child, sensors,
										 rows, columns,
										 connections_per_gene,
										 m, c, chromosomes);
				}
			}
		}
		else {

			/* does the child's main program contain a
			   reference to this ADF? */
			if (gprc_contains_ADFs(child,
								   0, m, rows, columns,
								   connections_per_gene,
								   sensors)==-1) {
				continue;
			}

			/* copy ADF_module from one parent or the other */
			p1 = gprc_contains_ADFs(parent1,
									0, m, rows, columns,
									connections_per_gene,
									sensors);
			p2 = gprc_contains_ADFs(parent2,
									0, m, rows, columns,
									connections_per_gene,
									sensors);

			if ((p1!=-1) || (p2!=-1)) {
				/* determine which parent to copy from */
				if (p1==-1) {
					parent = parent2;
				}
				else {
					if (p2==-1) {
						parent = parent1;
					}
					else {
						if (rand_num(&parent1->random_seed)%10000 >
							5000) {
							parent = parent2;
						}
						else {
							parent = parent1;
						}
					}
				}
				gprc_copy_ADF_module(parent, child, m,
									 rows, columns,
									 connections_per_gene,
									 sensors, actuators);
			}
		}
	}

	parent1_gene = parent1->genome[0].gene;
	parent2_gene = parent2->genome[0].gene;
	child_gene = child->genome[0].gene;

	/* actuators */
	n = rows*columns*GPRC_GENE_SIZE(connections_per_gene);
	for (i = 0; i < actuators; i++, n++) {
		if (rand_num(&parent1->random_seed)%10000>5000) {
			child_gene[n] = parent1_gene[n];
		}
		else {
			child_gene[n] = parent2_gene[n];
		}
	}

	/* cross over the sensor sources */
	if ((parent1->no_of_sensor_sources > 0) &&
		(parent2->no_of_sensor_sources > 0) &&
		(child->no_of_sensor_sources > 0)) {
		crossover_point = rand_num(&parent1->random_seed)%sensors;
		for (i = 0; i < sensors; i++) {
			if (i < sensors/2) {
				child->sensor_source[crossover_point] =
					parent1->sensor_source[crossover_point];
			}
			else {
				child->sensor_source[crossover_point] =
					parent2->sensor_source[crossover_point];
			}
			crossover_point++;
			if (crossover_point>=sensors) {
				crossover_point = 0;
			}
		}		
	}

	/* cross over the actuator destinations */
	if ((parent1->no_of_actuator_destinations > 0) &&
		(parent2->no_of_actuator_destinations > 0) &&
		(child->no_of_actuator_destinations > 0)) {
		crossover_point = rand_num(&parent2->random_seed)%actuators;
		for (i = 0; i < actuators; i++) {
			if (i < actuators/2) {
				child->actuator_destination[crossover_point] =
					parent1->actuator_destination[crossover_point];
			}
			else {
				child->actuator_destination[crossover_point] =
					parent2->actuator_destination[crossover_point];
			}
			crossover_point++;
			if (crossover_point>=actuators) {
				crossover_point = 0;
			}
		}		
	}	
}

/* Kills an individual with the given array index within
   the given environment.
   This really just swaps the pointers for maximum efficiency */
void gprc_death(gprc_environment * population,
				int victim_index)
{
	if ((victim_index < 0) ||
		(victim_index >= population->population_size)) {
		return;
	}

	if (population->population_size > 1) {
		gprc_copy(&population->individual[population->population_size-1],
				  &population->individual[victim_index],
				  population->rows, population->columns,
				  population->connections_per_gene,
				  population->sensors, population->actuators);
	}

	population->population_size--;
}

/* Two parents mate and produce an offspring.
   Here the array indexes of the parents are given and the child index
   is returned.  Indexes of parents and children are stored within
   the mating array */
int gprc_mate_environment(gprc_environment * population,
						  int parent1_index,
						  int parent2_index,
						  float mutation_prob,
						  int use_crossover,
						  int * instruction_set, int no_of_instructions)
{
	gprc_function *parent1, *parent2;

	/* has the maximum population size been reached? */
	if (population->population_size >=
		population->max_population_size) {
		return -1;
	}

	/* is the mating array full? */
	if (population->matings >=
		population->max_population_size) {
		printf("Too many matings.  Matings value should be cleared.\n");
		return -1;
	}

	/* get the parents */
	parent1 =
		&population->individual[parent1_index];
	parent2 =
		&population->individual[parent2_index];

	/* two parents mate */
	gprc_mate(parent1, parent2,
			  population->rows,
			  population->columns,
			  population->sensors,
			  population->actuators,
			  population->connections_per_gene,
			  population->min_value,
			  population->max_value,
			  population->integers_only,
			  mutation_prob, use_crossover,
			  population->chromosomes,
			  instruction_set, no_of_instructions, 0,
			  &population->individual[population->population_size]);

	/* age of the child is zero */
	(&population->individual[population->population_size])->age=0;

	/* store the indexes of the parents */
	population->mating[population->matings*3] = parent1_index;
	population->mating[population->matings*3+1] = parent2_index;
	/* store the index of the child */
	population->mating[population->matings*3+2] =
		population->population_size;
	/* increment the number of matings during this cycle */
	population->matings++;
	/* increment the size of the population */
	population->population_size++;

	/* return the array index of the child */
	return population->population_size-1;
}

/* two parents mate and produce a child */
void gprc_mate(gprc_function *parent1, gprc_function *parent2,
			   int rows, int columns,
			   int sensors, int actuators,
			   int connections_per_gene,
			   float min_value, float max_value,
			   int integers_only,
			   float mutation_prob,
			   int use_crossover,
			   int chromosomes,
			   int * instruction_set, int no_of_instructions,
			   int allocate_memory,
			   gprc_function *child)
{
	const int max_depth = 5;

	if (use_crossover > 0) {
		/* crossover two parents */
		gprc_crossover(parent1, parent2,
					   rows, columns,
					   sensors, actuators,
					   connections_per_gene,
					   chromosomes,
					   allocate_memory,
					   child);
	}
	else {
		if (allocate_memory > 0) {
			gprc_init(child,
					  rows, columns,
					  sensors, actuators,
					  connections_per_gene,
					  parent1->ADF_modules,
					  &parent1->random_seed);
		}

		/* clone one parent or the other */
		if (rand_num(&parent1->random_seed)%10000 > 5000) {
			/* clone of the first parent */
			gprc_copy(parent1, child,
					  rows, columns,
					  connections_per_gene,
					  sensors, actuators);
		}
		else {
			/* clone of the second parent */
			gprc_copy(parent2, child,
					  rows, columns,
					  connections_per_gene,
					  sensors, actuators);
		}
	}

	/* add mutations */
	gprc_mutate(child,
				rows, columns,
				sensors, actuators,
				connections_per_gene,
				chromosomes,
				mutation_prob, mutation_prob*0.5f,
				min_value, max_value,
				integers_only,
				instruction_set, no_of_instructions);

	/* which functions are used */
	gprc_used_functions(child,
						rows, columns,
						connections_per_gene,
						sensors, actuators);

	/* compress */
	gprc_compress_ADF(child,0,-1,
					  rows, columns,
					  connections_per_gene,
					  sensors, actuators,
					  min_value, max_value, max_depth, 1);

	/* update all ADF connections */
	gprc_update_ADF_modules(child, rows, columns,
							connections_per_gene,
							sensors);
}

/* Returns a fitness histogram for the given population */
static void gprc_fitness_histogram(gprc_population * population,
								   int *histogram,
								   int histogram_levels,
								   float *min_fitness,
								   float *max_fitness)
{
	int i, index;

	*min_fitness=999999;
	*max_fitness=-999999;

	/* get the minimum and maximum fitness values */
	for (i=0;i<population->size;i++) {
		if (population->fitness[i]>*max_fitness) {
			*max_fitness = population->fitness[i];
		}
		if ((population->fitness[i]<*min_fitness) &&
			(population->fitness[i]>0)) {
			*min_fitness = population->fitness[i];
		}
	}

	/* clear the histogram */
	memset((void*)histogram,'\0',histogram_levels*sizeof(int));

	if (*max_fitness <= *min_fitness) return;

	for (i = 0; i < population->size; i++) {
		if (population->fitness[i]>0) {
			index =
				(int)((population->fitness[i] -
					   *min_fitness)*(histogram_levels-1) /
					  (*max_fitness - *min_fitness));
			histogram[index]++;
		}
	}	
}

/* Returns a fitness histogram for the given population */
static void gprc_fitness_histogram_system(gprc_system * system,
										  int *histogram,
										  int histogram_levels,
										  float *min_fitness,
										  float *max_fitness)
{
	int i, j, index;
	gprc_population * population;

	*min_fitness=999999;
	*max_fitness=-999999;

	/* get the minimum and maximum fitness values */
	for (i = 0; i < system->size; i++) {
		population = &system->island[i];
		for (j=0;j<population->size;j++) {
			if (population->fitness[i]>*max_fitness) {
				*max_fitness = population->fitness[i];
			}
			if ((population->fitness[i]<*min_fitness) &&
				(population->fitness[i]>0)) {
				*min_fitness = population->fitness[i];
			}
		}
	}

	/* clear the histogram */
	memset((void*)histogram,'\0',histogram_levels*sizeof(int));

	if (*max_fitness <= *min_fitness) return;

	for (i = 0; i < system->size; i++) {
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			if (population->fitness[i]>0) {
				index =
					(int)((population->fitness[i] -
						   *min_fitness)*(histogram_levels-1) /
						  (*max_fitness - *min_fitness));			
				histogram[index]++;
			}
		}
	}
}

/* returns a value in the range 0.0 - 1.0 indicating the
   diversity within the population */
static float gprc_diversity(gprc_population * population)
{
	int i, hits=0, histogram[GPR_HISTOGRAM_LEVELS];
	float min_fitness=0, max_fitness=0;
	float average=0,variance=0;
	float occupied_fraction=0;

	gprc_fitness_histogram(population,
						   histogram, GPR_HISTOGRAM_LEVELS,
						   &min_fitness, &max_fitness);

	/* calculate the average and occupied fraction */
	for (i = 0; i < GPR_HISTOGRAM_LEVELS; i++) {
		if (histogram[i]>0) {
			average += histogram[i];
			occupied_fraction++;
			hits++;
		}
	}
	if (hits==0) return 0;

	average /= hits;
	occupied_fraction /= GPR_HISTOGRAM_LEVELS;

	/* calculate the variance */
	if (average>0) {
		for (i = 0; i < GPR_HISTOGRAM_LEVELS; i++) {
			if (histogram[i]>0) {
				variance += fabs(histogram[i] - average);
			}
		}
		if (hits>0) {
			variance = (variance/(float)hits)/average;
		}
	}

	return occupied_fraction * (1.0f/(1.0f+variance));
}

/* Produce the next generation.
   This assumes that fitness has already been evaluated */
void gprc_generation(gprc_population * population,
					 float elitism,
					 float mutation_prob,
					 int use_crossover, unsigned int * random_seed,
					 int * instruction_set, int no_of_instructions)
{
	int i, threshold;
	float diversity,mutation_prob_range;
	gprc_function * parent1, * parent2, * child;

	/* sort the population in order of fitness */
	gprc_sort(population);

	diversity = gprc_diversity(population);
	mutation_prob_range = (1.0f-mutation_prob)/2;
	mutation_prob +=
		mutation_prob_range -
		(mutation_prob_range*diversity);

	/* store the fitness history */
	population->history.tick++;
	if (population->history.tick >= population->history.interval) {
		population->history.log[population->history.index] =
			population->fitness[0];
		population->history.average[population->history.index] =
			gprc_average_fitness(population);
		population->history.diversity[population->history.index] =
			gprc_diversity(population)*100;
		population->history.index++;
		population->history.tick = 0;

		if (population->history.index >= GPR_MAX_HISTORY-2) {
			for (i = 0; i < GPR_MAX_HISTORY/2; i++) {
				population->history.log[i] =
					population->history.log[i*2];
				population->history.average[i] =
					population->history.average[i*2];
				population->history.diversity[i] =
					population->history.diversity[i*2];
			}
			population->history.index /= 2;
			population->history.interval *= 2;
		}
	}

	/* range checking */
	if ((elitism < 0.1f) || (elitism > 0.9f)) {
		elitism = 0.3f;
	}

	/* index setting the threshold for the fittest individuals */
	threshold = (int)((1.0f - elitism)*(population->size-1));

#pragma omp parallel for
	for (i = 0; i < population->size - threshold; i++) {
		/* randomly choose parents from the fittest
		   section of the population */
		parent1 =
			&population->individual[rand_num(random_seed)%threshold];
		parent2 =
			&population->individual[rand_num(random_seed)%threshold];

		/* produce a new child */
		child = &population->individual[threshold + i];
		gprc_mate(parent1, parent2,
				  population->rows, population->columns,
				  population->sensors, population->actuators,
				  population->connections_per_gene,
				  population->min_value, population->max_value,
				  population->integers_only,
				  mutation_prob,
				  use_crossover,
				  population->chromosomes,
				  instruction_set, no_of_instructions,
				  0, child);

		/* fitness not yet evaluated */
		population->fitness[threshold + i] = 0;

		/* reset the age of the child */
		child->age = 0;
	}

}

/* Produce the next generation for a system containing multiple
   sub-populations. This assumes that fitness has already
   been evaluated */
void gprc_generation_system(gprc_system * system,
							int migration_interval,
							float elitism,
							float mutation_prob,
							int use_crossover,
							unsigned int * random_seed,
							int * instruction_set,
							int no_of_instructions)
{
	int i, migrant_index;
	gprc_population *population1, *population2;
	int island1_index, island2_index;

#pragma omp parallel for
	for (i = 0; i < system->size; i++) {
		gprc_generation(&system->island[i],
						elitism,
						mutation_prob,
						use_crossover, random_seed,
						instruction_set, no_of_instructions);
	}

	/* sort by average fitness */
	gprc_sort_system(system);

	/* migrate individuals between islands */
	system->migration_tick--;
	if (system->migration_tick <= 0) {
		/* reset the counter */
		system->migration_tick = migration_interval;

		island1_index = rand_num(random_seed)%system->size; 
		island2_index = rand_num(random_seed)%system->size; 
		if ((island1_index != island2_index)) {
			/* migrate */
			population1 = &system->island[island1_index];
			population2 = &system->island[island2_index];

			/* pick a migrant */
			migrant_index = rand_num(random_seed)%population1->size;

			/* copy it to the island */
			gprc_copy(&population1->individual[migrant_index],
					  &population2->individual[population2->size-1],
					  population1->rows, population1->columns,
					  population1->connections_per_gene,
					  population1->sensors, population1->actuators);

			/* copy the fitness value */
			population2->fitness[population2->size-1] =
				population1->fitness[migrant_index];

			/* create a new random individual */
			population1->fitness[migrant_index] = 0;
			gprc_random(&population1->individual[migrant_index],
						population1->rows, population1->columns,
						population1->sensors, population1->actuators,
						population1->connections_per_gene,
						population1->min_value, population1->max_value,
						population1->integers_only,
						&population1->individual[migrant_index].random_seed,
						instruction_set, no_of_instructions);

			/* remove any ADF functions */
			gprc_remove_ADFs(&population1->individual[migrant_index],
							 population1->rows, population1->columns,
							 population1->connections_per_gene);

			/* ensure that any ADFs are valid */
			gprc_valid_ADFs(&population1->individual[migrant_index],
							population1->rows, population1->columns,
							population1->connections_per_gene,
							population1->sensors,
							population1->min_value,
							population1->max_value);

			/* update the used functions */
			gprc_used_functions(&population1->individual[migrant_index],
								population1->rows, population1->columns,
								population1->connections_per_gene,
								population1->sensors,
								population1->actuators);
		}
	}
}

/* save the given individual to file */
int gprc_save(gprc_function * f,
			  int rows, int columns,
			  int connections_per_gene,
			  int sensors, int actuators,
			  FILE * fp)
{
	int retval,m,act;
	float * gene;

	gprc_used_functions(f,rows,columns,
						connections_per_gene,
						sensors, actuators);

	retval = fwrite(&f->ADF_modules, sizeof(int), 1, fp);

	for (m = 0; m < f->ADF_modules+1; m++) {
		act = gprc_get_actuators(m,actuators);

		gene = f->genome[m].gene;

		retval = fwrite(gene, sizeof(float),
						(rows*columns*
						 GPRC_GENE_SIZE(connections_per_gene)) +
						act, fp);
	}

	retval = fwrite(&f->no_of_sensor_sources, sizeof(int), 1, fp);
	retval = fwrite(&f->no_of_actuator_destinations,
					sizeof(int), 1, fp);

	if (f->no_of_sensor_sources>0) {
		retval = fwrite(f->sensor_source, sizeof(int), sensors, fp);
	}

	if (f->no_of_actuator_destinations>0) {
		retval = fwrite(f->actuator_destination,
						sizeof(int), actuators, fp);
	}

	retval = fwrite(&f->random_seed, sizeof(unsigned int), 1, fp);
	return retval;
}

/* save a population */
void gprc_save_population(gprc_population * population,
						  FILE * fp)
{
	int i;

	fprintf(fp,"%d\n",population->size);
	fprintf(fp,"%d\n",population->rows);
	fprintf(fp,"%d\n",population->columns);
	fprintf(fp,"%d\n",population->sensors);
	fprintf(fp,"%d\n",population->actuators);
	fprintf(fp,"%d\n",population->connections_per_gene);
	fprintf(fp,"%.8f\n",population->min_value);
	fprintf(fp,"%.8f\n",population->max_value);
	fprintf(fp,"%d\n",population->integers_only);

	fprintf(fp,"%d\n",population->history.index);
	fprintf(fp,"%d\n",population->history.tick);
	fprintf(fp,"%d\n",population->history.interval);
	fprintf(fp,"%d\n",population->chromosomes);
	fprintf(fp,"%d\n",population->ADF_modules);
	for (i = 0; i < population->history.index; i++) {
		fprintf(fp,"%.10f\n",population->history.log[i]);
	}

	for (i = 0; i < population->size; i++) {
		gprc_save(&population->individual[i],
				  population->rows, population->columns,
				  population->connections_per_gene,
				  population->sensors, population->actuators,
				  fp);		
	}
}

/* save the system to file */
void gprc_save_system(gprc_system *system, FILE * fp)
{
	int i;

	/* save population parameters */
	fprintf(fp,"%d\n",system->size);
	fprintf(fp,"%d\n",system->migration_tick);
	for (i = 0; i < system->size; i++) {
		/* save the population */
		gprc_save_population((gprc_population*)&system->island[i],fp);
	}
}

/* load an individual from file */
int gprc_load(gprc_function * f,
			  int rows, int columns,
			  int connections_per_gene,
			  int sensors, int actuators,
			  FILE * fp)
{
	int retval,m,act;
	float * gene;

	retval = fread(&f->ADF_modules, sizeof(int), 1, fp);

	/* read the genome */
	for (m = 0; m < f->ADF_modules+1; m++) {
		act = gprc_get_actuators(m,actuators);

		gene = f->genome[m].gene;

		retval = fread(gene, sizeof(float),
					   (rows*columns*
						GPRC_GENE_SIZE(connections_per_gene)) +
					   act, fp);
	}

	/* read the number of sensor sources and actuator destinations */
	retval = fread(&f->no_of_sensor_sources, sizeof(int), 1, fp);
	retval = fread(&f->no_of_actuator_destinations, sizeof(int), 1, fp);

	if (f->no_of_sensor_sources>0) {
		/* create the array if necessary */
		if (f->sensor_source == 0) {
			f->sensor_source = (int*)malloc(sensors*sizeof(int));
		}
		/* read the sources */
		retval = fread(f->sensor_source, sizeof(int), sensors, fp);
	}

	if (f->no_of_actuator_destinations>0) {
		/* create the array if necessary */
		if (f->actuator_destination == 0) {
			f->actuator_destination =
				(int*)malloc(actuators*sizeof(int));
		}
		/* read the destinations */
		retval =
			fread(f->actuator_destination, sizeof(int), actuators, fp);
	}

	/* read the random seed */
	retval = fread(&f->random_seed, sizeof(unsigned int), 1, fp);

	/* calculate the function usage array */
	gprc_used_functions(f,rows,columns,
						connections_per_gene,
						sensors, actuators);

	return retval;
}

/* load a population */
void gprc_load_population(gprc_population * population,
						  FILE * fp,
						  int * instruction_set,
						  int no_of_instructions)
{
	char line[256];
	int i,ctr=0;
	int size=0,rows=0,columns=0,actuators=0,sensors=0;
	int connections_per_gene=0,integers_only=0;
	int chromosomes=0,ADF_modules=0;
	float min_value=0, max_value=0;
	unsigned int random_seed = 1234;
	int history_index=0,history_interval=0,history_tick=0;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				switch(ctr) {
				case 0: {
					size = atoi(line);
					break;
				}
				case 1: {
					rows = atoi(line);
					break;
				}
				case 2: {
					columns = atoi(line);
					break;
				}
				case 3: {
					sensors = atoi(line);
					break;
				}
				case 4: {
					actuators = atoi(line);
					break;
				}
				case 5: {
					connections_per_gene = atoi(line);
					break;
				}
				case 6: {
					min_value = atof(line);
					break;
				}
				case 7: {
					max_value = atof(line);
					break;
				}
				case 8: {
					integers_only = atoi(line);					
					break;
				}
					/* index in the fitness history */
				case 9: {
					history_index = atoi(line);
					break;
				}
					/* tick in the fitness history */
				case 10: {
					history_tick = atoi(line);
					break;
				}
					/* interval in the fitness history */
				case 11: {
					history_interval = atoi(line);
					break;
				}
				case 12: {
					chromosomes = atoi(line);
					break;
				}
				case 13: {
					ADF_modules = atoi(line);
					break;
				}

				}
				if (ctr==13) break;
				ctr++;
			}
		}
	}

	if (ctr==13) {
		gprc_init_population(population,
							 size,
							 rows, columns,
							 sensors, actuators,
							 connections_per_gene,
							 ADF_modules,
							 chromosomes,
							 min_value, max_value,
							 integers_only, &random_seed,
							 instruction_set, no_of_instructions);
	}

	/* load the fitness history */
	population->history.index = history_index;
	population->history.tick = history_tick;
	population->history.interval = history_interval;
	for (i = 0; i < history_index; i++) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line) > 0) {
				population->history.log[i] = atof(line);
			}
		}
	}

	/* load the individuals */
	for (i = 0; i < size; i++) {
		gprc_load(&population->individual[i],
				  population->rows, population->columns,
				  population->connections_per_gene,
				  population->sensors, population->actuators,
				  fp);		
	}
}

/* save the environment */
void gprc_save_environment(gprc_environment * population,
						   FILE * fp)
{
	int i;

	fprintf(fp,"%d\n",population->max_population_size);
	fprintf(fp,"%d\n",population->population_size);
	fprintf(fp,"%d\n",population->rows);
	fprintf(fp,"%d\n",population->columns);
	fprintf(fp,"%d\n",population->sensors);
	fprintf(fp,"%d\n",population->actuators);
	fprintf(fp,"%d\n",population->connections_per_gene);
	fprintf(fp,"%.8f\n",population->min_value);
	fprintf(fp,"%.8f\n",population->max_value);
	fprintf(fp,"%d\n",population->integers_only);

	fprintf(fp,"%d\n",population->chromosomes);
	fprintf(fp,"%d\n",population->ADF_modules);

	for (i = 0; i < population->population_size; i++) {
		gprc_save(&population->individual[i],
				  population->rows, population->columns,
				  population->connections_per_gene,
				  population->sensors, population->actuators,
				  fp);		
	}
}

/* load an environment */
void gprc_load_environment(gprc_environment * population,
						   FILE * fp,
						   int * instruction_set,
						   int no_of_instructions)
{
	char line[256];
	int i,ctr=0;
	int max_population_size=0,population_size=0;
	int rows=0,columns=0,actuators=0,sensors=0;
	int connections_per_gene=0,integers_only=0;
	int chromosomes=0,ADF_modules=0;
	float min_value=0, max_value=0;
	unsigned int random_seed = 1234;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				switch(ctr) {
				case 0: {
					max_population_size = atoi(line);
					break;
				}
				case 1: {
					population_size = atoi(line);
					break;
				}
				case 2: {
					rows = atoi(line);
					break;
				}
				case 3: {
					columns = atoi(line);
					break;
				}
				case 4: {
					sensors = atoi(line);
					break;
				}
				case 5: {
					actuators = atoi(line);
					break;
				}
				case 6: {
					connections_per_gene = atoi(line);
					break;
				}
				case 7: {
					min_value = atof(line);
					break;
				}
				case 8: {
					max_value = atof(line);
					break;
				}
				case 9: {
					integers_only = atoi(line);					
					break;
				}
				case 10: {
					chromosomes = atoi(line);
					break;
				}
				case 11: {
					ADF_modules = atoi(line);
					break;
				}

				}
				if (ctr==11) break;
				ctr++;
			}
		}
	}

	if (ctr==11) {
		gprc_init_environment(population,
							  max_population_size,
							  population_size,
							  rows, columns,
							  sensors, actuators,
							  connections_per_gene,
							  ADF_modules,
							  chromosomes,
							  min_value, max_value,
							  integers_only, &random_seed,
							  instruction_set, no_of_instructions);
	}

	/* load the individuals */
	for (i = 0; i < population_size; i++) {
		gprc_load(&population->individual[i],
				  population->rows, population->columns,
				  population->connections_per_gene,
				  population->sensors, population->actuators,
				  fp);		
	}
}

/* load a system from file */
void gprc_load_system(gprc_system * system,
					  FILE * fp,
					  int * instruction_set, int no_of_instructions)
{
	char line[256];
	int i,ctr=0,islands=0,population_per_island=10;
	int migration_tick=0;
	int rows=5, columns=5;
	int sensors=1, actuators=1;
	int connections_per_gene=2;
	int chromosomes=1, ADF_modules=1;
	float min_value=-10, max_value=10;
	int integers_only=0;
	unsigned int random_seed = 1234;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				switch(ctr) {
				case 0: {
					islands = atoi(line);
					break;
				}
				case 1: {
					migration_tick = atoi(line);
					break;
				}
				}
				if (ctr==1) break;
				ctr++;
			}
		}
	}

	if (islands==0) return;

	/* create a system.
	   It doesn't matter what the parameters are here, because
	   they will be overwritten later */
	gprc_init_system(system,
					 islands,
					 population_per_island,
					 rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 ADF_modules, chromosomes,
					 min_value, max_value,
					 integers_only, &random_seed,
					 instruction_set, no_of_instructions);

	/* set the current tick in the migration cycle */
	system->migration_tick = migration_tick;

	/* load each population */
	for (i = 0; i < islands; i++) {
		/* deallocate the existing population */
		gprc_free_population(&system->island[i]);
		/* load the population */
		gprc_load_population(&system->island[i], fp,
							 instruction_set, no_of_instructions);
	}
}

/* arduino setup */
static void gprc_arduino_setup(FILE * fp,
							   int no_of_digital_inputs,
							   int no_of_analog_inputs,
							   int no_of_digital_outputs,
							   int no_of_analog_outputs,
							   int baud_rate,
							   int ADF_modules)
{
	int i;

	fprintf(fp,"%s","void setup() {\n");

	fprintf(fp,"%s","  analogReference(DEFAULT);\n\n");

	fprintf(fp,"%s","  // initialize serial communication:\n");
	fprintf(fp,"  Serial.begin(%d);\n\n", baud_rate);

	for (i = 0; i < ADF_modules; i++) {
		fprintf(fp,"  genome[%d] = gene%d;\n",i,i);
		fprintf(fp,"  state[%d] = state%d;\n",i,i);
	}

	fprintf(fp,"%s","  // Clear the state array\n");
	fprintf(fp,"%s","  for (int i = 0; i < sensors + ");
	fprintf(fp,"%s","(rows*columns) + actuators; i++) {\n");
	fprintf(fp,     "    for (int j = 0; j < %d; j++) {\n",
			ADF_modules);
	fprintf(fp,"%s","      state[j][i] = 0;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n\n");

	if (no_of_digital_inputs>0) {
		fprintf(fp,"%s","  // Define digital inputs\n");
		fprintf(fp,"%s","  for (int i = 0; i < ");
		fprintf(fp,"%s","digital_inputs; i++) {\n");
		fprintf(fp,"%s","    if (dInput[i]>=0) {\n");
		fprintf(fp,"%s","      pinMode(dInput[i], INPUT);\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n\n");
	}
	if (no_of_digital_outputs>0) {
		fprintf(fp,"%s","  // Define digital outputs\n");
		fprintf(fp,"%s","  for (int i = 0; i < ");
		fprintf(fp,"%s","digital_inputs; i++) {\n");
		fprintf(fp,"%s","    pinMode(dOutput[i], OUTPUT);\n");
		fprintf(fp,"%s","  }\n\n");
	}
	if (no_of_analog_outputs>0) {
		fprintf(fp,"%s","  // Define analog (PWM) outputs\n");
		fprintf(fp,"%s","  for (int i = 0; i < ");
		fprintf(fp,"%s","analog_outputs; i++) {\n");
		fprintf(fp,"%s","    pinMode(aOutput[i], OUTPUT);\n");
		fprintf(fp,"%s","  }\n\n");
	}

	fprintf(fp,"%s","}\n\n");
}

/* C program setup */
static void gprc_c_setup(FILE * fp, int ADF_modules)
{
	int i;

	fprintf(fp,"%s","static void setup()\n{\n");
	fprintf(fp,"%s","  int i,j;\n\n");

	for (i = 0; i < ADF_modules; i++) {
		fprintf(fp,"  genome[%d] = gene%d;\n",i,i);
		fprintf(fp,"  state[%d] = state%d;\n",i,i);
	}
	fprintf(fp,"%s","\n");

	fprintf(fp,"%s","  /* Clear the state array */\n");
	fprintf(fp,"%s","  for (i = 0; i < sensors + ");
	fprintf(fp,"%s","(rows*columns) + actuators; i++) {\n");
	fprintf(fp,     "    for (j = 0; j < %d; j++) {\n", ADF_modules);
	fprintf(fp,"%s","      state[j][i] = 0;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n\n");
	fprintf(fp,"%s","}\n\n");
}

/* arduino main loop */
static void gprc_arduino_main(FILE * fp, int itterations)
{
	fprintf(fp,"%s","void loop() {\n");
	fprintf(fp,"%s","  get_inputs();\n");
	fprintf(fp,     "  for (int i = 0; i < %d; i++) {\n",itterations);
	fprintf(fp,"%s","    run(0);\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  set_outputs();\n\n");
	fprintf(fp,"%s","  // Increment the time step\n");
	fprintf(fp,"%s","  tick++;\n");
	fprintf(fp,"%s","  if (tick>32000) tick=0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* read arguments from stdin */
static void gpr_c_stdin_args(FILE * fp)
{
	fprintf(fp,"%s","static int read_stdin_args(char *argv[])\n{\n");
	fprintf(fp,"%s","  char argstr[1024], value[1024], *retval;\n");
	fprintf(fp,"%s","  int i,ctr=0,field_index=0;\n");
	fprintf(fp,"%s","  FILE * fp = fopen(\"stdin\",\"r\");\n\n");

	fprintf(fp,"%s","  if (!fp) return 0;\n");
	fprintf(fp,"%s","  while (!feof(fp)) {\n");
	fprintf(fp,"%s","    retval = fgets(argstr,1023,fp);\n");
	fprintf(fp,"%s","    if (retval) {\n");
	fprintf(fp,"%s","      for (i = 0; i < strlen(argstr); i++) {\n");
	fprintf(fp,"%s","        if ((argstr[i]==' ') ||\n");
	fprintf(fp,"%s","            (argstr[i]==',') ||\n");
	fprintf(fp,"%s","            (argstr[i]==';') ||\n");
	fprintf(fp,"%s","            (i == strlen(argstr)-1)) {\n");
	fprintf(fp,"%s","          if (i == strlen(argstr)-1) {\n");
	fprintf(fp,"%s","            value[ctr++] = argstr[i];\n");
	fprintf(fp,"%s","          }\n");
	fprintf(fp,"%s","          value[ctr] = 0;\n");
	fprintf(fp,"%s","          ctr = 0;\n");
	fprintf(fp,"%s","          argv[field_index] = ");
	fprintf(fp,"%s","(char*)malloc(strlen(value)+2);\n");
	fprintf(fp,"%s","          sprintf(argv[field_index++],\n");
	fprintf(fp,"%s","\"%s\",value);\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          value[ctr++] = argstr[i];\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  fclose(fp);\n");
	fprintf(fp,"%s","  return field_index;\n");
	fprintf(fp,"%s","}\n\n");
}

/* C program main loop */
static void gprc_c_main(FILE * fp, int itterations)
{
	fprintf(fp,"%s","int main(int argc, char* argv[])\n{\n");
	fprintf(fp,"%s","  char * stdin_args[1024];\n");
	fprintf(fp,"%s","  int i,no_of_stdin_args =\n");
	fprintf(fp,"%s","    read_stdin_args(stdin_args);\n\n");
	fprintf(fp,"%s","  setup();\n\n"); 
	fprintf(fp,"%s","  if (no_of_stdin_args>0) {\n");
	fprintf(fp,"%s","    if (no_of_stdin_args != sensors) {\n");
	fprintf(fp,"%s","      printf(\"Invalid number of arguments ");
	fprintf(fp,"%s","%d/%d\\n\",no_of_stdin_args,sensors);\n");
	fprintf(fp,"%s","      for (i = 0; i < no_of_stdin_args; i++) {\n");
	fprintf(fp,"%s","        free(stdin_args[i]);\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,"%s","      return -1;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","    get_inputs(no_of_stdin_args,stdin_args);\n");
	fprintf(fp,"%s","    for (i = 0; i < no_of_stdin_args; i++) {\n");
	fprintf(fp,"%s","      free(stdin_args[i]);\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  else {\n");
	fprintf(fp,"%s","    if (argc-1 != sensors) {\n");
	fprintf(fp,"%s","      printf(\"Invalid number of arguments ");
	fprintf(fp,"%s","%d/%d\\n\",argc-1,sensors);\n");
	fprintf(fp,"%s","      return -1;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","    get_inputs(argc-1,&argv[1]);\n"); 
	fprintf(fp,"%s","  }\n\n");


	fprintf(fp,     "  for (i = 0; i < %d; i++) {\n",itterations);
	fprintf(fp,"%s","    run(0);\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  set_outputs();\n\n");
	fprintf(fp,"%s","  /* Increment the time step */\n");
	fprintf(fp,"%s","  tick++;\n");
	fprintf(fp,"%s","  if (tick>32000) tick=0;\n");
	fprintf(fp,"%s","\n  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static void gprc_c_get_inputs(FILE * fp, int integers_only)
{
	fprintf(fp,"%s","void get_inputs(int argc, char* argv[])\n{\n");
	fprintf(fp,"%s","  int i;\n\n");
	fprintf(fp,"%s","  for (i = 0; i < argc; i++) {\n");
	if (integers_only <= 0) {
		fprintf(fp,"%s","    state[0][i] = atof(argv[i]);\n");
	}
	else {
		fprintf(fp,"%s","    state[0][i] = atoi(argv[i]);\n");
	}
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","}\n\n");
}

static void gprc_arduino_get_inputs(FILE * fp,
									int no_of_digital_inputs,
									int no_of_analog_inputs,
									int digital_high)
{
	fprintf(fp,"%s","void get_inputs() {\n");

	if (no_of_digital_inputs>0) {
		fprintf(fp,"%s",
				"  for (int i = 0; i < digital_inputs; i++) {\n");
		fprintf(fp,"%s","    if (dInput[i] >= 0) {\n");
		fprintf(fp,"%s","      if (digitalRead(dInput[i]) == HIGH) {\n");
		fprintf(fp,"        state[0][i] = %d;\n",digital_high);
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","      else {\n");
		fprintf(fp,"%s","        state[0][i] = 0;\n");
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","    else {\n");
		fprintf(fp,"%s","      // Input the current time step\n");
		fprintf(fp,"%s","      state[0][i] = tick;\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	if (no_of_analog_inputs>0) {
		fprintf(fp,"%s","  for (int i = 0; i < analog_inputs; i++) {\n");
		fprintf(fp,"%s","    if (aInput[i] >= 0) {\n");
		if (no_of_digital_inputs==0) {
			fprintf(fp,"%s",
					"      state[0][i] = analogRead(aInput[i]);\n");
		}
		else {
			fprintf(fp,"%s",
					"      state[0][i+digital_inputs] = ");
			fprintf(fp,"%s","analogRead(aInput[i]);\n");
		}
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","    else {\n");
		fprintf(fp,"%s","      // Input the current time step\n");
		if (no_of_digital_inputs==0) {
			fprintf(fp,"%s","      state[0][i] = tick;\n");
		}
		else {
			fprintf(fp,"%s",
					"      state[0][i+digital_inputs] = tick;\n");
		}
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	fprintf(fp,"%s","}\n\n");
}

static void gprc_c_run(FILE * fp,
					   int integers_only,
					   int connections_per_gene,
					   int ADF_modules)
{
	fprintf(fp,"%s","void run(int ADF_module)\n{\n");

	fprintf(fp,"%s",
			"  int row,col,n=0,i=0,j,k,g,ctr,src,dest,\n");
	fprintf(fp,"%s","block_from,block_to,no_of_args;\n");
	fprintf(fp,"%s","  int sens,act;\n");
	if (ADF_modules > 0) {
		fprintf(fp,"%s","  int call_ADF_module,itt;\n");
	}
	if (integers_only > 0) {
		fprintf(fp,"%s","  int * gp;\n\n");
	}
	else {
		fprintf(fp,"%s","  float * gp;\n\n");
	}

	fprintf(fp,"%s","  if (ADF_module == 0) {\n");
	fprintf(fp,"%s","    sens = sensors;\n");
	fprintf(fp,"%s","    act = actuators;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  else {\n");
	fprintf(fp,     "    sens = %d;\n",GPRC_MAX_ADF_MODULE_SENSORS);
	fprintf(fp,"%s","    act = 1;\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  for (col = 0; col < columns; col++) {\n");
	fprintf(fp,     "    for (row = 0; row < rows; row++, ");
	fprintf(fp,     "i++,n+=%d) {\n",
			GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,"%s","      if (genome[ADF_module][n] < 0) continue;\n");
	fprintf(fp,"%s","      gp = &genome[ADF_module][n];\n");
	fprintf(fp,     "      switch((int)gp[%d]) {\n",
			GPRC_GENE_FUNCTION_TYPE);

	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_ADF);
	if (ADF_modules > 0) {
		fprintf(fp,"%s","        if (ADF_module == 0) {\n");
		fprintf(fp,"%s","          call_ADF_module = 1 + ");
		fprintf(fp,     "(abs((int)gp[%d])%c%d);\n",
				GPRC_GENE_CONSTANT,
				'%',ADF_modules);

		fprintf(fp,"%s","          no_of_args = 1 + ");
		fprintf(fp,     "(((int)gp[%d])%c%d);\n",GPRC_INITIAL,'%',
				GPRC_MAX_ADF_MODULE_SENSORS);
		fprintf(fp,     "          if (no_of_args >= %d) {\n",
				connections_per_gene);
		fprintf(fp,     "            no_of_args = %d-1;\n",
				connections_per_gene);
		fprintf(fp,"%s","          }\n");

		fprintf(fp,"%s","          for (j = 0; j < ");
		fprintf(fp,"%s","no_of_args; j++) {\n");
		fprintf(fp,"%s","            state[call_ADF_module][j] = ");
		fprintf(fp,     "state[ADF_module][(int)gp[%d+j]];\n",
				1+GPRC_INITIAL);
		fprintf(fp,"%s","          }\n");

		fprintf(fp,"%s","          for (itt = 0; itt < 2; itt++) {\n");
		fprintf(fp,"%s","            run(call_ADF_module);\n");
		fprintf(fp,"%s","          }\n");

		fprintf(fp,"%s","          state[ADF_module][sens+i] =\n");
		fprintf(fp,"%s","            state[call_ADF_module]");
		fprintf(fp,"[%d+(rows*columns)];\n",
				GPRC_MAX_ADF_MODULE_SENSORS);
		fprintf(fp,"%s","        }\n");
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");

	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_GET);
	fprintf(fp,"%s","        j = abs((int)state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]] + (int)state[ADF_module]",
			GPRC_INITIAL);
	fprintf(fp,     "[(int)gp[%d]])\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","            %(rows*columns);\n");
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,"%s","state[ADF_module][sens+j];\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");

	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_SET);
	fprintf(fp,"%s","        j = abs((int)state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]])\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","            %(rows*columns);\n");
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "gp[%d]*state[ADF_module][(int)gp[%d]];\n",
			GPRC_GENE_CONSTANT, GPRC_INITIAL);
	fprintf(fp,"%s","        state[ADF_module][sens+j] = ");
	fprintf(fp,"%s","state[ADF_module][sens+i];\n\n");
	fprintf(fp,     "        if (state[ADF_module][sens+j] > %d) {\n",
			GPR_MAX_CONSTANT);
	fprintf(fp,     "          state[ADF_module][sens+j] = %d;\n",
			GPR_MAX_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,     "        if (state[ADF_module][sens+j] < %d) {\n",
			-GPR_MAX_CONSTANT);
	fprintf(fp,     "          state[ADF_module][sens+j] = %d;\n",
			-GPR_MAX_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");

	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_VALUE);
	fprintf(fp,     "        state[ADF_module][sens+i] = (int)gp[%d];\n",
			GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_SIGMOID);


	fprintf(fp,"%s","        no_of_args =\n");
	fprintf(fp,	    "          1 + (abs((int)gp[%d])%%(%d));\n",
			GPRC_GENE_CONSTANT, connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        for (j = 0; j < no_of_args; j++) {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] +=\n");
	fprintf(fp,     "            state[ADF_module][(int)gp[%d+j]]*",
			GPRC_INITIAL);
	fprintf(fp,     "            gp[j+%d];\n",
			GPRC_INITIAL+connections_per_gene);
	fprintf(fp,"%s","        }\n\n");

	fprintf(fp,"%s","        state[ADF_module][sens+i] =\n");
	fprintf(fp,"%s","          1.0f / (1.0f + ");
	fprintf(fp,"%s","exp(-state[ADF_module][sens+i]));\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_ADD);
	fprintf(fp,"%s","        no_of_args = 1 + ");
	fprintf(fp,     "(abs((int)gp[%d])%c%d);\n",
			GPRC_GENE_CONSTANT, '%', connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        for (j = 0; j < no_of_args; j++) {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] += ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d+j]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_SUBTRACT);
	fprintf(fp,"%s","        no_of_args = 1 + ");
	fprintf(fp,     "(abs((int)gp[%d])%c%d);\n",
			GPRC_GENE_CONSTANT, '%', connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        for (j = 0; j < no_of_args; j++) {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] -= ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d+j]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_NEGATE);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "-(int)state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_MULTIPLY);
	fprintf(fp,"%s","        no_of_args = 1 + ");
	fprintf(fp,     "(abs((int)gp[%d])%c%d);\n",
			GPRC_GENE_CONSTANT, '%', connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        for (j = 0; j < no_of_args; j++) {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] *= ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d+j]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_WEIGHT);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]] * gp[%d];\n",
			GPRC_INITIAL, GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_DIVIDE);
	fprintf(fp,"%s","        if (state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]] == 0) {\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]] / ",
			GPRC_INITIAL);
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			1+GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_MODULUS);
	fprintf(fp,"%s","        if ((int)state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]] == 0) {\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)state[ADF_module][(int)gp[%d]] %% ",
			GPRC_INITIAL);
	fprintf(fp,     "(int)state[ADF_module][(int)gp[%d]];\n",
			1+GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_FLOOR);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
				GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "floor(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_AVERAGE);
	fprintf(fp,"%s","        no_of_args = 1 + ");
	fprintf(fp,     "(abs((int)gp[%d])%c%d);\n",
			GPRC_GENE_CONSTANT, '%', connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        for (j = 1; j < no_of_args; j++) {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] += ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d+j]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        state[ADF_module][sens+i] /= ");
	fprintf(fp,"%s","no_of_args;\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_NOOP1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_NOOP2);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_NOOP3);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_NOOP4);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_GREATER_THAN);
	fprintf(fp,     "        if (state[ADF_module][(int)gp[%d]] > ",
			GPRC_INITIAL);
	fprintf(fp,     "state[ADF_module][(int)gp[%d]]) {\n",
			1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)gp[%d];\n", GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        }");
	fprintf(fp,"%s","        else {");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        }");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_LESS_THAN);
	fprintf(fp,     "        if (state[ADF_module][(int)gp[%d]] < ",
			GPRC_INITIAL);
	fprintf(fp,     "state[ADF_module][(int)gp[%d]]) {\n",
			1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)gp[%d];\n", GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_EQUALS);
	fprintf(fp,     "        if (state[ADF_module][(int)gp[%d]] == ",
			GPRC_INITIAL);
	fprintf(fp,     "state[ADF_module][(int)gp[%d]]) {\n",
			1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)gp[%d];\n", GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_AND);
	fprintf(fp,"%s","        if ((state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]]>0) &&\n",GPRC_INITIAL);
	fprintf(fp,"%s","            (state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]]>0)) {\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)gp[%d];\n", GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_OR);
	fprintf(fp,"%s","        if ((state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]]>0) ||\n",GPRC_INITIAL);
	fprintf(fp,"%s","            (state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]]>0)) {\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)gp[%d];\n", GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_XOR);
	fprintf(fp,"%s","        if ((state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]]>0) !=\n",GPRC_INITIAL);
	fprintf(fp,"%s","            (state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]]>0)) {\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)gp[%d];\n", GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_NOT);
	fprintf(fp,"%s","        if (((int)state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]]) !=\n",GPRC_INITIAL);
	fprintf(fp,"%s","            ((int)state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d]])) {\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","          state[ADF_module][sens+i] = ");
	fprintf(fp,     "(int)gp[%d];\n", GPRC_GENE_CONSTANT);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_HEBBIAN);
	fprintf(fp,"%s","        no_of_args =\n");
	fprintf(fp,     "          1 + (abs((int)gp[%d])%%(%d));\n",
			GPRC_GENE_CONSTANT, connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","        for (j = 0; j < no_of_args; j++) {\n");
	fprintf(fp,"%s","          state[ADF_module][sens+i] +=\n");
	fprintf(fp,     "            state[ADF_module][(int)gp[%d+j]] *\n",
			GPRC_INITIAL);
	fprintf(fp,     "            gp[j+%d];\n",
			GPRC_INITIAL+connections_per_gene);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        for (j = 0; j < no_of_args; j++) {\n");
	fprintf(fp,     "          gp[j+%d] +=\n",
			GPRC_INITIAL+connections_per_gene);
	fprintf(fp,"%s","            state[ADF_module][sens+i] * ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d+j]] * %f;\n",
			GPRC_INITIAL,GPR_HEBBIAN_LEARNING_RATE);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_EXP);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(int)exp(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(float)exp(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_SQUARE_ROOT);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,"%s","(int)sqrt(abs(state[ADF_module]");
		fprintf(fp,     "[(int)gp[%d]]));\n",GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,"%s","(float)sqrt(fabs(state[ADF_module]");
		fprintf(fp,     "[(int)gp[%d]]));\n",GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_ABS);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(int)abs(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(float)fabs(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_SINE);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,"%s","(int)(sin(state[ADF_module]");
		fprintf(fp,     "[(int)gp[%d]])*256);\n",GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,"%s","(float)(sin(state[ADF_module]");
		fprintf(fp,     "[(int)gp[%d]])*256);\n",GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_ARCSINE);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(int)asin(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(float)asin(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COSINE);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,"%s","(int)(cos(state[ADF_module]");
		fprintf(fp,     "[(int)gp[%d]])*256);\n",GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,"%s","(float)(cos(state[ADF_module]");
		fprintf(fp,     "[(int)gp[%d]])*256);\n",GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_ARCCOSINE);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(int)acos(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(float)acos(state[ADF_module][(int)gp[%d]]);\n",
				GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_POW);
	if (integers_only > 0) {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,"%s","(int)pow(state[ADF_module]");
		fprintf(fp,     "[(int)gp[%d]],state[ADF_module]",
				GPRC_INITIAL);
		fprintf(fp,     "[(int)gp[%d]]);\n",1+GPRC_INITIAL);
	}
	else {
		fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
		fprintf(fp,     "(float)pow(state[ADF_module][(int)gp[%d]],",
				GPRC_INITIAL);
		fprintf(fp,     "state[ADF_module][(int)gp[%d]]);\n",
				1+GPRC_INITIAL);
	}
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_MIN);
	fprintf(fp,"%s","        no_of_args = 1 + ");
	fprintf(fp,     "(abs((int)gp[%d])%c%d);\n",
			GPRC_GENE_CONSTANT, '%', connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        for (j = 1; j < no_of_args; j++) {\n");
	fprintf(fp,     "          if (state[ADF_module][(int)gp[%d+j]] < ",
			GPRC_INITIAL);
	fprintf(fp,"%s","state[ADF_module][sens+i]) {\n");
	fprintf(fp,"%s","            state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d+j]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","          }\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_MAX);
	fprintf(fp,"%s","        no_of_args = 1 + ");
	fprintf(fp,"(abs((int)gp[%d])%c%d);\n",
			GPRC_GENE_CONSTANT, '%',connections_per_gene-1);
	fprintf(fp,"%s","        state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        for (j = 1; j < no_of_args; j++) {\n");
	fprintf(fp,"%s","          if (state[ADF_module]");
	fprintf(fp,     "[(int)gp[%d+j]] > state[ADF_module][sens+i]) {\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","            state[ADF_module][sens+i] = ");
	fprintf(fp,     "state[ADF_module][(int)gp[%d+j]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","          }\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_FUNCTION);
	fprintf(fp,     "        if (((int)gp[%d] > sens) && ",
			GPRC_INITIAL);
	fprintf(fp,     "((int)gp[%d] > sens)) {\n",1+GPRC_INITIAL);
    fprintf(fp,     "          src = ((int)gp[%d]-sens)*%d;\n",
			GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
    fprintf(fp,     "          dest = ((int)gp[%d]-sens)*%d;\n",
			1+GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,"%s","          genome[ADF_module][dest] = ");
	fprintf(fp,"%s","genome[ADF_module][src];\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_CONSTANT);
	fprintf(fp,
			"        if (((int)gp[%d] > sens) && ",GPRC_INITIAL);
	fprintf(fp,     "((int)gp[%d] > sens)) {\n",1+GPRC_INITIAL);
    fprintf(fp,     "          src = ((int)gp[%d]-sens)*%d;\n",
			GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
    fprintf(fp,     "          dest = ((int)gp[%d]-sens)*%d;\n",
			1+GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,"%s",
			"          genome[ADF_module][dest+1] = ");
	fprintf(fp,"%s","genome[ADF_module][src+1];\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_STATE);
	fprintf(fp,     "        state[ADF_module][(int)gp[%d]] = ",
			1+GPRC_INITIAL);
	fprintf(fp,     "state[ADF_module][(int)gp[%d]];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_BLOCK);
	fprintf(fp,     "        block_from = (int)gp[%d];\n",
			GPRC_INITIAL);
	fprintf(fp,     "        block_to = (int)gp[%d];\n",1+GPRC_INITIAL);
	fprintf(fp,"%s","        if (block_from<block_to) {\n");
	fprintf(fp,     "          block_from = (int)gp[%d];\n",
			1+GPRC_INITIAL);
	fprintf(fp,     "          block_to = (int)gp[%d];\n",GPRC_INITIAL);
	fprintf(fp,"%s","        }\n");
	fprintf(fp,     "        k = block_to - %d;\n",GPR_BLOCK_WIDTH);
	fprintf(fp,
			"        for (j = block_from - %d;\n",GPR_BLOCK_WIDTH);
	fprintf(fp,
			"             j <= block_from + %d; j++,k++) {\n",
			GPR_BLOCK_WIDTH);
	fprintf(fp,"%s","          if ((j>sens) &&\n");
	fprintf(fp,"%s","              (k>sens) &&\n");
	fprintf(fp,"%s","              (j<i) && (k<i)) {\n");
	fprintf(fp,
			"            for (g = 0; g < %d; g++) {\n",
			GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,
			"              genome[ADF_module][(j-sens)*%d + g] =\n",
			GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,
			"                genome[ADF_module][(k-sens)*%d + g];\n",
			GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,"%s","            }\n");
	fprintf(fp,"%s","          }\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_CONNECTION1);
	fprintf(fp,     "        if (gp[%d] > sens) {\n",GPRC_INITIAL);
	fprintf(fp,     "          src = ((int)gp[%d] - sens) * %d;\n",
			GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,     "          gp[%d] = genome[ADF_module][src+2];\n",
			1+GPRC_INITIAL);
	fprintf(fp,"%s","        };\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_CONNECTION2);
	fprintf(fp,     "        if (gp[%d] > sens) {\n",1+GPRC_INITIAL);
	fprintf(fp,     "          src = ((int)gp[%d] - sens) * %d;\n",
			GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,     "          gp[%d] = genome[ADF_module][src+2];\n",
			GPRC_INITIAL);
	fprintf(fp,"%s","        };\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_CONNECTION3);
	fprintf(fp,     "        if (gp[%d] > sens) {\n",1+GPRC_INITIAL);
	fprintf(fp,     "          src = ((int)gp[%d] - sens) * %d;\n",
			GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,     "          gp[%d] = genome[ADF_module][src+%d];\n",
			GPRC_INITIAL,1+GPRC_INITIAL);
	fprintf(fp,"%s","        };\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      case %d: {\n",GPR_FUNCTION_COPY_CONNECTION4);
	fprintf(fp,     "        if (gp[%d] > sens) {\n",GPRC_INITIAL);
	fprintf(fp,     "          src = ((int)gp[%d] - sens) * %d;\n",
			GPRC_INITIAL,GPRC_GENE_SIZE(connections_per_gene));
	fprintf(fp,     "          gp[%d] = genome[ADF_module][src+%d];\n",
			1+GPRC_INITIAL,1+GPRC_INITIAL);
	fprintf(fp,"%s","        };\n");
	fprintf(fp,"%s","        break;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,"%s","      /* prevent values from going ");
	fprintf(fp,"%s","out of range */\n");
	fprintf(fp,"%s","      if ((isnan(state[ADF_module][sens+i])) ");
	fprintf(fp,"%s","|| (isinf(state[ADF_module][sens+i]))) {\n");
	fprintf(fp,"%s","        state[ADF_module][sens+i] = 0;\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      if (state[ADF_module][sens+i] > %d) {\n",
			GPR_MAX_CONSTANT);
	fprintf(fp,     "        state[ADF_module][sens+i] = %d;\n",
			GPR_MAX_CONSTANT);
	fprintf(fp,"%s","      }\n");
	fprintf(fp,     "      if (state[ADF_module][sens+i] < -%d) {\n",
			GPR_MAX_CONSTANT);
	fprintf(fp,     "        state[ADF_module][sens+i] = -%d;\n",
			GPR_MAX_CONSTANT);
	fprintf(fp,"%s","      }\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n\n");
	fprintf(fp,"%s","  /* set the actuator values */\n");
	fprintf(fp,"%s","  ctr = sens + i;\n");
	fprintf(fp,"%s","  for (i = 0; i < act; i++, ctr++, n++) {\n");
	fprintf(fp,"%s","    state[ADF_module][ctr] = ");
	fprintf(fp,"%s","state[ADF_module][(int)genome[ADF_module][n]];\n");
	fprintf(fp,"%s","  }\n");


	fprintf(fp,"%s","}\n\n");
}

static void gprc_arduino_set_outputs(FILE * fp,
									 int sensors,
									 int rows,
									 int columns,
									 int no_of_digital_outputs,
									 int no_of_analog_outputs,
									 int integers_only,
									 int digital_high)
{
	fprintf(fp,"%s","void set_outputs() {\n");

	if (no_of_digital_outputs>0) {
		fprintf(fp,"%s",
				"  for (int i = 0; i < digital_outputs; i++) {\n");
		fprintf(fp,"%s","    if (dOutput[i] >= 0) {\n");
		fprintf(fp,"      if (state[0][%d+(%d*%d)+i] > 0) {\n",
				sensors,rows,columns);
		fprintf(fp,"%s","        digitalWrite(dOutput[i],HIGH);\n");
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","      else {\n");
		fprintf(fp,"%s","        digitalWrite(dOutput[i],LOW);\n");
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	if (no_of_analog_outputs>0) {
		fprintf(fp,"%s",
				"  for (int i = 0; i < analog_outputs; i++) {\n");
		fprintf(fp,"%s","    if (aOutput[i] >= 0) {\n");
		if (no_of_digital_outputs==0) {
			fprintf(fp, "%s",
					"      analogWrite(aOutput[i],");
			fprintf(fp,"abs((int)state[0][%d+(%d*%d)+i])%%256);\n",
					sensors,rows,columns);
		}
		else {
			fprintf(fp, "%s",
					"      analogWrite(aOutput[i],abs((int)state[0]");
			fprintf(fp,
					"[%d+(%d*%d)+i+digital_outputs])%%256);\n",
					sensors,rows,columns);
		}
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	fprintf(fp,"%s","}\n\n");
}

static void gprc_c_set_outputs(FILE * fp, int integers_only)
{
	fprintf(fp,"%s","void set_outputs()\n{\n");
	fprintf(fp,"%s","  int i;\n\n");

	fprintf(fp,"%s","  for (i = 0; i < actuators; i++) {\n");
	fprintf(fp,"%s","    if (i > 0) printf(\" \");\n");
	if (integers_only<=0) {
		fprintf(fp,"%s","    printf(\"%.3f\",");
		fprintf(fp,"%s","state[0][sensors+(rows*columns)+i]);\n");
	}
	else {
		fprintf(fp,"%s","    printf(\"%d\",");
		fprintf(fp,"%s","state[0][sensors+(rows*columns)+i]);\n");
	}
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  printf(\"\\n\");\n");
	fprintf(fp,"%s","}\n\n");
}

/* saves a program suitable for use on an Arduino microcontroller */
void gprc_arduino(gprc_system * system,
				  gprc_function * f,
				  int baud_rate,
				  int digital_high,
				  int * digital_inputs, int no_of_digital_inputs,
				  int * analog_inputs, int no_of_analog_inputs,
				  int * digital_outputs, int no_of_digital_outputs,
				  int * analog_outputs, int no_of_analog_outputs,
				  int itterations,
				  int dynamic,
				  FILE * fp)
{
	gprc_population * population = &system->island[0];

	gprc_arduino_base(population->rows, population->columns,
					  population->connections_per_gene,
					  population->sensors, population->actuators,
					  population->ADF_modules,
					  population->integers_only,
					  f, baud_rate, digital_high,
					  digital_inputs, no_of_digital_inputs,
					  analog_inputs, no_of_analog_inputs,
					  digital_outputs, no_of_digital_outputs,
					  analog_outputs, no_of_analog_outputs,
					  itterations, dynamic, fp);
}

/* saves a program suitable for use on an Arduino microcontroller */
void gprc_arduino_base(int rows, int columns,
					   int connections_per_gene,
					   int sensors, int actuators,
					   int ADF_modules, int integers_only,
					   gprc_function * f,
					   int baud_rate,
					   int digital_high,
					   int * digital_inputs, int no_of_digital_inputs,
					   int * analog_inputs, int no_of_analog_inputs,
					   int * digital_outputs, int no_of_digital_outputs,
					   int * analog_outputs, int no_of_analog_outputs,
					   int itterations,
					   int dynamic,
					   FILE * fp)
{
	int row,col,i,index,ctr,m,act,w;
	float * gene;
	unsigned char * used;

	/* check that the number of inputs matches the number of sensors */
	if (no_of_digital_inputs+no_of_analog_inputs!=
		sensors) {
		fprintf(stderr,"%s",
				"Number of digital and analog inputs" \
				"should equal the number of sensors\n");
		return;
	}

	/* check that the number of outputs matches the
	   number of actuators */
	if (no_of_digital_outputs+no_of_analog_outputs!=
		actuators) {
		fprintf(stderr,"%s",
				"Number of digital and analog outputs " \
				"should equal the number of actuators\n");
		return;
	}

	/* comment header */
	fprintf(fp,"%s","// Cartesian Genetic Program\n");
	fprintf(fp,"%s","// Evolved using libgpr\n");
	fprintf(fp,"%s","// https://launchpad.net/libgpr\n\n");

	fprintf(fp,"const int sensors = %d;\n",sensors);
	fprintf(fp,"const int actuators = %d;\n",actuators);
	fprintf(fp,"const int rows = %d;\n",rows);
	fprintf(fp,"const int columns = %d;\n",columns);
	fprintf(fp,"const int ADF_modules = %d;\n",
			ADF_modules);
	fprintf(fp,"const int connections_per_gene = %d;\n",
			connections_per_gene);
	if (no_of_digital_inputs > 0) {
		fprintf(fp,"const int digital_inputs = %d;\n",
				no_of_digital_inputs);
	}
	if (no_of_analog_inputs > 0) {
		fprintf(fp,"const int analog_inputs = %d;\n",
				no_of_analog_inputs);
	}
	if (no_of_digital_outputs > 0) {
		fprintf(fp,"const int digital_outputs = %d;\n",
				no_of_digital_outputs);
	}
	if (no_of_analog_outputs > 0) {
		fprintf(fp,"const int analog_outputs = %d;\n",
				no_of_analog_outputs);
	}
	if (integers_only > 0) {
		fprintf(fp,"int * genome[%d];\n",ADF_modules+1);
		fprintf(fp,"int * state[%d];\n",ADF_modules+1);
	}
	else {
		fprintf(fp,"float * genome[%d];\n",ADF_modules+1);
		fprintf(fp,"float * state[%d];\n",ADF_modules+1);
	}
	fprintf(fp,"%s","int tick = 0;\n");
	fprintf(fp,"%s","\n");

	/* digital inputs */
	if (no_of_digital_inputs>0) {
		fprintf(fp,"%s","const int dInput[] = {");
		for (i = 0; i < no_of_digital_inputs; i++) {
			fprintf(fp,"%d",digital_inputs[i]);
			if (i<no_of_digital_inputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}

	/* analog inputs */
	if (no_of_analog_inputs>0) {
		fprintf(fp,"%s","const int aInput[] = {");
		for (i = 0; i < no_of_analog_inputs; i++) {
			fprintf(fp,"A%d",analog_inputs[i]);
			if (i<no_of_analog_inputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}

	/* digital outputs */
	if (no_of_digital_outputs>0) {
		fprintf(fp,"%s","int dOutput[] = {");
		for (i = 0; i < no_of_digital_outputs; i++) {
			fprintf(fp,"%d",digital_outputs[i]);
			if (i<no_of_digital_outputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}

	/* analog outputs */
	if (no_of_analog_outputs>0) {
		fprintf(fp,"%s","int aOutput[] = {");
		for (i = 0; i < no_of_analog_outputs; i++) {
			fprintf(fp,"%d",analog_outputs[i]);
			if (i<no_of_analog_outputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}
	fprintf(fp,"%s","\n");

	/* The genome */
	for (m = 0; m < ADF_modules+1; m++) {
		gene = f->genome[m].gene;
		/*state = f->genome[m].state;*/
		used = f->genome[m].used;
		act = gprc_get_actuators(m,actuators);

		fprintf(fp,"// Genome for ADF_module %d\n",m);
		if (integers_only > 0) {
			fprintf(fp,"int gene%d[] = {",m);
		}
		else {
			fprintf(fp,"float gene%d[] = {",m);
		}
		index = 0;
		ctr=gprc_get_sensors(m,sensors);
		for (col = 0; col < columns; col++) {
			for (row = 0;
				 row < rows;
				 row++, ctr++, index +=
					 GPRC_GENE_SIZE(connections_per_gene)) {
				/* function */
				if ((dynamic > 0) || (used[ctr] != 0)) {
					fprintf(fp,"%d,", (int)gene[index]);
					if (integers_only > 0) {
						/* constant */
						fprintf(fp,"%d,", (int)gene[index+1]);
					}
					else {
						/* constant */
						fprintf(fp,"%.3f,", gene[index+1]);
					}
					for (w = 2; w < GPRC_INITIAL; w++) {
						if (integers_only > 0) {
							fprintf(fp,"%d,", (int)gene[index+2+w]);
						}
						else {
							fprintf(fp,"%.3f,", gene[index+2+w]);
						}
					}
					/* connections */
					for (w = 0; w < GPRC_WEIGHTS_PER_CONNECTION; w++) {
						for (i = 0;
							 i < connections_per_gene; i++) {
							fprintf(fp,"%d,",
									(int)gene[index + GPRC_INITIAL + i +
											  (w*connections_per_gene)]);
						}
					}
				}
				else {
					/* this function isn't used */
					fprintf(fp,"%d,%d,", -1,-1);
					for (w = 2; w < GPRC_INITIAL; w++) {
						fprintf(fp,"%d,", -1);
					}
					for (w = 0; w < GPRC_WEIGHTS_PER_CONNECTION; w++) {
						for (i = 0; i <
								 connections_per_gene; i++) {
							fprintf(fp,"%d,", -1);
						}
					}
				}
			}
		}
		for (i = 0;i < act; i++, index++) {
			fprintf(fp,"%d",(int)gene[index]);
			if (i < act-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n\n");

		/* The state array */
		fprintf(fp,"// State array for ADF_module %d\n",m);
		if (integers_only>0) {
			fprintf(fp,"%s","int ");
		}
		else {
			fprintf(fp,"%s","float ");
		}
		fprintf(fp,"state%d[%d];\n\n",m,
				gprc_get_sensors(m,sensors) +
				(rows*columns) +
				act);
	}

	gprc_arduino_get_inputs(fp,
							no_of_digital_inputs,
							no_of_analog_inputs,
							digital_high);
	gprc_c_run(fp, integers_only,
			   connections_per_gene,
			   ADF_modules);
	gprc_arduino_set_outputs(fp, sensors, rows, columns,
							 no_of_digital_outputs,
							 no_of_analog_outputs,
							 integers_only,
							 digital_high);

	gprc_arduino_setup(fp,
					   no_of_digital_inputs,
					   no_of_analog_inputs,
					   no_of_digital_outputs,
					   no_of_analog_outputs,
					   baud_rate,
					   ADF_modules);
	gprc_arduino_main(fp,itterations);
}

/* saves as a standard C program */
void gprc_c_program(gprc_system * system,
					gprc_function * f,
					int itterations,
					int dynamic,
					FILE * fp)
{
	gprc_population * population = &system->island[0];

	gprc_c_program_base(population->rows, population->columns,
						population->connections_per_gene,
						population->sensors, population->actuators,
						population->ADF_modules,
						population->integers_only,
						f, itterations, dynamic, fp);
}

/* saves as a standard C program */
void gprc_c_program_base(int rows, int columns,
						 int connections_per_gene,
						 int sensors, int actuators,
						 int ADF_modules, int integers_only,
						 gprc_function * f,
						 int itterations,
						 int dynamic,
						 FILE * fp)
{
	int row,col,i,index,ctr,m,act,w;
	float * gene/*, * state*/;
	unsigned char * used;

	/* comment header */
	fprintf(fp,"%s","/* Cartesian Genetic Program\n");
	fprintf(fp,"%s","   Evolved using libgpr\n");
	fprintf(fp,"%s","   https://launchpad.net/libgpr\n\n");

	fprintf(fp,"%s","   To compile:\n");
	fprintf(fp,"%s","gcc -Wall -std=c99 -pedantic -o ");
	fprintf(fp,"%s","agent agent.c -lm\n*/\n\n");

	fprintf(fp,"%s","#include<stdio.h>\n");
	fprintf(fp,"%s","#include<stdlib.h>\n");
	fprintf(fp,"%s","#include<string.h>\n");
	fprintf(fp,"%s","#include<math.h>\n\n");

	fprintf(fp,"const int sensors = %d;\n",sensors);
	fprintf(fp,"const int actuators = %d;\n",actuators);
	fprintf(fp,"const int rows = %d;\n",rows);
	fprintf(fp,"const int columns = %d;\n",columns);
	fprintf(fp,"const int connections_per_gene = %d;\n",
			connections_per_gene);
	if (integers_only > 0) {
		fprintf(fp,"int * genome[%d];\n",ADF_modules+1);
		fprintf(fp,"int * state[%d];\n",ADF_modules+1);
	}
	else {
		fprintf(fp,"float * genome[%d];\n",ADF_modules+1);
		fprintf(fp,"float * state[%d];\n",ADF_modules+1);
	}
	fprintf(fp,"%s","int tick = 0;\n\n");

	/* The genome */
	for (m = 0; m < ADF_modules+1; m++) {
		gene = f->genome[m].gene;
		/*state = f->genome[m].state;*/
		used = f->genome[m].used;
		act = gprc_get_actuators(m,actuators);

		fprintf(fp,"/* Genome for ADF_module %d */\n",m);
		if (integers_only > 0) {
			fprintf(fp,"int gene%d[] = {",m);
		}
		else {
			fprintf(fp,"float gene%d[] = {",m);
		}
		index = 0;
		ctr = sensors;
		for (col = 0; col < columns; col++) {
			for (row = 0;
				 row < rows;
				 row++,ctr++,index+=
					 GPRC_GENE_SIZE(connections_per_gene)) {
				/* function */
				if ((dynamic>0) || (used[ctr]!=0)) {
					fprintf(fp,"%d,", (int)gene[index]);
					if (integers_only>0) {
						/* constant */
						fprintf(fp,"%d,", (int)gene[index+1]);
					}
					else {
						/* constant */
						fprintf(fp,"%.3f,", gene[index+1]);
					}
					for (w = 2; w < GPRC_INITIAL; w++) {
						if (integers_only > 0) {
							fprintf(fp,"%d,", (int)gene[index+2+w]);
						}
						else {
							fprintf(fp,"%.3f,", gene[index+2+w]);
						}
					}
					/* connections */
					for (w = 0; w < GPRC_WEIGHTS_PER_CONNECTION; w++) {
						for (i = 0;
							 i < connections_per_gene; i++) {
							fprintf(fp,"%d,",
									(int)gene[index + GPRC_INITIAL + i +
											  (w*connections_per_gene)]);
						}
					}
				}
				else {
					/* this function isn't used */
					fprintf(fp,"%d,%d,", -1,-1);
					for (w = 2; w < GPRC_INITIAL; w++) {
						fprintf(fp,"%d,", -1);
					}
					for (w = 0; w < GPRC_WEIGHTS_PER_CONNECTION; w++) {
						for (i = 0; i <
								 connections_per_gene; i++) {
							fprintf(fp,"%d,", -1);
						}
					}
				}
			}
		}
		for (i = 0;i < act; i++, index++) {
			fprintf(fp,"%d",(int)gene[index]);
			if (i < act-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n\n");

		/* The state array */
		fprintf(fp,"/* State array for ADF_module %d */\n",m);
		if (integers_only>0) {
			fprintf(fp,"%s","int ");
		}
		else {
			fprintf(fp,"%s","float ");
		}
		fprintf(fp,"state%d[%d];\n\n",m,
				gprc_get_sensors(m,sensors) +
				(rows*columns) +
				act);
	}

	gprc_c_get_inputs(fp, integers_only);
	gprc_c_run(fp, integers_only,
			   connections_per_gene,
			   ADF_modules);
	gprc_c_set_outputs(fp,
					   integers_only);

	gprc_c_setup(fp,ADF_modules);
	gpr_c_stdin_args(fp);
	gprc_c_main(fp,itterations);
}

/* creates an instruction set suitable for
   cartesian genetic programming */
int gprc_default_instruction_set(int * instruction_set)
{
	int i;

	for (i = GPR_FUNCTION_VALUE;
		 i < GPR_FUNCTION_TYPES_CARTESIAN; i++) {
		instruction_set[i - GPR_FUNCTION_VALUE] = i;
	}
	instruction_set[GPR_FUNCTION_TYPES_CARTESIAN -
					GPR_FUNCTION_VALUE] = GPR_FUNCTION_ADF;
	return GPR_FUNCTION_TYPES_CARTESIAN - GPR_FUNCTION_VALUE + 1;
}

/* creates an instruction set containing mathematical functions
   suitable for producing an equation */
int gprc_equation_instruction_set(int * instruction_set)
{
	instruction_set[0] = GPR_FUNCTION_VALUE;
	instruction_set[1] = GPR_FUNCTION_ADD;
	instruction_set[2] = GPR_FUNCTION_SUBTRACT;
	instruction_set[3] = GPR_FUNCTION_NEGATE;
	instruction_set[4] = GPR_FUNCTION_MULTIPLY;
	instruction_set[5] = GPR_FUNCTION_DIVIDE;
	instruction_set[6] = GPR_FUNCTION_MODULUS;
	instruction_set[7] = GPR_FUNCTION_AVERAGE;
	instruction_set[8] = GPR_FUNCTION_NOOP1;
	instruction_set[9] = GPR_FUNCTION_NOOP2;
	instruction_set[10] = GPR_FUNCTION_NOOP3;
	instruction_set[11] = GPR_FUNCTION_NOOP4;
	instruction_set[12] = GPR_FUNCTION_NOOP1;
	instruction_set[13] = GPR_FUNCTION_NOOP2;
	instruction_set[14] = GPR_FUNCTION_NOOP3;
	instruction_set[15] = GPR_FUNCTION_NOOP4;
	instruction_set[16] = GPR_FUNCTION_NOOP1;
	instruction_set[17] = GPR_FUNCTION_NOOP2;
	instruction_set[18] = GPR_FUNCTION_NOOP3;
	instruction_set[19] = GPR_FUNCTION_NOOP4;
	instruction_set[20] = GPR_FUNCTION_EXP;
	instruction_set[21] = GPR_FUNCTION_SQUARE_ROOT;
	instruction_set[22] = GPR_FUNCTION_ABS;
	instruction_set[23] = GPR_FUNCTION_SINE;
	instruction_set[24] = GPR_FUNCTION_ARCSINE;
	instruction_set[25] = GPR_FUNCTION_COSINE;
	instruction_set[26] = GPR_FUNCTION_ARCCOSINE;
	instruction_set[27] = GPR_FUNCTION_POW;
	instruction_set[28] = GPR_FUNCTION_ADF;
	return 29;
}

/* creates an instruction set containing mathematical functions
   suitable for producing an equation */
int gprc_equation_dynamic_instruction_set(int * instruction_set)
{
	instruction_set[0] = GPR_FUNCTION_VALUE;
	instruction_set[1] = GPR_FUNCTION_ADD;
	instruction_set[2] = GPR_FUNCTION_SUBTRACT;
	instruction_set[3] = GPR_FUNCTION_NEGATE;
	instruction_set[4] = GPR_FUNCTION_MULTIPLY;
	instruction_set[5] = GPR_FUNCTION_DIVIDE;
	instruction_set[6] = GPR_FUNCTION_MODULUS;
	instruction_set[7] = GPR_FUNCTION_AVERAGE;
	instruction_set[8] = GPR_FUNCTION_NOOP1;
	instruction_set[9] = GPR_FUNCTION_NOOP2;
	instruction_set[10] = GPR_FUNCTION_NOOP3;
	instruction_set[11] = GPR_FUNCTION_NOOP4;
	instruction_set[12] = GPR_FUNCTION_NOOP1;
	instruction_set[13] = GPR_FUNCTION_NOOP2;
	instruction_set[14] = GPR_FUNCTION_NOOP3;
	instruction_set[15] = GPR_FUNCTION_NOOP4;
	instruction_set[16] = GPR_FUNCTION_NOOP1;
	instruction_set[17] = GPR_FUNCTION_NOOP2;
	instruction_set[18] = GPR_FUNCTION_NOOP3;
	instruction_set[19] = GPR_FUNCTION_NOOP4;
	instruction_set[20] = GPR_FUNCTION_EXP;
	instruction_set[21] = GPR_FUNCTION_SQUARE_ROOT;
	instruction_set[22] = GPR_FUNCTION_ABS;
	instruction_set[23] = GPR_FUNCTION_SINE;
	instruction_set[24] = GPR_FUNCTION_ARCSINE;
	instruction_set[25] = GPR_FUNCTION_COSINE;
	instruction_set[26] = GPR_FUNCTION_ARCCOSINE;
	instruction_set[27] = GPR_FUNCTION_POW;
	instruction_set[28] = GPR_FUNCTION_ADF;
	instruction_set[29] = GPR_FUNCTION_GET;
	instruction_set[30] = GPR_FUNCTION_SET;
	return 31;
}

/* creates an instruction set suitable for
   dynamic cartesian genetic programming */
int gprc_dynamic_instruction_set(int * instruction_set)
{
	instruction_set[0] = GPR_FUNCTION_VALUE;
	instruction_set[1] = GPR_FUNCTION_ADD;
	instruction_set[2] = GPR_FUNCTION_SUBTRACT;
	instruction_set[3] = GPR_FUNCTION_NEGATE;
	instruction_set[4] = GPR_FUNCTION_MULTIPLY;
	instruction_set[5] = GPR_FUNCTION_WEIGHT;
	instruction_set[6] = GPR_FUNCTION_DIVIDE;
	instruction_set[7] = GPR_FUNCTION_MODULUS;
	instruction_set[8] = GPR_FUNCTION_FLOOR;
	instruction_set[9] = GPR_FUNCTION_AVERAGE;
	instruction_set[10] = GPR_FUNCTION_NOOP1;
	instruction_set[11] = GPR_FUNCTION_NOOP2;
	instruction_set[12] = GPR_FUNCTION_NOOP3;
	instruction_set[13] = GPR_FUNCTION_NOOP4;
	instruction_set[14] = GPR_FUNCTION_GREATER_THAN;
	instruction_set[15] = GPR_FUNCTION_LESS_THAN;
	instruction_set[16] = GPR_FUNCTION_EQUALS;
	instruction_set[17] = GPR_FUNCTION_AND;
	instruction_set[18] = GPR_FUNCTION_OR;
	instruction_set[19] = GPR_FUNCTION_XOR;
	instruction_set[20] = GPR_FUNCTION_NOT;
	instruction_set[21] = GPR_FUNCTION_EXP;
	instruction_set[22] = GPR_FUNCTION_SQUARE_ROOT;
	instruction_set[23] = GPR_FUNCTION_ABS;
	instruction_set[24] = GPR_FUNCTION_SINE;
	instruction_set[25] = GPR_FUNCTION_ARCSINE;
	instruction_set[26] = GPR_FUNCTION_COSINE;
	instruction_set[27] = GPR_FUNCTION_ARCCOSINE;
	instruction_set[28] = GPR_FUNCTION_POW;
	instruction_set[29] = GPR_FUNCTION_MIN;
	instruction_set[30] = GPR_FUNCTION_MAX;
	instruction_set[31] = GPR_FUNCTION_COPY_FUNCTION;
	instruction_set[32] = GPR_FUNCTION_COPY_CONSTANT;
	instruction_set[33] = GPR_FUNCTION_COPY_STATE;
	instruction_set[34] = GPR_FUNCTION_COPY_BLOCK;
	instruction_set[35] = GPR_FUNCTION_COPY_CONNECTION1;
	instruction_set[36] = GPR_FUNCTION_COPY_CONNECTION2;
	instruction_set[37] = GPR_FUNCTION_COPY_CONNECTION3;
	instruction_set[38] = GPR_FUNCTION_COPY_CONNECTION4;
	instruction_set[39] = GPR_FUNCTION_HEBBIAN;
	instruction_set[40] = GPR_FUNCTION_ADF;
	instruction_set[41] = GPR_FUNCTION_GET;
	instruction_set[42] = GPR_FUNCTION_SET;
	return 43;
}

/* an instruction set for associative learning */
int gprc_associative_instruction_set(int * instruction_set)
{
	instruction_set[0] = GPR_FUNCTION_VALUE;
	instruction_set[1] = GPR_FUNCTION_HEBBIAN;
	instruction_set[2] = GPR_FUNCTION_SIGMOID;
	instruction_set[3] = GPR_FUNCTION_MULTIPLY;
	instruction_set[4] = GPR_FUNCTION_SUBTRACT;
	instruction_set[5] = GPR_FUNCTION_NOOP1;
	instruction_set[6] = GPR_FUNCTION_NOOP2;
	return 7;
}

/* creates a larger instruction set suitable for
   cartesian genetic programming */
int gprc_advanced_instruction_set(int * instruction_set)
{
	instruction_set[0] = GPR_FUNCTION_VALUE;
	instruction_set[1] = GPR_FUNCTION_ADD;
	instruction_set[2] = GPR_FUNCTION_SUBTRACT;
	instruction_set[3] = GPR_FUNCTION_NEGATE;
	instruction_set[4] = GPR_FUNCTION_MULTIPLY;
	instruction_set[5] = GPR_FUNCTION_DIVIDE;
	instruction_set[6] = GPR_FUNCTION_MODULUS;
	instruction_set[7] = GPR_FUNCTION_FLOOR;
	instruction_set[8] = GPR_FUNCTION_AVERAGE;
	instruction_set[9] = GPR_FUNCTION_NOOP1;
	instruction_set[10] = GPR_FUNCTION_NOOP2;
	instruction_set[11] = GPR_FUNCTION_NOOP3;
	instruction_set[12] = GPR_FUNCTION_NOOP4;
	instruction_set[13] = GPR_FUNCTION_GREATER_THAN;
	instruction_set[14] = GPR_FUNCTION_LESS_THAN;
	instruction_set[15] = GPR_FUNCTION_EQUALS;
	instruction_set[16] = GPR_FUNCTION_AND;
	instruction_set[17] = GPR_FUNCTION_OR;
	instruction_set[18] = GPR_FUNCTION_XOR;
	instruction_set[19] = GPR_FUNCTION_NOT;
	instruction_set[20] = GPR_FUNCTION_EXP;
	instruction_set[21] = GPR_FUNCTION_SQUARE_ROOT;
	instruction_set[22] = GPR_FUNCTION_ABS;
	instruction_set[23] = GPR_FUNCTION_SINE;
	instruction_set[24] = GPR_FUNCTION_ARCSINE;
	instruction_set[25] = GPR_FUNCTION_COSINE;
	instruction_set[26] = GPR_FUNCTION_ARCCOSINE;
	instruction_set[27] = GPR_FUNCTION_POW;
	instruction_set[28] = GPR_FUNCTION_MIN;
	instruction_set[29] = GPR_FUNCTION_MAX;
	instruction_set[30] = GPR_FUNCTION_HEBBIAN;
	instruction_set[31] = GPR_FUNCTION_ADF;
	return 32;
}

/* uses gnuplot to plot the fitness history for the given population */
int gprc_plot_history(gprc_population * population,
					  int history_type,
					  char * filename, char * title,
					  int image_width, int image_height)
{
	int index,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float value, min_fitness = 0;
	float max_fitness = 0.01f;

	sprintf(data_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 0; index < population->history.index; index++) {
		switch(history_type) {
		case GPR_HISTORY_FITNESS: {
			value = population->history.log[index];
			if (value<0) value = 0;
			break;
		}
		case GPR_HISTORY_AVERAGE: {
			value = population->history.average[index];
			break;
		}
		case GPR_HISTORY_DIVERSITY: {
			value = population->history.diversity[index];
			break;
		}
		}


		fprintf(fp,"%d    %.10f\n",
				index*population->history.interval,value);
		/* record the maximum fitnes value */
		if (value > max_fitness) {
			max_fitness = value;
		}
		if ((index==0) ||
			(value < min_fitness)) {
			min_fitness = value;
		}
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [0:%d]\n",
			population->history.index*population->history.interval);
	fprintf(fp,"set yrange [%f:%f]\n",min_fitness,max_fitness*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Generation\"\n");
	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set ylabel \"Fitness\"\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set ylabel \"Average Fitness\"\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set ylabel \"Population Diversity\"\n");
		break;
	}
	}

	fprintf(fp,"%s","set grid\n");
	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set key right top\n");
		break;
	}
	}

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"plot \"%s\" using 1:2 notitle with lines\n",
			data_filename);
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* uses gnuplot to plot the fitness history for the given system */
int gprc_plot_history_system(gprc_system * sys,
							 int history_type,
							 char * filename, char * title,
							 int image_width, int image_height)
{
	int index,i,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float min_value = 0;
	float max_value = 0.0001f;
	float value;

	sprintf(data_filename,"%s%s",
			GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",
			GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 0; index < sys->island[0].history.index; index++) {
		fprintf(fp,"%d", index*sys->island[0].history.interval);
		for (i = 0; i < sys->size; i++) {
			switch(history_type) {
			case GPR_HISTORY_FITNESS: {
				value = sys->island[i].history.log[index];
				if (value<0) value=0;
				break;
			}
			case GPR_HISTORY_AVERAGE: {
				value = sys->island[i].history.average[index];
				break;
			}
			case GPR_HISTORY_DIVERSITY: {
				value = sys->island[i].history.diversity[index];
				break;
			}
			}

			fprintf(fp,"    %.10f",value);
			/* record the maximum value */
			if (value > max_value) {
				max_value = value;
			}
			/* record the minimum value */
			if (((index==0) && (i==0)) ||
				(value < min_value)) {
				min_value = value;
			}
		}
		fprintf(fp,"%s","\n");
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [0:%d]\n",
			sys->island[0].history.index*
			sys->island[0].history.interval);
	fprintf(fp,"set yrange [%f:%f]\n",min_value,max_value*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Generation\"\n");
	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set ylabel \"Fitness\"\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set ylabel \"Average Fitness\"\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set ylabel \"Population Diversity\"\n");
		break;
	}
	}
	fprintf(fp,"%s","set grid\n");

	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set key right top\n");
		break;
	}
	}

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"%s","plot");

	for (i = 0; i < sys->size; i++) {
		fprintf(fp,
				" \"%s\" using 1:%d title \"Island %d\" with lines",
				data_filename, (i+2), i+1);
		if (i < sys->size-1) {
			fprintf(fp,"%s",",");
		}
	}
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* uses gnuplot to plot the fitness histogram
   for the given population */
int gprc_plot_fitness(gprc_population * population,
					  char * filename, char * title,
					  int image_width, int image_height)
{
	int index,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float min_fitness = 0;
	float max_fitness = 0.01f;
	int histogram_max=1;
	int histogram[GPR_HISTOGRAM_LEVELS];

	/* create the histogram */
	gprc_fitness_histogram(population, (int*)histogram,
						   GPR_HISTOGRAM_LEVELS,
						   &min_fitness, &max_fitness);

	if (max_fitness <= min_fitness) return 0;

	sprintf(data_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	for (index = 1; index < GPR_HISTOGRAM_LEVELS; index++) {
		if (index==1) {
			min_fitness = population->fitness[index];
			max_fitness = population->fitness[index]+0.001f;
		}
		else {
			if (population->fitness[index] > max_fitness) {
				max_fitness = population->fitness[index];
			}
			if (population->fitness[index] < min_fitness) {
				min_fitness = population->fitness[index];
			}
		}
	}

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 1; index < GPR_HISTOGRAM_LEVELS; index++) {
		fprintf(fp,"%f    %d\n",
				min_fitness + (index*(max_fitness-min_fitness)/
							   GPR_HISTOGRAM_LEVELS),
				histogram[index]);
		if (histogram[index]>histogram_max) {
			histogram_max = histogram[index];
		}
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [%f:%f]\n",min_fitness,max_fitness);
	fprintf(fp,"set yrange [0:%d]\n",histogram_max*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Fitness\"\n");
	fprintf(fp,"%s","set ylabel \"Instances\"\n");
	fprintf(fp,"%s","set grid\n");
	fprintf(fp,"%s","set key right top\n");

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"plot \"%s\" using 1:2 notitle with lines\n",
			data_filename);
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* uses gnuplot to plot the fitness histogram for the given system */
int gprc_plot_fitness_system(gprc_system * sys,
							 char * filename, char * title,
							 int image_width, int image_height)
{
	int index,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float min_fitness = 0;
	float max_fitness = 0.01f;
	int histogram_max=1;
	int histogram[GPR_HISTOGRAM_LEVELS];

	/* create the histogram */
	gprc_fitness_histogram_system(sys, (int*)histogram,
								  GPR_HISTOGRAM_LEVELS,
								  &min_fitness, &max_fitness);

	if (max_fitness <= min_fitness) return 0;

	sprintf(data_filename,"%s%s",
			GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",
			GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 0; index < GPR_HISTOGRAM_LEVELS; index++) {
		fprintf(fp,"%f    %d\n",
				min_fitness + (index*(max_fitness-min_fitness)/
							   GPR_HISTOGRAM_LEVELS),
				histogram[index]);
		if (histogram[index]>histogram_max) {
			histogram_max = histogram[index];
		}
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [%f:%f]\n",min_fitness,max_fitness);
	fprintf(fp,"set yrange [0:%d]\n",histogram_max*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Fitness\"\n");
	fprintf(fp,"%s","set ylabel \"Instances\"\n");
	fprintf(fp,"%s","set grid\n");
	fprintf(fp,"%s","set key right top\n");

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"plot \"%s\" using 1:2 notitle with lines\n",
			data_filename);
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* returns zero if the two functions are the same */
int gprc_functions_are_equal(gprc_function * f1,
							 gprc_function * f2,
							 int rows, int columns,
							 int connections_per_gene,
							 int modules,
							 int sensors)
{
	int i,j;
	int max =
		sensors +
		((rows*columns)*
		 GPRC_GENE_SIZE(connections_per_gene));

	for (i = sensors; i < max; i++) {
		for (j = 0; j < 1+modules; j++) {
			if (f1->genome[j].gene[i] !=
				f2->genome[j].gene[i]) {
				return -1;
			}
		}
	}
	
	return 0;
}

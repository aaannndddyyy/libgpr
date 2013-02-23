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

#include "gprcm.h"

/* use the morphology generator to specify the functions
   and constants within the main program */
static void gprcm_morphology(gprcm_function * f,
							 int rows, int columns,
							 int sensors, int actuators,
							 int connections_per_gene,
							 int ADF_modules,
							 int integers_only,
							 float min_value, float max_value,
							 int * instruction_set,
							 int no_of_instructions)
{
	int row, col, m, index, n;
	gprc_function * morphology = &f->morphology;
	gprc_function * program = &f->program;
	float constant_value, dropout_prob = 0, v;
	int function_type, con, max_con, previous_values, dynamic = 0;
	int row_centre = rows / 2;
	int col_centre = columns / 2;

	/* the number of connections which can be specified by the
	   morphology generator */
	max_con = GPRCM_MORPHOLOGY_ACTUATORS-2;
	if (connections_per_gene < max_con) {
		max_con = connections_per_gene;
	}

	gprc_clear_state(morphology,
					 GPRCM_MORPHOLOGY_ROWS,
					 GPRCM_MORPHOLOGY_COLUMNS,
					 GPRCM_MORPHOLOGY_SENSORS,
					 GPRCM_MORPHOLOGY_ACTUATORS);

	/* for the main program and each ADF */
	/*for (m = 0; m < program->ADF_modules+1; m++) {*/
	for (m = 0; m < 1; m++) {
		n = 0;
		/* for every column within the Cartesian grid */
		for (col = 0; col < columns; col++) {
			previous_values = 
				(col*rows) + gprc_get_sensors(m, sensors);
			/* for every row within the Cartesian grid */
			for (row = 0; row < rows;
				 row++, n += GPRC_GENE_SIZE(connections_per_gene)) {
				/* set the sensors */
				gprc_set_sensor(morphology, 0, row - row_centre);
				gprc_set_sensor(morphology, 1, col - col_centre);
				gprc_set_sensor(morphology, 2, m);

				/* run the morphology generator */
				if (integers_only <= 0) {
					gprc_run_float(morphology, 0,
								   GPRCM_MORPHOLOGY_ROWS,
								   GPRCM_MORPHOLOGY_COLUMNS,
								   GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
								   GPRCM_MORPHOLOGY_SENSORS,
								   GPRCM_MORPHOLOGY_ACTUATORS,
								   dropout_prob, dynamic, 0);
				}
				else {
					gprc_run_int(morphology, 0,
								 GPRCM_MORPHOLOGY_ROWS,
								 GPRCM_MORPHOLOGY_COLUMNS,
								 GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
								 GPRCM_MORPHOLOGY_SENSORS,
								 GPRCM_MORPHOLOGY_ACTUATORS,
								 dropout_prob, dynamic, 0);
				}

				/* get the function type from the
				   morphology generator */
				v =
					gprc_get_actuator(morphology, 0,
									  GPRCM_MORPHOLOGY_ROWS,
									  GPRCM_MORPHOLOGY_COLUMNS,
									  GPRCM_MORPHOLOGY_SENSORS);
				index =	(abs((int)v)%(no_of_instructions+1))-1;

				if ((index < 0) ||
					(previous_values == 0)) {
					continue;
				}

				function_type = instruction_set[index];

				if (function_type == GPR_FUNCTION_ADF) {
					continue;
				}

				/* set the function type for a gene within
				   the main program */
				program->genome[m].gene[n+GPRC_GENE_FUNCTION_TYPE] =
					function_type;

				/* get the constant value type from the
				   morphology generator */
				v =
					gprc_get_actuator(morphology, 1,
									  GPRCM_MORPHOLOGY_ROWS,
									  GPRCM_MORPHOLOGY_COLUMNS,
									  GPRCM_MORPHOLOGY_SENSORS);
				constant_value = min_value +
					fmod(fabs(v), (max_value - min_value));

				/* set the constant value for a gene
				   within the main program*/				
				if (integers_only < 1) {
					program->genome[m].gene[n+GPRC_GENE_CONSTANT] =
						constant_value;
				}
				else {
					program->genome[m].gene[n+GPRC_GENE_CONSTANT] =
						(int)constant_value;
				}

				/* get the connection */
				for (con = 0; con < max_con; con++) {
					index =
						(int)gprc_get_actuator(morphology, 2+con,
											   GPRCM_MORPHOLOGY_ROWS,
											   GPRCM_MORPHOLOGY_COLUMNS,
											   GPRCM_MORPHOLOGY_SENSORS);

					if (index >= 0) {						
						program->genome[m].gene[n+GPRC_INITIAL+con] =
							index % previous_values;
					}
					else {
						program->genome[m].gene[n+GPRC_INITIAL+con] =
							(int)program->genome[m].gene[n+GPRC_INITIAL+con] % previous_values;
					}
				}
			}
		}
	}

	gprc_unique_outputs(program, rows, columns,
						connections_per_gene,
						sensors, actuators,
						&program->random_seed);
	
	gprc_valid_logical_operators(program, rows, columns,
								 connections_per_gene,
								 sensors,
								 &program->random_seed);
	
	gprc_valid_ADFs(program, rows, columns,
					connections_per_gene,
					sensors, min_value, max_value);	
}

/* returns an instruction set used by the morphology generator */
int gprcm_morphology_instruction_set(int * instruction_set)
{
	return gprc_equation_instruction_set(instruction_set);
}

void gprcm_init(gprcm_function * f,
				int rows, int columns, int sensors, int actuators,
				int connections_per_gene, int ADF_modules,
				unsigned int * random_seed)
{
	/* create an instruction set for the morphology generator */
	f->morphology_no_of_instructions =
		gprcm_morphology_instruction_set(f->morphology_instruction_set);

	/* create the morphology generator with a fixed architecture */
	gprc_init(&f->morphology,
			  GPRCM_MORPHOLOGY_ROWS,
			  GPRCM_MORPHOLOGY_COLUMNS,
			  GPRCM_MORPHOLOGY_SENSORS,
			  GPRCM_MORPHOLOGY_ACTUATORS,
			  GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
			  0, random_seed);

	/* create the program */
	gprc_init(&f->program,
			  rows, columns, sensors, actuators,
			  connections_per_gene, ADF_modules,
			  random_seed);
}

/* free memory */
void gprcm_free(gprcm_function * f)
{
	gprc_free(&f->program);
	gprc_free(&f->morphology);
}

/* initialise sensor sources for the given system */
void gprcm_init_sensor_sources(gprcm_system * system,
							   int no_of_sensor_sources,
							   unsigned int * random_seed)
{
	int i,j,k;
	gprcm_function * f;
	gprcm_population * population;

	if (no_of_sensor_sources <= 0) return;

	for (i = 0; i < system->size; i++) {
		/* for every island in the system */
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			f = &population->individual[j];
			(&f->program)->no_of_sensor_sources = no_of_sensor_sources;
			/* create the array */
			(&f->program)->sensor_source =
				(int*)malloc(population->sensors*sizeof(int));
			/* set random values */
			for (k = 0; k < population->sensors; k++) {
				(&f->program)->sensor_source[k] =
					rand_num(random_seed)%no_of_sensor_sources;
			}
		}
	}
}

/* initialise actuator_detinations for the given system */
void gprcm_init_actuator_destinations(gprcm_system * system,
									  int no_of_actuator_destinations,
									  unsigned int * random_seed)
{
	int i,j,k;
	gprcm_function * f;
	gprcm_population * population;

	if (no_of_actuator_destinations <= 0) return;

	for (i = 0; i < system->size; i++) {
		/* for every island in the system */
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			f = &population->individual[j];
			(&f->program)->no_of_actuator_destinations =
				no_of_actuator_destinations;
			/* create the array */
			(&f->program)->actuator_destination =
				(int*)malloc(population->actuators*sizeof(int));
			/* initial random values */
			for (k = 0; k < population->actuators; k++) {
				(&f->program)->actuator_destination[k] =
					rand_num(random_seed)%no_of_actuator_destinations;
			}
		}
	}
}

/* clears the state of an individual */
void gprcm_clear_state(gprcm_function * f,
					   int rows, int columns,
					   int sensors, int actuators)
{
	/* clear the morphology generator state */
	gprc_clear_state(&f->morphology,
					 GPRCM_MORPHOLOGY_ROWS,
					 GPRCM_MORPHOLOGY_COLUMNS,
					 GPRCM_MORPHOLOGY_SENSORS,
					 GPRCM_MORPHOLOGY_ACTUATORS);

	/* clear the main program state */
	gprc_clear_state(&f->program,
					 rows, columns,
					 sensors, actuators);
}

/* creates an initial random state for an individual */
void gprcm_random(gprcm_function * f,
				  int rows, int columns,
				  int sensors, int actuators,
				  int connections_per_gene,
				  float min_value, float max_value,
				  int integers_only, unsigned int * random_seed,
				  int * instruction_set, int no_of_instructions)
{
	/* randomise the morphology */
	gprc_random(&f->morphology,
				GPRCM_MORPHOLOGY_ROWS, GPRCM_MORPHOLOGY_COLUMNS,
				GPRCM_MORPHOLOGY_SENSORS, GPRCM_MORPHOLOGY_ACTUATORS,
				GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,				
				min_value, max_value,
				integers_only, random_seed,
				(int*)f->morphology_instruction_set,
				f->morphology_no_of_instructions);

	/* randomise the main program */
	gprc_random(&f->program, rows, columns,
				sensors, actuators,
				connections_per_gene,
				min_value, max_value,
				integers_only, random_seed,
				instruction_set, no_of_instructions);

	/* overlay morphology onto the main program */
	gprcm_morphology(f, rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 f->program.ADF_modules, integers_only,
					 min_value, max_value,
					 instruction_set, no_of_instructions);
}

/* validates an individual */
int gprcm_validate(gprcm_function * f,
				   int rows, int columns,
				   int sensors, int actuators,
				   int connections_per_gene,
				   int integers_only,
				   int * instruction_set,
				   int no_of_instructions)
{
	int result;

	result = 
		gprc_validate(&f->morphology,
					  GPRCM_MORPHOLOGY_ROWS,
					  GPRCM_MORPHOLOGY_COLUMNS,
					  GPRCM_MORPHOLOGY_SENSORS,
					  GPRCM_MORPHOLOGY_ACTUATORS,
					  GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
					  integers_only,
					  (int*)f->morphology_instruction_set,
					  f->morphology_no_of_instructions);

	if (result == GPR_VALIDATE_OK) {
		return gprc_validate(&f->program,
							 rows, columns,
							 sensors, actuators,
							 connections_per_gene,
							 integers_only,
							 instruction_set,
							 no_of_instructions);
	}
	return result;
}

void gprcm_run(gprcm_function * f, gprcm_population * population,
			   float dropout_prob, int dynamic,
			   float (*custom_function)(float,float,float))
{
	if (population->integers_only<=0) {
		gprc_run_float(&f->program, 0,
					   population->rows, population->columns,
					   population->connections_per_gene,
					   population->sensors, population->actuators,
					   dropout_prob, dynamic, (*custom_function));
	}
	else {
		gprc_run_int(&f->program, 0,
					 population->rows, population->columns,
					 population->connections_per_gene,
					 population->sensors, population->actuators,
					 dropout_prob, dynamic, (*custom_function));
	}
}

void gprcm_run_environment(gprcm_function * f,
						   gprcm_environment * population,
						   float dropout_prob, int dynamic,
						   float (*custom_function)(float,float,float))
{
	if (population->integers_only<=0) {
		gprc_run_float(&f->program, 0,
					   population->rows, population->columns,
					   population->connections_per_gene,
					   population->sensors, population->actuators,
					   dropout_prob, dynamic, (*custom_function));
	}
	else {
		gprc_run_int(&f->program, 0,
					 population->rows, population->columns,
					 population->connections_per_gene,
					 population->sensors, population->actuators,
					 dropout_prob, dynamic, (*custom_function));
	}
}

void gprcm_init_population(gprcm_population * population,
						   int size,
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
	gprc_function * f;

	population->individual =
		(gprcm_function*)malloc(size*sizeof(gprcm_function));
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
		gprcm_init(&population->individual[i],
				   rows, columns, sensors, actuators,
				   connections_per_gene, ADF_modules, random_seed);

		/* initialise individuals randomly */
		gprcm_random(&population->individual[i],
					 rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 min_value, max_value,
					 integers_only, random_seed,
					 instruction_set, no_of_instructions);

		/* clear the fitness value */
		population->fitness[i] = 0;

		f = &(&population->individual[i])->program;

		/* remove any ADF functions */
		gprc_remove_ADFs(f, rows, columns,
						 connections_per_gene);

		/* ensure that any ADFs are valid */
		gprc_valid_ADFs(f, rows, columns,
						connections_per_gene,
						sensors,
						min_value, max_value);

		/* update the used functions for the main program */
		gprc_used_functions(f, rows, columns,
							connections_per_gene,
							sensors, actuators);
	}
}

/* initialise an environment */
void gprcm_init_environment(gprcm_environment * population,
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
	gprc_function * f;

	population->individual =
		(gprcm_function*)malloc(max_population_size*
								sizeof(gprcm_function));
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
		gprcm_init(&population->individual[i],
				   rows, columns, sensors, actuators,
				   connections_per_gene, ADF_modules, random_seed);

		/* initialise individuals randomly */
		gprcm_random(&population->individual[i],
					 rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 min_value, max_value,
					 integers_only, random_seed,
					 instruction_set, no_of_instructions);

		f = &(&population->individual[i])->program;

		/* remove any ADF functions */
		gprc_remove_ADFs(f, rows, columns,
						 connections_per_gene);

		/* ensure that any ADFs are valid */
		gprc_valid_ADFs(f, rows, columns,
						connections_per_gene,
						sensors,
						min_value, max_value);

		/* update the used functions for the main program */
		gprc_used_functions(f, rows, columns,
							connections_per_gene,
							sensors, actuators);
	}
}

/* free memory for the given population */
void gprcm_free_population(gprcm_population * population)
{
	for (int i = 0; i < population->size; i++) {
		gprcm_free(&population->individual[i]);
	}
	free(population->individual);
	free(population->fitness);
}

/* free memory for the given environment population */
void gprcm_free_environment(gprcm_environment * population)
{
	for (int i = 0; i < population->max_population_size; i++) {
		gprcm_free(&population->individual[i]);
	}
	free(population->individual);
	free(population->mating);
}

/* copy from the source to the destination */
void gprcm_copy(gprcm_function * source, gprcm_function * dest,
				int rows, int columns, int connections_per_gene,
				int sensors, int actuators)
{
	gprc_copy(&source->morphology, &dest->morphology,
			  GPRCM_MORPHOLOGY_ROWS,
			  GPRCM_MORPHOLOGY_COLUMNS,
			  GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
			  GPRCM_MORPHOLOGY_SENSORS,
			  GPRCM_MORPHOLOGY_ACTUATORS);

	gprc_copy(&source->program, &dest->program,
			  rows, columns, connections_per_gene,
			  sensors, actuators);
}

/* Evaluates the fitness of all individuals in the population.
   Here we use openmp to speed up the process, since each
   evaluation is independent */
void gprcm_evaluate(gprcm_population * population,
					int time_steps, int reevaluate,
					float (*evaluate_program)
					(int,gprcm_population*,int,int))
{
	int i;

#pragma omp parallel for
	for (i = 0; i < population->size; i++) {
		if ((population->fitness[i]==0) ||
			(reevaluate>0)) {
			int s;
			gprc_function * f = &(&population->individual[i])->program;
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
		(&(&population->individual[i])->program)->age++;
		if ((&(&population->individual[i])->program)->age > GPR_MAX_AGE) {
			population->fitness[i] = 0;
		}
	}
}

/* returns the highest fitness value */
float gprcm_best_fitness(gprcm_population * population)
{
	return population->fitness[0];
}

/* returns the lowest fitness value */
float gprcm_worst_fitness(gprcm_population * population)
{
	return population->fitness[population->size-1];
}

/* returns the average fitness of the population */
float gprcm_average_fitness(gprcm_population * population)
{
	int i;
	float av = 0;

	for (i = 0; i < population->size; i++) {
		av += population->fitness[i];
	}
	return av / population->size;
}

/* returns the fittest individual in the population */
gprcm_function * gprcm_best_individual(gprcm_population * population)
{
	return &population->individual[0];
}

/* set a sensor to the given value */
void gprcm_set_sensor(gprcm_function * f, int index, float value)
{
	gprc_set_sensor(&f->program, index, value);
}

/* returns the value of a sensor */
float gprcm_get_sensor(gprcm_function * f, int index)
{
	return gprc_get_sensor(&f->program, index);
}

/* returns the sensor source identifier for the given sensor */
int gprcm_get_sensor_source(gprcm_function * f, int index)
{
	return gprc_get_sensor_source(&f->program, index);
}

/* returns the value of an actuator */
float gprcm_get_actuator(gprcm_function * f, int index,
						 int rows, int columns, int sensors)
{
	return gprc_get_actuator(&f->program, index,
							 rows, columns, sensors);
}

/* returns the actuator destination identifier
   for the given actuator */
int gprcm_get_actuator_destination(gprcm_function * f, int index)
{
	return gprc_get_actuator_destination(&f->program, index);
}

/* sorts individuals in order of fitness */
void gprcm_sort(gprcm_population * population)
{
	int i, j, best;
	float max, temp_fitness;
	gprcm_function temp_individual;

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
			population->individual[i] =	population->individual[best];
			population->individual[best] = temp_individual;
		}
	}
}

/* two parents mate and produce a child */
void gprcm_mate(gprcm_function *parent1, gprcm_function *parent2,
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
				int ADF_modules,
				gprcm_function *child)
{	
	gprc_mate(&parent1->morphology, &parent2->morphology,
			  GPRCM_MORPHOLOGY_ROWS,
			  GPRCM_MORPHOLOGY_COLUMNS,
			  GPRCM_MORPHOLOGY_SENSORS,
			  GPRCM_MORPHOLOGY_ACTUATORS,
			  GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
			  min_value, max_value,
			  integers_only, mutation_prob,
			  0, 1,
			  (int*)parent1->morphology_instruction_set,
			  parent1->morphology_no_of_instructions,
			  allocate_memory,
			  &child->morphology);

	gprc_mate(&parent1->program, &parent2->program,
			  rows, columns,
			  sensors, actuators,
			  connections_per_gene,
			  min_value, max_value,
			  integers_only,
			  mutation_prob,
			  use_crossover,
			  chromosomes,
			  instruction_set, no_of_instructions,
			  allocate_memory,
			  &child->program);

	gprcm_morphology(child, rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 parent1->program.ADF_modules, integers_only,
					 min_value, max_value,
					 instruction_set, no_of_instructions);
}

/* Returns a fitness histogram for the given population */
static void gprcm_fitness_histogram(gprcm_population * population,
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

/* returns a value in the range 0.0 - 1.0 indicating the
   diversity within the population */
static float gprcm_diversity(gprcm_population * population)
{
	int i, hits=0, histogram[GPR_HISTOGRAM_LEVELS];
	float min_fitness=0, max_fitness=0;
	float average=0,variance=0;
	float occupied_fraction=0;

	gprcm_fitness_histogram(population,
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
	if (average > 0) {
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
void gprcm_generation(gprcm_population * population,
					  float elitism,
					  float mutation_prob,
					  int use_crossover, unsigned int * random_seed,
					  int * instruction_set, int no_of_instructions)
{
	int i, threshold;
	float diversity,mutation_prob_range;
	gprcm_function * parent1, * parent2, * child;

	/* sort the population in order of fitness */
	gprcm_sort(population);

	diversity = gprcm_diversity(population);
	mutation_prob_range = (1.0f - mutation_prob) / 2;
	mutation_prob +=
		mutation_prob_range -
		(mutation_prob_range*diversity);

	/* store the fitness history */
	population->history.tick++;
	if (population->history.tick >= population->history.interval) {
		population->history.log[population->history.index] =
			population->fitness[0];
		population->history.average[population->history.index] =
			gprcm_average_fitness(population);
		population->history.diversity[population->history.index] =
			gprcm_diversity(population)*100;
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
		gprcm_mate(parent1, parent2,
				   population->rows, population->columns,
				   population->sensors, population->actuators,
				   population->connections_per_gene,
				   population->min_value, population->max_value,
				   population->integers_only,
				   mutation_prob, use_crossover,
				   population->chromosomes,
				   instruction_set, no_of_instructions,
				   0, population->ADF_modules, child);

		/* fitness not yet evaluated */
		population->fitness[threshold + i] = 0;

		/* reset the age of the child */
		(&child->program)->age = 0;
	}
}

/* save the given individual to file */
int gprcm_save(gprcm_function * f,
			   int rows, int columns,
			   int connections_per_gene,
			   int sensors, int actuators,
			   FILE * fp)
{
	gprc_save(&f->morphology,
			  GPRCM_MORPHOLOGY_ROWS,
			  GPRCM_MORPHOLOGY_COLUMNS,
			  GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,				
			  GPRCM_MORPHOLOGY_SENSORS,
			  GPRCM_MORPHOLOGY_ACTUATORS,
			  fp);

	return gprc_save(&f->program,
					 rows, columns,
					 connections_per_gene,
					 sensors, actuators,
					 fp);
}

/* load an individual from file */
int gprcm_load(gprcm_function * f,
			   int rows, int columns,
			   int connections_per_gene,
			   int sensors, int actuators,
			   FILE * fp)
{
	gprc_load(&f->morphology,
			  GPRCM_MORPHOLOGY_ROWS,
			  GPRCM_MORPHOLOGY_COLUMNS,
			  GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,				
			  GPRCM_MORPHOLOGY_SENSORS,
			  GPRCM_MORPHOLOGY_ACTUATORS,
			  fp);

	return gprc_load(&f->program,
					 rows, columns,
					 connections_per_gene,
					 sensors, actuators,
					 fp);
}

/* save a population */
void gprcm_save_population(gprcm_population * population,
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
		gprcm_save(&population->individual[i],
				   population->rows, population->columns,
				   population->connections_per_gene,
				   population->sensors, population->actuators,
				   fp);		
	}
}

/* save the environment */
void gprcm_save_environment(gprcm_environment * population,
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
		gprcm_save(&population->individual[i],
				   population->rows, population->columns,
				   population->connections_per_gene,
				   population->sensors, population->actuators,
				   fp);		
	}
}

/* load a population */
void gprcm_load_population(gprcm_population * population,
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
		gprcm_init_population(population,
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
		gprcm_load(&population->individual[i],
				   population->rows, population->columns,
				   population->connections_per_gene,
				   population->sensors, population->actuators,
				   fp);		
	}
}

/* load an environment */
void gprcm_load_environment(gprcm_environment * population,
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
		gprcm_init_environment(population,
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
		gprcm_load(&population->individual[i],
				   population->rows, population->columns,
				   population->connections_per_gene,
				   population->sensors, population->actuators,
				   fp);		
	}
}

void print_gprcm(gprcm_function * f,
				 int ADF_module,
				 int rows, int columns,
				 int sensors, int actuators,
				 int connections_per_gene,
				 int integers_only)
{
	print_gprc(&f->program,
			   ADF_module,
			   rows, columns,
			   sensors, actuators,
			   connections_per_gene,
			   integers_only);
}

void gprcm_arduino(gprcm_system * system,
				   gprcm_function * f,
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
	gprcm_population * population = &system->island[0];

    gprc_arduino_base(population->rows, population->columns,
					  population->connections_per_gene,
					  population->sensors, population->actuators,
					  population->ADF_modules,
					  population->integers_only,
					  &f->program, baud_rate, digital_high,
					  digital_inputs, no_of_digital_inputs,
					  analog_inputs, no_of_analog_inputs,
					  digital_outputs, no_of_digital_outputs,
					  analog_outputs, no_of_analog_outputs,
					  itterations, dynamic, fp);
}

/* saves as a standard C program */
void gprcm_c_program(gprcm_system * system,
					 gprcm_function * f,
					 int itterations,
					 int dynamic,
					 FILE * fp)
{
	gprcm_population * population = &system->island[0];

	gprc_c_program_base(population->rows, population->columns,
						population->connections_per_gene,
						population->sensors, population->actuators,
						population->ADF_modules,
						population->integers_only,
						&f->program, itterations, dynamic, fp);
}

/* initialise a system which contains multiple sub-populations */
void gprcm_init_system(gprcm_system * system,
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
		(gprcm_population*)malloc(islands*sizeof(gprcm_population));
	system->fitness = (float*)malloc(islands*sizeof(float));

	/* clear the fitness values */
	for (i = 0; i < islands; i++) {

		/* create a population for the island */
		gprcm_init_population(&system->island[i],
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
void gprcm_free_system(gprcm_system * system)
{
	for (int i = 0; i < system->size; i++) {
		gprcm_free_population(&system->island[i]);
	}
	free(system->island);
	free(system->fitness);
}

/* evaluates a system containing multiple sub-populations */
void gprcm_evaluate_system(gprcm_system * system,
						   int time_steps, int reevaluate,
						   float (*evaluate_program)
						   (int,gprcm_population*,int,int))
{
	int i;

#pragma omp parallel for
	for (i = 0; i < system->size; i++) {
		/* evaluate the island population */
		gprcm_evaluate(&system->island[i],
					   time_steps, reevaluate,
					   (*evaluate_program));
		/* set the average fitness */
		system->fitness[i] = gprcm_average_fitness(&system->island[i]);
	}
}

/* Produce the next generation for a system containing multiple
   sub-populations. This assumes that fitness has already
   been evaluated */
void gprcm_generation_system(gprcm_system * system,
							 int migration_interval,
							 float elitism,
							 float mutation_prob,
							 int use_crossover,
							 unsigned int * random_seed,
							 int * instruction_set,
							 int no_of_instructions)
{
	int i, migrant_index;
	gprcm_population *population1, *population2;
	int island1_index, island2_index;

#pragma omp parallel for
	for (i = 0; i < system->size; i++) {
		gprcm_generation(&system->island[i],
						 elitism,
						 mutation_prob,
						 use_crossover, random_seed,
						 instruction_set, no_of_instructions);
	}

	/* sort by average fitness */
	gprcm_sort_system(system);

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
			gprcm_copy(&population1->individual[migrant_index],
					   &population2->individual[population2->size-1],
					   population1->rows, population1->columns,
					   population1->connections_per_gene,
					   population1->sensors, population1->actuators);

			/* copy the fitness value */
			population2->fitness[population2->size-1] =
				population1->fitness[migrant_index];

			/* create a new random individual */
			population1->fitness[migrant_index] = 0;
			gprcm_random(&population1->individual[migrant_index],
						 population1->rows, population1->columns,
						 population1->sensors, population1->actuators,
						 population1->connections_per_gene,
						 population1->min_value, population1->max_value,
						 population1->integers_only,
						 &((&population1->individual[migrant_index].program)->random_seed),
						 instruction_set, no_of_instructions);

			/* remove any ADF functions */
			gprc_remove_ADFs(&(&population1->individual[migrant_index])->program,
							 population1->rows, population1->columns,
							 population1->connections_per_gene);

			/* ensure that any ADFs are valid */
			gprc_valid_ADFs(&(&population1->individual[migrant_index])->program,
							population1->rows, population1->columns,
							population1->connections_per_gene,
							population1->sensors,
							population1->min_value,
							population1->max_value);

			/* update the used functions */
			gprc_used_functions(&(&population1->individual[migrant_index])->program,
								population1->rows, population1->columns,
								population1->connections_per_gene,
								population1->sensors,
								population1->actuators);
		}
	}
}

/* sorts populations in order of average fitness */
void gprcm_sort_system(gprcm_system * system)
{
	int i, j, best;
	float max, temp_fitness;
	gprcm_population temp_island;

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

/* returns the highest fitness value for the given system */
float gprcm_best_fitness_system(gprcm_system * system)
{
	return gprcm_best_fitness(&system->island[0]);
}

/* returns the fittest individual in the given system */
gprcm_function * gprcm_best_individual_system(gprcm_system * system)
{
	return gprcm_best_individual(&system->island[0]);
}

/* load a system from file */
void gprcm_load_system(gprcm_system * system,
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
	gprcm_init_system(system,
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
		gprcm_free_population(&system->island[i]);
		/* load the population */
		gprcm_load_population(&system->island[i], fp,
							  instruction_set, no_of_instructions);
	}
}

/* save the system to file */
void gprcm_save_system(gprcm_system *system, FILE * fp)
{
	int i;

	/* save population parameters */
	fprintf(fp,"%d\n",system->size);
	fprintf(fp,"%d\n",system->migration_tick);
	for (i = 0; i < system->size; i++) {
		/* save the population */
		gprcm_save_population((gprcm_population*)&system->island[i],fp);
	}
}

int gprcm_default_instruction_set(int * instruction_set)
{
	return gprc_default_instruction_set(instruction_set);
}

int gprcm_equation_instruction_set(int * instruction_set)
{
	return gprc_equation_instruction_set(instruction_set);
}

int gprcm_equation_dynamic_instruction_set(int * instruction_set)
{
	return gprc_equation_dynamic_instruction_set(instruction_set);
}

int gprcm_advanced_instruction_set(int * instruction_set)
{
	return gprc_advanced_instruction_set(instruction_set);
}

int gprcm_dynamic_instruction_set(int * instruction_set)
{
	return gprc_dynamic_instruction_set(instruction_set);
}

int gprcm_associative_instruction_set(int * instruction_set)
{
	return gprc_associative_instruction_set(instruction_set);
}

/* uses gnuplot to plot the fitness history for the given population */
int gprcm_plot_history(gprcm_population * population,
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
int gprcm_plot_history_system(gprcm_system * sys,
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

/* returns the median fitness value */
float gprcm_median_fitness(gprcm_population * population)
{
	return population->fitness[population->size/2];
}

/* uses gnuplot to plot the fitness histogram
   for the given population */
int gprcm_plot_fitness(gprcm_population * population,
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
	gprcm_fitness_histogram(population, (int*)histogram,
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

/* Returns a fitness histogram for the given population */
static void gprcm_fitness_histogram_system(gprcm_system * system,
										   int *histogram,
										   int histogram_levels,
										   float *min_fitness,
										   float *max_fitness)
{
	int i, j, index;
	gprcm_population * population;

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

/* uses gnuplot to plot the fitness histogram for the given system */
int gprcm_plot_fitness_system(gprcm_system * sys,
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
	gprcm_fitness_histogram_system(sys, (int*)histogram,
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

void gprcm_dot(gprcm_function * f, gprcm_population * population,
			   int main_program,
			   char * sensor_names[],
			   char * actuator_names[],
			   FILE * fp)
{
	int m;
	gprc_function * f2 = &f->program;

	if (main_program < 1) {
		f2 = &f->morphology;
	}

	/* ensure that all ADF values are valid */
	gprc_valid_ADFs(f2,
					population->rows, population->columns,
					population->connections_per_gene,
					population->sensors,
					population->min_value,
					population->max_value);

	fprintf(fp,"%s","digraph graphname {\n");
	for (m = 0; m < f2->ADF_modules+1; m++) {

		if ((m == 0) ||
			(gprc_contains_ADFs(f2, 0, m,
								population->rows, population->columns,
								population->connections_per_gene,
								population->sensors)>-1)) {
			fprintf(fp,"%s","subgraph {\n");
			gprc_dot_label(f2, m,
						   population->rows,
						   population->columns,
						   population->connections_per_gene,
						   population->sensors,
						   population->actuators,
						   sensor_names,
						   actuator_names,
						   fp);
			gprc_dot_links(f2, m,
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

/* returns the number of sensors for the given ADF_module */
int gprcm_get_sensors(int ADF_module, int sensors)
{
	return gprc_get_sensors(ADF_module, sensors);
}

/* get the number of actuators for the given ADF_module */
int gprcm_get_actuators(int ADF_module, int actuators)
{
	return gprc_get_actuators(ADF_module, actuators);
}

/* Tries to convert code within the given module into
   an automatically defined function */
int gprcm_compress_ADF(gprcm_function * f,
					   int ADF_module,
					   int start_index,
					   int rows, int columns,
					   int connections_per_gene,
					   int sensors, int actuators,
					   float min_value, float max_value,
					   int max_depth,
					   int random_termination)
{
	return gprc_compress_ADF(&f->program,
							 ADF_module, start_index,
							 rows, columns,
							 connections_per_gene,
							 sensors, actuators,
							 min_value, max_value,
							 max_depth,
							 random_termination);
}

/* detects a subgraph up to a certain maximum depth */
int gprcm_get_subgraph(gprcm_function * f,
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
	return gprc_get_subgraph(&f->program, ADF_module,
							 index, parent_index, connection,
							 rows, columns,
							 connections_per_gene,
							 sensors, depth,
							 max_depth, max_genes,
							 genes, no_of_genes,
							 no_of_inputs,
							 random_termination);
}

/* Two parents mate and produce an offspring.
   Here the array indexes of the parents are given and the child index
   is returned.  Indexes of parents and children are stored within
   the mating array */
int gprcm_mate_environment(gprcm_environment * population,
						   int parent1_index,
						   int parent2_index,
						   float mutation_prob,
						   int use_crossover,
						   int * instruction_set,
						   int no_of_instructions)
{
	gprcm_function *parent1, *parent2;

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
	gprcm_mate(parent1, parent2,
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
			   population->ADF_modules,
			   &population->individual[population->population_size]);

	/* age of the child is zero */
	(&(&population->individual[population->population_size])->program)->age=0;

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

/* Kills an individual with the given array index within
   the given environment.
   This really just swaps the pointers for maximum efficiency */
void gprcm_death(gprcm_environment * population,
				 int victim_index)
{
	if ((victim_index < 0) ||
		(victim_index >= population->population_size)) {
		return;
	}

	if (population->population_size > 1) {
		gprcm_copy(&population->individual[population->population_size-1],
				  &population->individual[victim_index],
				  population->rows, population->columns,
				  population->connections_per_gene,
				  population->sensors, population->actuators);
	}

	population->population_size--;
}

/* returns zero if the two functions are the same */
int gprcm_functions_are_equal(gprcm_function * f1,
							  gprcm_function * f2,
							  int rows, int columns,
							  int connections_per_gene,
							  int modules, int sensors)
{
	if (gprc_functions_are_equal(&f1->morphology,
								 &f2->morphology,
								 GPRCM_MORPHOLOGY_ROWS,
								 GPRCM_MORPHOLOGY_COLUMNS,
								 GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
								 0, GPRCM_MORPHOLOGY_SENSORS) !=0) {
		return -1;
	}
	return gprc_functions_are_equal(&f1->program,
									&f2->program,
									rows, columns,
									connections_per_gene,
									modules, sensors);
}

/* show the population within an image */
void gprcm_show_population(unsigned char * img,
						   int img_width, int img_height, int bpp,
						   gprcm_population * population)
{
	int ix, iy, i=0;
	int tx, ty, bx, by, dimension;

	if (population->size == 0) return;

	dimension = (int)sqrt(population->size);

	for (iy = 0; iy < dimension; iy++) {
		for (ix = 0; ix < dimension; ix++, i++) {
			if (i >= population->size) break;

			tx = ix * img_width / dimension;
			ty = iy * img_height / dimension;
			bx = (ix+1) * img_width / dimension;
			by = (iy+1) * img_height / dimension;

			gprc_show_genome(img, img_width, img_height, bpp,
							 tx, ty, bx, by,
							 &(&population->individual[i])->program,
							 population->rows, population->columns,
							 population->sensors,
							 population->actuators,
							 population->connections_per_gene,
							 0);
		}
	}
}

/* show the ebnvironment population within an image */
void gprcm_show_environment(unsigned char * img,
							int img_width, int img_height, int bpp,
							gprcm_environment * population)
{
	int ix, iy, i=0;
	int tx, ty, bx, by, dimension;

	if (population->population_size == 0) return;

	dimension = (int)sqrt(population->population_size);

	for (iy = 0; iy < dimension; iy++) {
		for (ix = 0; ix < dimension; ix++, i++) {
			if (i >= population->population_size) break;

			tx = ix * img_width / dimension;
			ty = iy * img_height / dimension;
			bx = (ix+1) * img_width / dimension;
			by = (iy+1) * img_height / dimension;

			gprc_show_genome(img, img_width, img_height, bpp,
							 tx, ty, bx, by,
							 &(&population->individual[i])->program,
							 population->rows, population->columns,
							 population->sensors,
							 population->actuators,
							 population->connections_per_gene,
							 0);
		}
	}
}

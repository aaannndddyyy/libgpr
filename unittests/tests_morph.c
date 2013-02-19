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

#include "tests_morph.h"

static void show_validation_message(int retval)
{
	switch(retval) {
	case GPR_VALIDATE_OK: {
		break;
	}
	case GPR_VALIDATE_FUNCTION_LESS_THAN_MIN: {
		printf("\nFunction type less than minimum\n");
		break;
	}
	case GPR_VALIDATE_FUNCTION_TYPE_TOO_LARGE: {
		printf("\nFunction type out of range\n");
		break;
	}
	case GPR_VALIDATE_CONNECTION_LESS_THAN_ZERO: {
		printf("\nConnection index less than zero\n");
		break;
	}
	case GPR_VALIDATE_CONNECTION_TOO_LARGE: {
		printf("\nConnection index out of range\n");
		break;
	}
	case GPR_VALIDATE_FUNCTION_TYPE_NOT_IN_SET:{
		printf("\nFunction type not in set\n");
		break;
	}
	case GPR_VALIDATE_ACTUATOR_CONNECTION_LESS_THAN_ZERO: {
		printf("\nActuator connection index less than zero\n");
		break;
	}
	case GPR_VALIDATE_ACTUATOR_CONNECTION_OUT_OF_RANGE: {
		printf("\nActuator connection index out of range\n");
		break;
	}
	case GPR_VALIDATE_STATE_VALUE_OUT_OF_RANGE: {
		printf("\nState value out of range -%d -> %d\n",
			   GPR_MAX_CONSTANT, GPR_MAX_CONSTANT);
		break;
	}
	case GPR_VALIDATE_ADF_WITHIN_ADF: {
		printf("\nADF occurs within ADF\n");
		break;
	}
	case GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE: {
		printf("\nADF index out of range\n");
		break;
	}
	case GPR_VALIDATE_ADF_CONNECTIONS_OUT_OF_RANGE: {
		printf("\nADF number of connections out of range\n");
		break;
	}
	case GPR_VALIDATE_ADF_NO_OF_ARGS: {
		printf("\nADF no of arguments do not match\n");
		break;
	}
	default: {
		printf("\nValidation result: %d\n",retval);
		break;
	}
	}
}

static void test_gprcm_init()
{
	gprcm_function f;
	int rows=10, columns=20, sensors=8, actuators=4, i;
	int connections_per_gene=2;
	int modules = 1;
	unsigned int random_seed = 123;
	int instruction_set[64];
	int no_of_instructions;

	printf("test_gprcm_init...");	


	for (i = 0; i < 64; i++) {
		instruction_set[i] = 9999;
	}

	/* create a morphology instruction set */
	no_of_instructions =
		gprcm_morphology_instruction_set(instruction_set);

	assert(no_of_instructions > 0);
	for (i = 0; i < no_of_instructions; i++) {
		assert(instruction_set[i] != 9999);
	}

	/* create an individual */
	gprcm_init(&f,
			   rows, columns, sensors, actuators,
			   connections_per_gene, modules,
			   &random_seed);

	/* check that arrays are not null */
	assert((&f)->program.genome[0].gene != 0);
	assert((&f)->program.genome[0].state != 0);
	assert((&f)->program.genome[0].used != 0);
	assert((&f)->program.temp_genes != 0);

	assert((&f)->morphology.genome[0].gene != 0);
	assert((&f)->morphology.genome[0].state != 0);
	assert((&f)->morphology.genome[0].used != 0);
	assert((&f)->morphology.temp_genes != 0);

	assert((&f)->program.genome[1].gene != 0);
	assert((&f)->program.genome[1].state != 0);
	assert((&f)->program.genome[1].used != 0);
	assert((&f)->program.temp_genes != 0);

	/* check that the instruction sets are the same */
	assert(no_of_instructions ==
		   f.morphology_no_of_instructions);
	for (i = 0; i < no_of_instructions; i++) {
		assert(instruction_set[i] ==
			   f.morphology_instruction_set[i]);
	}

	/* free memory */
	gprcm_free(&f);

	printf("Ok\n");
}

static void test_gprcm_random()
{
	gprcm_function f;
	int rows=10, columns=20, sensors=8, actuators=4;
	int connections_per_gene=2,retval, trial,i;
	float min_value=-10, max_value=10;
	int instruction_set[64], no_of_instructions=0;
	int integers_only=0;
	int modules = 1;
	unsigned int random_seed = 123;

	printf("test_gprcm_random...");

	/* create an instruction set */
	for (i = 0; i < 64; i++) {
		instruction_set[i] = 999;
	}
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions > 0);
	assert(no_of_instructions < 64);
	for (i = 0; i < no_of_instructions; i++) {
		assert(instruction_set[i] != 999);
	}

	for (integers_only = 0; integers_only < 2; integers_only++) {

		for (trial = 0; trial < 100; trial++) {

			/* create an individual */
			gprcm_init(&f,
					   rows, columns, sensors, actuators,
					   connections_per_gene,
					   modules, &random_seed);

			assert((&f)->program.genome[0].gene != 0);
			assert((&f)->program.genome[0].state != 0);
			assert((&f)->program.genome[0].used != 0);

			assert((&f)->morphology.genome[0].gene != 0);
			assert((&f)->morphology.genome[0].state != 0);
			assert((&f)->morphology.genome[0].used != 0);

			assert((&f)->program.genome[1].gene != 0);
			assert((&f)->program.genome[1].state != 0);
			assert((&f)->program.genome[1].used != 0);

			/* randomize it */
			gprcm_random(&f, rows, columns,
						 sensors, actuators,
						 connections_per_gene,
						 min_value, max_value,
						 integers_only, &random_seed,
						 modules,
						 instruction_set, no_of_instructions);

			/* validate the result */
			retval = gprcm_validate(&f,	rows, columns,
									sensors, actuators,
									connections_per_gene,
									integers_only,
									instruction_set,
									no_of_instructions);
			show_validation_message(retval);
			assert(retval==GPR_VALIDATE_OK);

			/* free memory */
			gprcm_free(&f);

		}
	}

	printf("Ok\n");
}

static void test_gprcm_run()
{
	gprcm_function f;
	int rows=5, columns=10, sensors=8, actuators=8;
	int connections_per_gene=2, tick, i, j, ctr, trial;
	int chromosomes=2;
	int modules=2;
	int active_actuators=0;
	float min_value=-10, max_value=10;
	int integers_only=0;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;
	gprcm_population population;
	float dropout_rate = 0.0f;

	printf("test_gprcm_run...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprcm_init_population(&population,
						  2,
						  rows, columns,
						  sensors, actuators,
						  connections_per_gene,
						  modules,
						  chromosomes,
						  min_value, max_value,
						  integers_only, &random_seed,
						  instruction_set, no_of_instructions);

	for (trial = 0; trial < 10; trial++) {

		/* create an individual */
		gprcm_init(&f,
				   rows, columns, sensors, actuators,
				   connections_per_gene, modules,
				   &random_seed);

		for (j = 0; j < 1+modules; j++) {
			assert((&f)->program.genome[j].gene != 0);
			assert((&f)->program.genome[j].state != 0);
			assert((&f)->program.genome[j].used != 0);
		}

		assert((&f)->morphology.genome[0].gene != 0);
		assert((&f)->morphology.genome[0].state != 0);
		assert((&f)->morphology.genome[0].used != 0);

		/* randomize it */
		gprcm_random(&f,
					 rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 min_value, max_value,
					 integers_only, &random_seed,
					 modules,
					 instruction_set, no_of_instructions);

		/* clear the state */
		gprcm_clear_state(&f,rows,columns,sensors,actuators);
		for (j = 0; j < 1+modules; j++) {
			for (i = 0; i < gprcm_get_sensors(j,sensors) +
					 (rows*columns) +
					 gprcm_get_actuators(j, actuators); i++) {
				assert(f.program.genome[j].state[i]==0);
			}
		}

		for (tick = 0; tick < 1000; tick++) {
			/* some random sensor values */
			for (i = 0; i < sensors; i++) {
				gprcm_set_sensor(&f,i,rand_num(&random_seed)%256);
			}

			/* run the program */
			gprcm_run(&f, &population, dropout_rate, 0, 0);
		}

		/* count the number of non-zero states */
		ctr=0;
		for (j = 0; j < 1+modules; j++) {
			for (i = gprcm_get_sensors(j, sensors);
				 i < gprcm_get_sensors(j, sensors) + (rows*columns);
				 i++) {
				if (f.program.genome[j].state[i] != 0) ctr++;
			}
		}
		if (ctr==0) {
			printf("\nGene state values are all zero\n");
		}
		assert(ctr>0);

		/* check that at least some actuators have non-zero values */
		for (j = 0; j < 1+modules; j++) {
			for (i = gprcm_get_sensors(j,sensors) + (rows*columns);
				 i < gprcm_get_sensors(j, sensors) + (rows*columns) +
					 gprcm_get_actuators(j, actuators); i++) {
				if (f.program.genome[j].state[i] != 0) {
					active_actuators++;
				}
			}
		}

		/* free memory */
		gprcm_free(&f);
	}
	if (active_actuators==0) {
		printf("\nActuator values are all zero\n");
	}
	assert(active_actuators>0);

	gprcm_free_population(&population);

	printf("Ok\n");
}

static void test_gprcm_environment()
{
	int result, i, n, population_size = 32;
	int max_population_size = 64;
	int rows=8, columns=10, sensors=8, actuators=8;
	int connections_per_gene=8;
	int chromosomes=2;
	int modules=0;
	float min_value=-10, max_value=10;
	int integers_only=0;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;
	gprcm_environment population, population2;
	int parent1_index, parent2_index, child_index, victim_index;
	float mutation_prob = 0.2f;
	FILE * fp;
	char filename[256];

	printf("test_gprcm_environment...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprcm_init_environment(&population,
						   max_population_size,
						   population_size,
						   rows, columns,
						   sensors, actuators,
						   connections_per_gene,
						   modules,
						   chromosomes,
						   min_value, max_value,
						   integers_only, &random_seed,
						   instruction_set, no_of_instructions);

	/* clear the number of matings */
	population.matings = 0;
	for (i = 0; i < 33; i++) {
		/* randomly select parents */
		parent1_index =
			rand_num(&random_seed)%population.population_size;
		parent2_index =
			rand_num(&random_seed)%population.population_size;
		/* mate */
		child_index =
			gprcm_mate_environment(&population,
								   parent1_index, parent2_index,
								   mutation_prob,1,
								   instruction_set, no_of_instructions);
		/* check the child indexes */
		if (i < 32) {
			n = population.matings-1;
			assert(child_index == population.population_size-1);
			assert(population.mating[n*3] == parent1_index);
			assert(population.mating[n*3+1] == parent2_index);
			assert(population.mating[n*3+2] == child_index);			
		}
		else {
			assert(child_index == -1);
		}
	}
	assert(population.population_size == 64);

	/* some individuals die */
	for (i = 0; i < 10; i++) {
		victim_index =
			rand_num(&random_seed)%population.population_size;
		gprcm_death(&population, victim_index);
	}
	/* check reduced population size */
	if (population.population_size != 54) {
		printf("\nPopulation size: %d\n",
			   population.population_size);
	}
	assert(population.population_size == 54);

	/* save to file */
	sprintf(filename,"%stemp_gprcm_env.dat",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"w");
	assert(fp);
	gprcm_save_environment(&population, fp);
	fclose(fp);

	/* load the file */
	fp = fopen(filename,"r");
	assert(fp);
	gprcm_load_environment(&population2, fp,
						   instruction_set, no_of_instructions);
	fclose(fp);

	/* compare the populations */
	assert(population.max_population_size ==
		   population2.max_population_size);
	assert(population.population_size ==
		   population2.population_size);
	assert(population2.matings == 0);
	for (i = 0; i < population2.population_size; i++) {
		result =
			gprcm_functions_are_equal(&population.individual[i],
									  &population2.individual[i],
									  rows, columns,
									  connections_per_gene,
									  modules, sensors);
		if (result != 0) {
			printf("\nresult %d: %d\n", i, result);
		}
		assert(result == 0);
	}

	/* free the memory */
	gprcm_free_environment(&population);
	gprcm_free_environment(&population2);

	printf("Ok\n");
}

static void test_gprcm_run_dynamic()
{
	gprcm_function f,f2;
	int rows=20, columns=20, sensors=8, actuators=8;
	int connections_per_gene=2, tick, i, j, ctr, trial;
	int chromosomes=2;
	int modules=1;
	int retval;
	float min_value=-10, max_value=10;
	int integers_only=0;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;
	gprcm_population population;
	float dropout_rate = 0.0f;

	printf("test_gprcm_run_dynamic...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_dynamic_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprcm_init_population(&population,
						  2,
						  rows, columns,
						  sensors, actuators,
						  connections_per_gene,
						  modules,
						  chromosomes,
						  min_value, max_value,
						  integers_only, &random_seed,
						  instruction_set, no_of_instructions);

	for (integers_only=0; integers_only<2; integers_only++) {

		for (trial = 0; trial < 10; trial++) {

			/* create individuals */
			gprcm_init(&f,
					   rows, columns, sensors, actuators,
					   connections_per_gene, modules,
					   &random_seed);
			gprcm_init(&f2,
					   rows, columns, sensors, actuators,
					   connections_per_gene, modules,
					   &random_seed);

			assert((&f)->program.genome[0].gene != 0);
			assert((&f)->program.genome[0].state != 0);
			assert((&f)->program.genome[0].used != 0);

			assert((&f)->morphology.genome[0].gene != 0);
			assert((&f)->morphology.genome[0].state != 0);
			assert((&f)->morphology.genome[0].used != 0);

			assert((&f)->program.genome[1].gene != 0);
			assert((&f)->program.genome[1].state != 0);
			assert((&f)->program.genome[1].used != 0);

			/* randomize it */
			gprcm_random(&f,
						 rows, columns,
						 sensors, actuators,
						 connections_per_gene,
						 min_value, max_value,
						 integers_only, &random_seed,
						 modules,
						 instruction_set, no_of_instructions);

			/* clear the state */
			gprcm_clear_state(&f,rows,columns,sensors,actuators);
			for (j = 0; j < 1+modules; j++) {
				for (i = 0; i < gprcm_get_sensors(j,sensors) +
						 (rows*columns) +
						 gprcm_get_actuators(j,actuators); i++) {
					assert(f.program.genome[j].state[i]==0);
				}
			}

			/* validate the original */
			retval = gprcm_validate(&f,
									rows, columns,
									sensors, actuators,
									connections_per_gene,
									integers_only,
									instruction_set,
									no_of_instructions);
			show_validation_message(retval);
			assert(retval==GPR_VALIDATE_OK);

			/* take a copy of the original program */
			gprcm_copy(&f, &f2, rows, columns, connections_per_gene,
					   sensors, actuators);

			/* validate the copy */
			retval = gprcm_validate(&f2,
									rows, columns,
									sensors, actuators,
									connections_per_gene,
									integers_only,
									instruction_set,
									no_of_instructions);
			show_validation_message(retval);
			assert(retval==GPR_VALIDATE_OK);

			/* check that the copy has the same genome */
			ctr=0;
			for (i = 0;
				 i < (rows*columns)*
					 GPRC_GENE_SIZE(connections_per_gene); i++) {
				for (j = 0; j < 1+modules; j++) {
					if (f.program.genome[j].gene[i] !=
						f2.program.genome[j].gene[i]) ctr++;
				}
			}
			if (ctr!=0) {
				printf("\nGenome copy failed\n");
			}
			assert(ctr==0);


			for (tick = 0; tick < 1000; tick++) {
				/* some random sensor values */
				for (i = 0; i < sensors; i++) {
					gprcm_set_sensor(&f,i,rand_num(&random_seed)%256);
				}

				/* run the program */
				gprcm_run(&f, &population, dropout_rate, 1, 0);
			}

			/* count the number of changed functions */
			ctr=0;
			for (i = 0;
				 i < (rows*columns)*
					 GPRC_GENE_SIZE(connections_per_gene); i++) {
				if (f.program.genome[0].gene[i] !=
					f2.program.genome[0].gene[i]) ctr++;
			}
			if (ctr==0) {
				printf("\nGenome values have not changed\n");
			}
			assert(ctr>0);

			/* check that the altered genome is valid */
			retval = gprcm_validate(&f,
									rows, columns,
									sensors, actuators,
									connections_per_gene,
									integers_only,
									instruction_set,
									no_of_instructions);
			show_validation_message(retval);
			assert(retval==GPR_VALIDATE_OK);

			/* free memory */
			gprcm_free(&f);
			gprcm_free(&f2);
		}
	}

	gprcm_free_population(&population);

	printf("Ok\n");
}

static void test_gprcm_copy()
{
	gprcm_function f,f2;
	int rows=10, columns=20, sensors=8, actuators=4;
	int i,connections_per_gene=2;
	float min_value=-10, max_value=10;
	int integers_only=0;
	int modules=1;
	int instruction_set[64], no_of_instructions=0;
	unsigned int random_seed = 123;

	printf("test_gprcm_copy...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create an individual */
	gprcm_init(&f,
			   rows, columns, sensors, actuators,
			   connections_per_gene, modules,
			   &random_seed);

	assert((&f)->program.genome[0].gene != 0);
	assert((&f)->program.genome[0].state != 0);
	assert((&f)->program.genome[0].used != 0);

	assert((&f)->morphology.genome[0].gene != 0);
	assert((&f)->morphology.genome[0].state != 0);
	assert((&f)->morphology.genome[0].used != 0);

	assert((&f)->program.genome[1].gene != 0);
	assert((&f)->program.genome[1].state != 0);
	assert((&f)->program.genome[1].used != 0);

	/* randomize it */
	gprcm_random(&f,
				 rows, columns,
				 sensors, actuators,
				 connections_per_gene,
				 min_value, max_value,
				 integers_only, &random_seed,
				 modules,
				 instruction_set, no_of_instructions);

	/* create a second individual */
	gprcm_init(&f2,
			   rows, columns, sensors, actuators,
			   connections_per_gene, modules,
			   &random_seed);

	/* copy from the first individual to the second */
	gprcm_copy(&f, &f2,
			   rows, columns, connections_per_gene,
			   sensors, actuators);

	/* check that the values match */
	for (i = 0;
		 i <	(rows*columns*
				 GPRC_GENE_SIZE(connections_per_gene)) +
			 actuators;i++) {
		assert(f.program.genome[0].gene[i] ==
			   f2.program.genome[0].gene[i]);
	}

	/* free memory */
	gprcm_free(&f);
	gprcm_free(&f2);

	printf("Ok\n");
}

static void test_gprcm_mutate()
{
	gprcm_function f,f2;
	int rows=10, columns=10, sensors=8, actuators=4, trial;
	int i,diff,connections_per_gene=2, percent_diff;
	int chromosomes=2;
	int modules=1;
	float min_value=-10, max_value=10;
	float mutation_prob=0.5f;
	int expected_differences, retval;
	int integers_only=0;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gprcm_mutate...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	for (integers_only = 0; integers_only < 2; integers_only++) {

		/* since this includes randomness try a few times */
		for (trial = 0; trial < 100; trial++) {

			/* create an individual */
			gprcm_init(&f,
					   rows, columns, sensors, actuators,
					   connections_per_gene, modules,
					   &random_seed);

			assert((&f)->program.genome[0].gene != 0);
			assert((&f)->program.genome[0].state != 0);
			assert((&f)->program.genome[0].used != 0);

			assert((&f)->morphology.genome[0].gene != 0);
			assert((&f)->morphology.genome[0].state != 0);
			assert((&f)->morphology.genome[0].used != 0);

			assert((&f)->program.genome[1].gene != 0);
			assert((&f)->program.genome[1].state != 0);
			assert((&f)->program.genome[1].used != 0);

			/* randomize it */
			gprcm_random(&f,
						 rows, columns,
						 sensors, actuators,
						 connections_per_gene,
						 min_value, max_value,
						 integers_only, &random_seed,
						 modules,
						 instruction_set, no_of_instructions);

			/* create a second individual */
			gprcm_init(&f2,
					   rows, columns, sensors, actuators,
					   connections_per_gene, modules,
					   &random_seed);

			/* copy from the first individual to the second */
			gprcm_copy(&f, &f2,
					   rows, columns, connections_per_gene,
					   sensors, actuators);

			/* check that the values match */
			for (i = 0;
				 i < (rows*columns*
					  GPRC_GENE_SIZE(connections_per_gene)) +
					 actuators; i++) {
				assert(f.program.genome[0].gene[i] ==
					   f2.program.genome[0].gene[i]);
			}

			/* now mutate the second individual */
			gprcm_mutate(&f2,
						 rows, columns,
						 sensors, actuators,
						 connections_per_gene,
						 chromosomes,
						 mutation_prob,0,
						 min_value, max_value,
						 integers_only,
						 instruction_set, no_of_instructions);

			/* count the differences */
			expected_differences = mutation_prob*(rows*columns);
			diff = 0;
			for (i = 0;
				 i < (rows*columns*
					  GPRC_GENE_SIZE(connections_per_gene)) +
					 gprcm_get_actuators(0,actuators); i++) {
				if (f.program.genome[0].gene[i] !=
					f2.program.genome[0].gene[i]) diff++;
			}
			for (i = 0;
				 i < (rows*columns*
					  GPRC_GENE_SIZE(connections_per_gene)) +
					 gprcm_get_actuators(1,actuators); i++) {
				if (f.program.genome[1].gene[i] !=
					f2.program.genome[1].gene[i]) diff++;
			}
			percent_diff = diff*100/expected_differences;
			if (abs(percent_diff) > 150) {
				printf("Number of mutations larger " \
					   "than expected (%d%%)\n",
					   percent_diff);
			}
			if (abs(percent_diff) < 50) {
				printf("Number of mutations smaller " \
					   "than expected (%d%%)\n",
					   percent_diff);
			}
			assert(percent_diff <= 150);
			assert(percent_diff >= 50);

			/* validate the result */
			retval = gprcm_validate(&f2,
									rows, columns,
									sensors, actuators,
									connections_per_gene,
									integers_only,
									instruction_set,
									no_of_instructions);
			show_validation_message(retval);
			assert(retval==GPR_VALIDATE_OK);

			/* free memory */
			gprcm_free(&f);
			gprcm_free(&f2);

		}
	}

	printf("Ok\n");
}


static void test_gprcm_crossover()
{
	gprcm_function parent1,parent2,child;
	int rows=10, columns=20, sensors=8, actuators=4;
	int i,connections_per_gene=2;
	int chromosomes=5;
	int modules=1;
	float min_value=-10, max_value=10;
	int instruction_set[64], no_of_instructions=0;
	int parent1_hits, parent2_hits;
	int integers_only=0;
	unsigned int random_seed = 123;

	printf("test_gprcm_crossover...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create the first parent */
	gprcm_init(&parent1,
			   rows, columns, sensors, actuators,
			   connections_per_gene, modules,
			   &random_seed);

	assert((&parent1)->program.genome[0].gene != 0);
	assert((&parent1)->program.genome[0].state != 0);
	assert((&parent1)->program.genome[0].used != 0);

	assert((&parent1)->morphology.genome[0].gene != 0);
	assert((&parent1)->morphology.genome[0].state != 0);
	assert((&parent1)->morphology.genome[0].used != 0);

	assert((&parent1)->program.genome[1].gene != 0);
	assert((&parent1)->program.genome[1].state != 0);
	assert((&parent1)->program.genome[1].used != 0);

	/* randomize the first parent */
	gprcm_random(&parent1,
				 rows, columns,
				 sensors, actuators,
				 connections_per_gene,
				 min_value, max_value,
				 integers_only, &random_seed,
				 modules,
				 instruction_set, no_of_instructions);

	/* create a second individual */
	gprcm_init(&parent2,
			   rows, columns, sensors, actuators,
			   connections_per_gene, modules,
			   &random_seed);

	assert((&parent2)->program.genome[0].gene != 0);
	assert((&parent2)->program.genome[0].state != 0);
	assert((&parent2)->program.genome[0].used != 0);

	assert((&parent2)->morphology.genome[0].gene != 0);
	assert((&parent2)->morphology.genome[0].state != 0);
	assert((&parent2)->morphology.genome[0].used != 0);

	assert((&parent2)->program.genome[1].gene != 0);
	assert((&parent2)->program.genome[1].state != 0);
	assert((&parent2)->program.genome[1].used != 0);

	/* randomize the second parent */
	gprcm_random(&parent2,
				 rows, columns,
				 sensors, actuators,
				 connections_per_gene,
				 min_value, max_value,
				 integers_only, &random_seed,
				 modules,
				 instruction_set, no_of_instructions);


	for (i = 0;
		 i < (rows*columns*
			  GPRC_GENE_SIZE(connections_per_gene)); i++) {
		parent1.program.genome[0].gene[i] = 1;
		parent1.program.genome[1].gene[i] = 1;
	}

	for (i = 0;
		 i < (rows*columns*
			  GPRC_GENE_SIZE(connections_per_gene)); i++) {
		parent2.program.genome[0].gene[i] = 2;
		parent2.program.genome[1].gene[i] = 2;
	}

	/* produce the child */
	gprcm_crossover(&parent1, &parent2,
					rows, columns,
					sensors, actuators,
					connections_per_gene,
					chromosomes,
					1, &child);

	/* check that there is some contribution from each parent */
	parent1_hits=0;
	parent2_hits=0;
	for (i = 0;
		 i < (rows*columns*
			  GPRC_GENE_SIZE(connections_per_gene)); i++) {
		if ((child.program.genome[0].gene[i]==1) ||
			(child.program.genome[1].gene[i]==1)) {
			parent1_hits++;
		}
		if ((child.program.genome[0].gene[i]==2) ||
			(child.program.genome[1].gene[i]==2)) {
			parent2_hits++;
		}
	}
	assert(parent1_hits>0);
	assert(parent2_hits>0);

	/* free memory */
	gprcm_free(&parent1);
	gprcm_free(&parent2);
	gprcm_free(&child);

	printf("Ok\n");
}

static void test_gprcm_sort()
{
	int population_size = 1000;
	int rows=10, columns=20, sensors=8, actuators=4;
	int i,connections_per_gene=2;
	int chromosomes=2;
	int modules=1;
	float min_value=-10, max_value=10;
	gprcm_population population;
	float f;
	gprcm_function * individual, * prev_individual;
	int integers_only=0;
	int instruction_set[64], no_of_instructions=0;
	unsigned int random_seed = 123;

	printf("test_gprcm_sort...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprcm_init_population(&population,
						  population_size,
						  rows, columns,
						  sensors, actuators,
						  connections_per_gene,
						  modules,
						  chromosomes,
						  min_value, max_value,
						  integers_only, &random_seed,
						  instruction_set, no_of_instructions);

	/* set random fitness values */
	for (i = 0; i < population_size; i++) {
		f = (float)(rand_num(&random_seed)%10000);
		(&population)->fitness[i] = f;
		individual = &(population.individual[i]);
		individual->program.genome[0].state[0] = f;
	}

	/* sort by fitness */
	gprcm_sort(&population);

	/* check the sort order */	
	for (i = 1; i < population_size; i++) {
		individual = &(population.individual[i]);
		prev_individual = &(population.individual[i-1]);
		assert(population.fitness[i-1] >= population.fitness[i]);
		assert(prev_individual->program.genome[0].state[0] >=
			   individual->program.genome[0].state[0]);
		assert(population.fitness[i-1] ==
			   prev_individual->program.genome[0].state[0]);
		assert(population.fitness[i] ==
			   individual->program.genome[0].state[0]);
	}	

	/* free memory */
	gprcm_free_population(&population);

	printf("Ok\n");
}

static void test_gprcm_sort_system()
{
	int population_per_island = 256;
	int rows = 10, columns = 20, sensors = 8, actuators = 4;
	int i, j, connections_per_gene=2;
	int chromosomes = 2;
	int modules = 1;
	float min_value = -10, max_value = 10;
	gprcm_system system;
	gprcm_population * population;
	float f;
	int integers_only = 0;
	int islands = 4;
	int instruction_set[64], no_of_instructions=0;
	unsigned int random_seed = 123;

	printf("test_gprcm_sort_system...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a system */
	gprcm_init_system(&system,
					  islands,
					  population_per_island,
					  rows, columns,
					  sensors, actuators,
					  connections_per_gene,
					  modules,
					  chromosomes,
					  min_value, max_value,
					  integers_only, &random_seed,
					  instruction_set, no_of_instructions);

	/* set random fitness values */
	for (i = 0; i < islands; i++) {
		f = (float)(rand_num(&random_seed)%10000);
		(&system)->fitness[i] = f;
		population = &system.island[i];
		for (j = 0; j < population->size; j++) {
			population->fitness[j] =
				(float)(rand_num(&random_seed)%10000);
		}
	}

	/* sort by fitness */
	gprcm_sort_system(&system);

	/* check the sort order */	
	for (i = 1; i < islands; i++) {
		assert(system.fitness[i-1] >= system.fitness[i]);
	}	

	/* free memory */
	gprcm_free_system(&system);

	printf("Ok\n");
}

/* A test evaluation function.
   This tests how close the output is to the equation y = 3x^2 + 2x - 5 */
static float test_evaluate_program(int time_steps,
								   gprcm_population * population,
								   int individual_index,
								   int custom_command)
{
	int t,x,i;
	float result,fitness=0, reference;
	float dropout_rate = 0.3f;
	gprcm_function * f = &population->individual[individual_index];

	/* clear the state */
	gprcm_clear_state(f,
					  population->rows, population->columns,
					  population->sensors, population->actuators);

	/* for each time step */
	for (t = 0; t < time_steps; t++) {
		x = t+1;
		/* sensor set to the current time step */
		for (i = 0; i < population->sensors; i++) {
			gprcm_set_sensor(f,i,x);
		}
		/* run the program */
		gprcm_run(f, population, dropout_rate, 0, 0);
		/* observe the output */
		result = gprcm_get_actuator(f,0,
									population->rows,
									population->columns,
									population->sensors);
		/* the target value */
		reference = (3*x*x) + (2*x) - 5;
		/* how close is the output to the target equation? */
		if (fabs(result - reference)<100) {
			fitness += 100 - fabs(result - reference);
		}
		
	}
	fitness /= (float)time_steps;
	return fitness;
}

static void test_gprcm_generation()
{
	int population_size = 512;
	int rows = 6, columns = 10, sensors = 5, actuators = 3;
	int connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int chromosomes = 2;
	int modules=1;
	float min_value = -5, max_value = 5;
	gprcm_population population;
	float elitism = 0.3f;
	float mutation_prob = 0.5f;
	int retval, i, gen, time_steps = 10;
	int integers_only;
	int use_crossover = 1;
	unsigned int random_seed = 123;
	char dot_filename[128];
	FILE * fp;
	int instruction_set[64], no_of_instructions=0;
	char * sensor_names[] = {
		"Sensor 0", "Sensor 1", "Sensor 2", "Sensor 3", "Sensor 4"
	};
	char * actuator_names[] = {
		"Actuator 0", "Actuator 1", "Actuator 2"
	};
	char command_str[256];

	printf("test_gprcm_generation...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	for (integers_only = 0; integers_only < 2; integers_only++) {
		/* create a population */
		gprcm_init_population(&population,
							  population_size,
							  rows, columns,
							  sensors, actuators,
							  connections_per_gene,
							  modules,
							  chromosomes,
							  min_value, max_value,
							  integers_only, &random_seed,
							  instruction_set, no_of_instructions);

		/* validate the individuals */
		for (i = 0; i < population_size; i++) {
			retval =
				gprcm_validate(&population.individual[i],
							   rows, columns,
							   sensors, actuators,
							   connections_per_gene,
							   integers_only,
							   instruction_set,
							   no_of_instructions);
			show_validation_message(retval);
			assert(retval==GPR_VALIDATE_OK);
		}

		for (gen = 0; gen < 20; gen++) {

			/* evaluate each individual */
			gprcm_evaluate(&population,
						   time_steps,0,
						   (*test_evaluate_program));

			/* produce the next generation */
			gprcm_generation(&population,
							 elitism,
							 mutation_prob,
							 use_crossover,
							 &random_seed,
							 instruction_set, no_of_instructions);

			/* validate the individuals */
			for (i = 0; i < population_size; i++) {
				retval = gprcm_validate(&population.individual[i],
										rows, columns,
										sensors, actuators,
										connections_per_gene,
										integers_only,
										instruction_set,
										no_of_instructions);
				show_validation_message(retval);
				assert(retval==GPR_VALIDATE_OK);

			}

		}

		assert(population.history.index == gen);

		/* validate the best individual */
		retval = gprcm_validate(gprcm_best_individual(&population),
								rows, columns,
								sensors, actuators,
								connections_per_gene,
								integers_only,
								instruction_set,
								no_of_instructions);
		show_validation_message(retval);
		assert(retval==GPR_VALIDATE_OK);

		/* save a dot file for the best individual */
		sprintf(dot_filename,"%slibgpr_fittest.dot",GPR_TEMP_DIRECTORY);
		fp = fopen(dot_filename,"w");
		assert(fp);
		gprcm_dot(gprcm_best_individual(&population), &population,
				  1, sensor_names, actuator_names,
				  fp);
		fclose(fp);
		sprintf(command_str,"rm %s",dot_filename);
		assert(system(command_str)==0);

		/* free memory */
		gprcm_free_population(&population);

	}

	printf("Ok\n");
}

static void test_gprcm_generation_system()
{
	int population_per_island = 256;
	int islands = 4;
	int migration_interval=10;
	int rows = 5, columns = 10, sensors = 1, actuators = 1;
	int connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int chromosomes=2;
	int modules=1;
	float min_value = -5, max_value = 5;
	gprcm_system sys;
	float elitism = 0.3f;
	float mutation_prob = 0.5f;
	int gen, time_steps = 10;
	/*float start_fitness = 0;*/
	int integers_only;
	int use_crossover = 1;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;
	char source_filename[256],compile_command[256];
	char binary_filename[256];
	FILE * fp;

	printf("test_gprcm_generation_system...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	for (integers_only = 0; integers_only < 2; integers_only++) {
		/* create a population */
		gprcm_init_system(&sys, islands,
						  population_per_island,
						  rows, columns,
						  sensors, actuators,
						  connections_per_gene,
						  modules,
						  chromosomes,
						  min_value, max_value,
						  integers_only, &random_seed,
						  instruction_set, no_of_instructions);

		for (gen = 0; gen < 20; gen++) {

			/* evaluate each population */
			gprcm_evaluate_system(&sys,
								  time_steps,0,
								  (*test_evaluate_program));

			/* produce the next generation */
			gprcm_generation_system(&sys,
									migration_interval,
									elitism,
									mutation_prob,
									use_crossover, &random_seed,
									instruction_set, no_of_instructions);

			/*if (gen == 0) start_fitness = gprcm_best_fitness_system(&system);*/
		}

		/*if (gprcm_best_fitness_system(&system) <= start_fitness) {
		  printf("\nIntegers_only = %d\n", integers_only);
		  }
		
		  assert(gprcm_best_fitness_system(&system) > start_fitness);*/

		/* save the best individual as a C program */
		sprintf(source_filename,"%sagent.c",GPR_TEMP_DIRECTORY);
		fp = fopen(source_filename,"w");
		assert(fp);
		gprcm_c_program(&sys,
						gprcm_best_individual_system(&sys),
						1, 0, fp);
		fclose(fp);

		/* compile the result */
		sprintf(binary_filename,"%sagent",GPR_TEMP_DIRECTORY);
		sprintf(compile_command,
				"gcc -Wall -std=c99 -pedantic -o %s %s -lm",
				binary_filename,source_filename);
		assert(system(compile_command)==0);

		/* free memory */
		gprcm_free_system(&sys);

	}

	printf("Ok\n");
}

static void test_gprcm_save_load()
{
	gprcm_population population, population2;
	gprcm_function *f1, *f2;
	int population_size = 1000;
	int rows = 20, columns = 40, sensors = 5, actuators = 5;
	int connections_per_gene = 2;
	int chromosomes=2;
	int modules=1;
	float min_value = -5, max_value = 5;
	int integers_only = 0, i, j;
	char filename[256];
	FILE * fp;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gprcm_save_load...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprcm_init_population(&population,
						  population_size,
						  rows, columns,
						  sensors, actuators,
						  connections_per_gene,
						  modules,
						  chromosomes,
						  min_value, max_value,
						  integers_only, &random_seed,
						  instruction_set, no_of_instructions);

	sprintf(filename,"%stestpopulation.dat",GPR_TEMP_DIRECTORY);

	/* save to file */
	fp  = fopen(filename,"w");
	assert(fp!=0);
	gprcm_save_population(&population,fp);
	fclose(fp);

	/* load from file */
	fp  = fopen(filename,"r");
	assert(fp!=0);
	gprcm_load_population(&population2,fp,
						  instruction_set, no_of_instructions);
	fclose(fp);

	/* check population parameters */
	assert(population.history.index==population2.history.index);
	assert(population.history.tick==population2.history.tick);
	assert(population.history.interval==population2.history.interval);
	assert(population.size==population2.size);
	assert(population.rows==population2.rows);
	assert(population.connections_per_gene==
		   population2.connections_per_gene);
	assert(population.columns==population2.columns);
	assert(population.sensors==population2.sensors);
	assert(population.actuators==population2.actuators);
	assert(population.integers_only==population2.integers_only);
	assert(abs((int)(population.min_value*1000) -
			   (int)(population2.min_value*1000)) < 2);
	assert(abs((int)(population.max_value*1000) -
			   (int)(population2.max_value*1000)) < 2);

	/* check that individuals are the same */
	for (i = 0; i < population.size; i++) {
		f1 = &population.individual[i];
		f2 = &population2.individual[i];
		for (j = 0;
			 j < (rows*columns*GPRC_GENE_SIZE(connections_per_gene))+
				 actuators; j++) {
			assert(f1->program.genome[0].gene[j] ==
				   f2->program.genome[0].gene[j]);
		}
		assert(f1->program.random_seed == f2->program.random_seed);
		for (j = 0;
			 j < (rows*columns) + sensors + actuators; j++) {
			if (f1->program.genome[0].used[j] !=
				f2->program.genome[0].used[j]) {
				printf("\n%d  %d %d\n",
					   j,(int)f1->program.genome[0].used[j],
					   (int)f2->program.genome[0].used[j]);
			}
			assert(f1->program.genome[0].used[j] ==
				   f2->program.genome[0].used[j]);
		}
	}

	gprcm_free_population(&population);
	gprcm_free_population(&population2);

	printf("Ok\n");
}

static void test_gprcm_save_load_system()
{
	int islands=4;
	gprcm_system system, system2;
	gprcm_population *p1, *p2;
	gprcm_function *f1, *f2;
	int population_per_island = 256;
	int rows = 20, columns = 40, sensors = 5, actuators = 5;
	int connections_per_gene = 2;
	int chromosomes=2;
	int m,modules=1;
	float min_value = -5, max_value = 5;
	int integers_only = 0, i, j, k;
	char filename[256];
	FILE * fp;
	unsigned int random_seed = 123;
	int no_of_sensor_sources = 120;
	int no_of_actuator_destinations = 64;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gprcm_save_load_system...");

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprcm_init_system(&system,islands,
					  population_per_island,
					  rows, columns,
					  sensors, actuators,
					  connections_per_gene,
					  modules,
					  chromosomes,
					  min_value, max_value,
					  integers_only, &random_seed,
					  instruction_set, no_of_instructions);

	gprcm_init_sensor_sources(&system,
							  no_of_sensor_sources,
							  &random_seed);
	gprcm_init_actuator_destinations(&system,
									 no_of_actuator_destinations,
									 &random_seed);
	system.migration_tick=2;

	sprintf(filename,"%stestsystem.dat",GPR_TEMP_DIRECTORY);

	/* save to file */
	fp  = fopen(filename,"w");
	assert(fp!=0);
	gprcm_save_system(&system,fp);
	fclose(fp);

	/* load from file */
	fp  = fopen(filename,"r");
	assert(fp!=0);
	gprcm_load_system(&system2,fp,
					  instruction_set, no_of_instructions);
	fclose(fp);

	/* check system parameters */
	assert(system.size==system2.size);
	assert(system.migration_tick==system2.migration_tick);

	/* check that populations are the same */
	for (i = 0; i < islands; i++) {
		p1 = &system.island[i];
		p2 = &system2.island[i];

		for (k = 0; k < p1->size; k++) {
			f1 = &p1->individual[k];
			f2 = &p2->individual[k];

			for (m = 0; m < modules+1; m++) {
				for (j = 0;
					 j < (rows*columns*
						  GPRC_GENE_SIZE(connections_per_gene))+
						 gprcm_get_actuators(m,actuators); j++) {
					assert(f1->program.genome[m].gene[j] ==
						   f2->program.genome[m].gene[j]);
				}
			}
			assert(f1->program.no_of_sensor_sources ==
				   f2->program.no_of_sensor_sources);
			assert(f1->program.no_of_actuator_destinations ==
				   f2->program.no_of_actuator_destinations);
			assert(f1->program.random_seed ==
				   f2->program.random_seed);

			for (m = 0; m < modules+1; m++) {
				for (j = 0;
					 j < (rows*columns) +
						 gprcm_get_sensors(m,sensors) +
						 gprcm_get_actuators(m,actuators); j++) {
					if (f1->program.genome[m].used[j] !=
						f2->program.genome[m].used[j]) {
						printf("\n%d  %d %d\n",
							   j,(int)f1->program.genome[m].used[j],
							   (int)f2->program.genome[m].used[j]);
					}
					assert(f1->program.genome[m].used[j] ==
						   f2->program.genome[m].used[j]);
				}
			}

			for (j = 0; j < sensors; j++) {
				assert(f1->program.sensor_source[j] ==
					   f2->program.sensor_source[j]);
				assert(f1->program.sensor_source[j] >= 0);
				assert(f1->program.sensor_source[j] <
					   no_of_sensor_sources);
				assert(f2->program.sensor_source[j] >= 0);
				assert(f2->program.sensor_source[j] <
					   no_of_sensor_sources);
			}
			for (j = 0; j < actuators; j++) {
				assert(f1->program.actuator_destination[j] ==
					   f2->program.actuator_destination[j]);
				assert(f1->program.actuator_destination[j] >= 0);
				assert(f1->program.actuator_destination[j] <
					   no_of_actuator_destinations);
				assert(f2->program.actuator_destination[j] >= 0);
				assert(f2->program.actuator_destination[j] <
					   no_of_actuator_destinations);
			}
		}
	}

	gprcm_free_system(&system);
	gprcm_free_system(&system2);

	printf("Ok\n");
}

static void set_function(gprcm_function * f,
						 int ADF_module,
						 int row, int col,						 
						 int function_type,
						 float value,
						 int connection_sensor,
						 int connection1_row, int connection1_col,
						 int connection2_row, int connection2_col,
						 int rows, int columns,
						 int connections_per_gene,
						 int sensors)
{
	int ctr = GPRC_INITIAL;
	int sens = gprcm_get_sensors(ADF_module,sensors);
	int index = sens + ((col*rows) + row);
	int n = ((col*rows) + row)*GPRC_GENE_SIZE(connections_per_gene);

	f->program.genome[ADF_module].used[index] = 1;
	f->program.genome[ADF_module].gene[n] = function_type;
	f->program.genome[ADF_module].gene[n+1] = value;
	if (connection_sensor>-1) {
		f->program.genome[ADF_module].used[connection_sensor] = 1;
		f->program.genome[ADF_module].gene[n+ctr] = connection_sensor;
		ctr++;
	}
	if (connection1_col>-1) {
		f->program.genome[ADF_module].gene[n+ctr] =
			sens + ((connection1_col*rows) + connection1_row);
		ctr++;
	}
	if (connection2_col>-1) {
		f->program.genome[ADF_module].gene[n+ctr] =
			sens + ((connection2_col*rows) + connection2_row);
		ctr++;
	}
}

static void test_gprcm_compress_ADF()
{
	gprcm_function f;
	gprcm_population population;
	int i,itt,rows=9, columns=20, sensors=2, actuators=1;	
	int connections_per_gene=10,index,n;
	int chromosomes=3;
	int no_of_genes, no_of_inputs, retval, modules = 2;
	int start_row=1;
	int start_col=columns-1;
	int instruction_set[64], no_of_instructions=0;
	float value_before, value_after, min_value = -1;
	float max_value = 1;
	int integers_only=0;
	unsigned int random_seed = 123;
	char before_filename[256];
	char after_filename[256];
	FILE * fp;
	char * sensor_names[] = {
		"Sensor 0", "Sensor 1"
	};
	char * actuator_names[] = {
		"Actuator 0"
	};

	printf("test_gprcm_compress_ADF...");	

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprcm_init_population(&population,
						  2,
						  rows, columns,
						  sensors, actuators,
						  connections_per_gene,
						  modules,
						  chromosomes,
						  min_value, max_value,
						  integers_only, &random_seed,
						  instruction_set, no_of_instructions);

	f = (&population)->individual[0];

	/* check that arrays are not null */
	assert((&f)->program.genome[0].gene != 0);
	assert((&f)->program.genome[0].state != 0);
	assert((&f)->program.genome[0].used != 0);
	assert((&f)->program.temp_genes != 0);

	assert((&f)->morphology.genome[0].gene != 0);
	assert((&f)->morphology.genome[0].state != 0);
	assert((&f)->morphology.genome[0].used != 0);
	assert((&f)->morphology.temp_genes != 0);

	assert((&f)->program.genome[1].gene != 0);
	assert((&f)->program.genome[1].state != 0);
	assert((&f)->program.genome[1].used != 0);
	assert((&f)->program.temp_genes != 0);

	for (index = 0; index < gprcm_get_sensors(0,sensors); index++) {
		(&f)->program.genome[0].used[index] = 0;
	}
	n = 0;
	for (index = 0; index < rows*columns; index++,
			 n += GPRC_GENE_SIZE(connections_per_gene)) {
		(&f)->program.genome[0].gene[n] = GPR_FUNCTION_VALUE;
		(&f)->program.genome[0].used[index] = 0;
	}

	/* make a simple program */
	set_function(&f, 0, start_row, start_col,
				 GPR_FUNCTION_MIN,1,
				 -1,
				 start_row, start_col-1,
				 start_row+1, start_col-1,
				 rows, columns, connections_per_gene,
				 sensors);

	set_function(&f, 0, start_row, start_col-1,
				 GPR_FUNCTION_MULTIPLY,0,
				 -1,
				 start_row, start_col-2,
				 -1, -1,
				 rows, columns, connections_per_gene,
				 sensors);

	set_function(&f, 0, start_row+1, start_col-1,
				 GPR_FUNCTION_ADD,0,
				 -1,
				 start_row+1, start_col-2,
				 -1, -1,
				 rows, columns, connections_per_gene,
				 sensors);

	set_function(&f, 0, start_row, start_col-2,
				 GPR_FUNCTION_SUBTRACT,1,
				 0,
				 start_row, start_col-3,
				 -1, -1,
				 rows, columns, connections_per_gene,
				 sensors);

	set_function(&f, 0, start_row, start_col-3,
				 GPR_FUNCTION_VALUE,7.3f,
				 0,
				 -1, -1,
				 -1, -1,
				 rows, columns, connections_per_gene,
				 sensors);

	set_function(&f, 0, start_row+1, start_col-2,
				 GPR_FUNCTION_NOOP1,0,
				 1,
				 -1, -1,
				 -1, -1,
				 rows, columns, connections_per_gene,
				 sensors);

	f.program.genome[0].gene[rows*columns*
					 GPRC_GENE_SIZE(connections_per_gene)] =
		sensors + ((start_col*rows) + start_row);
	f.program.genome[0].used[sensors + (rows*columns)] = 1;

	sprintf(before_filename,"%stemp_before.dot",GPR_TEMP_DIRECTORY);
	sprintf(after_filename,"%stemp_after.dot",GPR_TEMP_DIRECTORY);
	fp = fopen(before_filename,"w");
	assert(fp);
	gprcm_dot(&f, &population, 1,
			  sensor_names, actuator_names, fp);
	fclose(fp);


	no_of_genes=0;
	no_of_inputs=0;
	gprcm_get_subgraph(&f, 0, sensors + ((start_col*rows)+start_row),
					   0, 0, rows, columns,
					   connections_per_gene,
					   sensors, 0,
					   10, GPRC_MAX_ADF_GENES,
					   (&f)->program.temp_genes, &no_of_genes,
					   &no_of_inputs, 0);
	if ((no_of_genes!=8) || (no_of_inputs!=2)) {
		printf("\nno_of_genes %d\nno_of_inputs %d\n",
			   no_of_genes, no_of_inputs);
	}
	assert(no_of_genes==8);
	assert(no_of_inputs==2);

	/* run the program */
	for (itt = 0; itt < 4; itt++) {
		for (i = 0; i < sensors; i++) {
			gprcm_set_sensor(&f,i,i*2);
		}
		gprcm_run(&f, &population, 0, 0, 0);
	}
	value_before = gprcm_get_actuator(&f, 0,
									  rows, columns, sensors);

	/* compress into an ADF */
	retval = gprcm_compress_ADF(&f, 0,
								sensors + ((start_col*rows)+start_row),
								rows, columns,
								connections_per_gene,
								sensors, actuators,
								min_value, max_value,
								10,0);
	switch(retval) {
	case -1: {
		printf("\nNo ADF modules\n");
		break;
	}
	case -2: {
		printf("\nNo used genes were found\n");
		break;
	}
	case -3: {
		printf("\nNumber of genes less than minimum\n");
		break;
	}
	case -4: {
		printf("\nNumber of genes greater than maximum\n");
		break;
	}
	case -5: {
		printf("\nNumber of sensors greater than maximum\n");
		break;
	}
	case -6: {
		printf("\nNumber of sensors is zero\n");
		break;
	}
	}
    assert(retval>-1);	

	fp = fopen(after_filename,"w");
	assert(fp);
	gprcm_dot(&f, &population, 1,
			  sensor_names, actuator_names, fp);
	fclose(fp);

	/* count the number of used genes */
	n=0;
	for (index = 0; index < rows*columns; index++) {
		if ((&f)->program.genome[1].used[index]!=0) n++;
	}
	if (n != 7) {
		printf("Number of used genes in the module = %d\n",n);
	}
	assert(n == 7);

	/* run the program again */
	for (itt = 0;itt < 4; itt++) {
		for (i = 0; i < sensors; i++) {
			gprcm_set_sensor(&f,i,i*2);
		}
		gprcm_run(&f, &population, 0, 0, 0);
	}
	value_after = gprcm_get_actuator(&f, 0,
									 rows, columns, sensors);

	if (fabs(value_before-value_after) >= 0.01f) {
		printf("before %.2f  after %.2f\n",
			   value_before, value_after);
	}
	assert(fabs(value_before-value_after) < 0.01f);

	/* free memory */
	gprcm_free_population(&population);

	printf("Ok\n");
}

int run_tests_morphology()
{
	printf("\nRunning Morphological Cartesian tests\n");

	test_gprcm_init();
	test_gprcm_random();
	test_gprcm_copy();
	test_gprcm_run();
	test_gprcm_run_dynamic();
	test_gprcm_mutate();
	test_gprcm_crossover();
	test_gprcm_sort();
	test_gprcm_sort_system();
	test_gprcm_generation();
	test_gprcm_generation_system();
	test_gprcm_save_load();
	test_gprcm_save_load_system();
	test_gprcm_compress_ADF();
	test_gprcm_environment();

	printf("All Morphological Cartesian tests completed\n");
	return 1;
}

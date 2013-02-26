/*
  Predicting concrete slump from various ingredients
  See https://en.wikipedia.org/wiki/Concrete_slump_test
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

#include <stdio.h>
#include <time.h>
#include "libgpr/globals.h"
#include "libgpr/gprcm.h"

#define MAX_EXAMPLES 200
#define MAX_TEST_EXAMPLES 20
#define MAX_FIELDS   12

#define RUN_STEPS  2

/* the required measurement accuracy in cm */
#define REQUIRED_ACCURACY_CM       1.0

/* required accuracy for 28 day hardness */
#define REQUIRED_ACCURACY_HARDNESS 2.0

float concrete_data[MAX_EXAMPLES*MAX_FIELDS];
float test_data[MAX_TEST_EXAMPLES*MAX_FIELDS];
float * current_data_set;

int no_of_examples = 0;
int fields_per_example = 0;

/* create a test data set from the original data.
   The test data can be used to calculate a final fitness
   value, because it was not seen during training and so
   provides an indication of how well the system has generalised */
static int create_test_data(float * training_data,
							int * no_of_training_examples,
							int fields_per_example,
							float * test_data)
{
	int i,j,k,index;
	int no_of_test_examples = 0;
	unsigned int random_seed = (unsigned int)time(NULL);

	for (i = 0; i < MAX_TEST_EXAMPLES; i++) {
		/* pick an example from the loaded data set */
		index = rand_num(&random_seed)%(*no_of_training_examples);

		/* increase the number of test examples */
		for (j = 0; j < fields_per_example; j++) {
			test_data[no_of_test_examples*fields_per_example + j] =
				training_data[index*fields_per_example + j];
		}
		no_of_test_examples++;

		/* reshuffle the original data set */
		for (j = index+1; j < (*no_of_training_examples); j++) {
			for (k = 0; k < fields_per_example; k++) {
				training_data[(j-1)*fields_per_example + k] = 
					training_data[j*fields_per_example + k];
			}
		}
		/* decrease the number of training data examples */
		*no_of_training_examples = *no_of_training_examples - 1;
	}

	return no_of_test_examples;
}
							

static int load_data(char * filename, float * training_data,
					 int max_examples,
					 int * fields_per_example)
{
	int i, field_number, ctr, examples_loaded = 0;
	FILE * fp;
	char line[2000],valuestr[256],*retval;
	float value;
	int training_data_index = 0;

	fp = fopen(filename,"r");
	if (!fp) return 0;

	while (!feof(fp)) {
		retval = fgets(line,1999,fp);
		if (retval) {
			if (strlen(line)>0) {
				field_number = 0;
				ctr = 0;
				for (i = 0; i < strlen(line); i++) {
					if ((line[i]==',') ||
						(i==strlen(line)-1)) {
						if (i==strlen(line)-1) {
							valuestr[ctr++]=line[i];
						}
						valuestr[ctr]=0;
						ctr=0;

						/* get the value from the string */
						value = 0;
						if (valuestr[0]!='?') {
							if ((valuestr[0]>='0') &&
								(valuestr[0]<='9')) {
								value = atof(valuestr);
							}
						}

						/* insert value into the array */
						training_data[training_data_index] = value;
						field_number++;
						training_data_index++;
					}
					else {
						/* update the value string */
						valuestr[ctr++]=line[i];
					}
				}
				*fields_per_example = field_number;
				examples_loaded++;
				if (examples_loaded >= max_examples) {
					fclose(fp);
					return examples_loaded;
				}
			}
		}
	}

	fclose(fp);

	return examples_loaded;
}

static float evaluate_features(int trials,
							   gprcm_population * population,
							   int individual_index,
							   int custom_command)
{
	int i,j,itt;
	float error,diff=0,v,reference,fitness;
	float dropout_rate = 0.0f;
	gprcm_function * f = &population->individual[individual_index];

	if (custom_command!=0) dropout_rate=0;

	for (i = 0; i < trials; i++) {
		/* clear the state */
		gprcm_clear_state(f,
						  population->rows, population->columns,
						  population->sensors, population->actuators);

		for (j=1;j<fields_per_example-3;j++) {
			gprcm_set_sensor(f,j-1,
							 current_data_set[i*fields_per_example+j]);
		}
		for (itt = 0; itt < RUN_STEPS; itt++) {
			/* run the program */
			gprcm_run(f, population, dropout_rate, 0, 0);
		}
		/* how close is the output to the actual slump? */
		for (j=0;j<3;j++) {
			reference =
				0.01f + fabs(current_data_set[i*fields_per_example+
											  fields_per_example-3+j]);
			v = 0.01f + fabs(gprcm_get_actuator(f,j,
												population->rows,
												population->columns,
												population->sensors));

			error = fabs(v - reference)/reference;
			diff += error*error;
		}
	}
	diff = (float)sqrt(diff/(float)(trials*3))*100;
	fitness = 100 - diff;
	if (fitness<0) fitness=0;
	return fitness;
}

static void slump_test()
{
	int islands = 4;
	int migration_interval = 500;
	int population_per_island = 64;
	int rows = 9, columns = 16;
	int i, gen=0;
	int connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int modules = 0;
	int chromosomes=3;
	gprcm_system sys;
	float min_value = -100;
	float max_value = 100;
	float elitism = 0.2f;
	float mutation_prob = 0.1f;
	int trials = 100;
	int use_crossover = 1;
	unsigned int random_seed = (unsigned int)time(NULL);
	int sensors=4, actuators=3;
	int integers_only = 0;
	int no_of_test_examples;
	float test_performance;
	FILE *fp;
	char compile_command[256];
	int instruction_set[64], no_of_instructions=0;
	char * sensor_names[] = {
		"Cement",
		"Slag",
		"Fly Ash",
		"Water",
		"Sand",
		"Coarse Aggr.",
		"Fine Aggr."
	};
	char * actuator_names[] = { "Slump", "Flow", "Compressive Strength" };

	/* load the data */
	no_of_examples =
		load_data("slump_test.data", concrete_data, MAX_EXAMPLES,
				  &fields_per_example);

	/* create a test data set */
	no_of_test_examples = create_test_data(concrete_data, &no_of_examples,
										   fields_per_example,
										   test_data);

	sensors = fields_per_example-4;
	trials = no_of_examples;

	printf("Number of training examples: %d\n",no_of_examples);
	printf("Number of test examples: %d\n",no_of_test_examples);
	printf("Number of fields: %d\n",fields_per_example);

	/* create an instruction set */
	no_of_instructions = gprcm_dynamic_instruction_set(instruction_set);

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

	gpr_xmlrpc_server("server.rb","concreteslump",3573,
					  "./agent",
					  sensors, actuators);

	test_performance = 0;
	while (test_performance < 99) {
		/* use the training data */
		current_data_set = concrete_data;

		/* evaluate each individual */
		gprcm_evaluate_system(&sys,
							  trials,0,
							  (*evaluate_features));
		/* produce the next generation */
		gprcm_generation_system(&sys,
								migration_interval,
								elitism,
								mutation_prob,
								use_crossover, &random_seed,
								instruction_set, no_of_instructions);

		/* evaluate the test data set */
		current_data_set = test_data;
		test_performance = evaluate_features(no_of_test_examples,
											 &sys.island[0],
											 0,1);

		/* show the best fitness value calculated from the test data set */
		printf("Generation %05d  Fitness %.2f/%.2f%% ",
			   gen, gprcm_best_fitness(&sys.island[0]),test_performance);
		for (i = 0; i < islands; i++) {
			printf("  %.3f",gprcm_average_fitness(&sys.island[i]));
		}
		printf("\n");

		if (((gen % 100 == 0) && (gen>0)) || (test_performance > 99)) {
			gprcm_draw_population("population.png",
								  640, 640, &sys.island[0]);

			gprcm_plot_history_system(&sys,
									  GPR_HISTORY_FITNESS,
									  "fitness.png",
									  "Concrete Slump Performance",
									  640, 480);

			gprcm_plot_history_system(&sys,
									  GPR_HISTORY_AVERAGE,
									  "fitness_average.png",
									  "Concrete Slump Average Performance",
									  640, 480);

			gprcm_plot_history_system(&sys,
									  GPR_HISTORY_DIVERSITY,
									  "diversity.png",
									  "Concrete Slump Diversity",
									  640, 480);

			gprcm_plot_fitness(&sys.island[0],
							   "fitness_histogram.png",
							   "Fitness Histogram",
							   640, 480);

			fp = fopen("agent.c","w");
			if (fp) {
				/* save the best program */
				gprcm_c_program(&sys,
								gprcm_best_individual_system(&sys),
								RUN_STEPS, 0, fp);
				fclose(fp);

				/* compile the program */
				sprintf(compile_command,
						"gcc -Wall -std=c99 -pedantic -o agent agent.c -lm");
				assert(system(compile_command)==0);
			}

			fp = fopen("fittest.dot","w");
			if (fp) {
				gprcm_dot(gprcm_best_individual_system(&sys),
						  &sys.island[0], 1,
						  sensor_names,  actuator_names,
						  fp);
				fclose(fp);
			}
		}

		if (test_performance > 99) break;
		gen++;
	}

	/* free memory */
	gprcm_free_system(&sys);
}

int main(int argc, char* argv[])
{	
	slump_test();
	return 1;
}


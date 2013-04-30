/*
  Leaf image classification

  This example classifies a particular species of leaf out of a total
  of 100 species based upon 64 shape, texture and margin features.
  If you wanted to devise a system capable of classifying multiple
  species then repeat the process for different species and
  run the resulting best programs in parallel, with the highest
  output value winning.

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

#define MAX_EXAMPLES      1600
#define MAX_TEST_EXAMPLES 600
#define MAX_FIELDS        (1+64)
#define SAMPLES_PER_CLASS 16
#define SPECIES           100

#define INITIAL_FIELDS    1

#define RUN_STEPS         2

float shape_feature_data[MAX_EXAMPLES*MAX_FIELDS];
float margin_feature_data[MAX_EXAMPLES*MAX_FIELDS];
float texture_feature_data[MAX_EXAMPLES*MAX_FIELDS];
float shape_test_data[MAX_TEST_EXAMPLES*MAX_FIELDS];
float margin_test_data[MAX_TEST_EXAMPLES*MAX_FIELDS];
float texture_test_data[MAX_TEST_EXAMPLES*MAX_FIELDS];
float * current_shape_data_set;
float * current_margin_data_set;
float * current_texture_data_set;

int no_of_shape_examples = 0;
int no_of_margin_examples = 0;
int no_of_texture_examples = 0;

int shape_fields_per_example = 0;
int margin_fields_per_example = 0;
int texture_fields_per_example = 0;

/* index of the current species being trained, in the range 0-99 */
int species_index = 1;

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
	int classification = -1;
	unsigned int random_seed = (unsigned int)time(NULL);

	for (i = 0; i < MAX_TEST_EXAMPLES; i++) {
		
		/* ensure that both positive and negative examples exist in the
		   test set */
		if (i < SAMPLES_PER_CLASS/2) { /* half the species examples */
			while (classification != species_index) {
				index = rand_num(&random_seed)%(*no_of_training_examples);
				classification =
					(int)training_data[index*fields_per_example];
			}
		}
		else {
			while (classification == species_index) {
				index = rand_num(&random_seed)%(*no_of_training_examples);
				classification =
					(int)training_data[index*fields_per_example];
			}
		}

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

static int load_data(char * filename,
					 float * training_data,
					 int max_examples,
					 int * fields_per_example)
{
	int i, field_number, ctr, examples_loaded = 0;
	FILE * fp;
	char line[2000], valuestr[256], *retval;
	float value;
	int training_data_index = 0;

	fp = fopen(filename,"r");
	if (!fp) return 0;

	while (!feof(fp)) {
		retval = fgets(line,1999,fp);
		if (retval) {
			if (strlen(line) > 0) {
				field_number = 0;
				ctr = 0;
				for (i = 0; i < strlen(line); i++) {
					if ((line[i] == ',') ||
						(i == strlen(line)-1)) {
						if (i == strlen(line)-1) {
							valuestr[ctr++] = line[i];
						}
						valuestr[ctr] = 0;
						ctr = 0;

						/* get the value from the string */
						value = 0;
						if (valuestr[0] != '?') {
							if ((valuestr[0] >= '0') &&
								(valuestr[0] <= '9')) {
								value = atof(valuestr);
							}
						}

						if (training_data_index%MAX_FIELDS == 0) {
							value = (int)(training_data_index /
										  (MAX_FIELDS*
										   SAMPLES_PER_CLASS));
						}

						/* insert value into the array */
						training_data[training_data_index] = value;
						field_number++;
						training_data_index++;
					}
					else {
						/* update the value string */
						valuestr[ctr++] = line[i];
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
							   int mode)
{
	int i,j,itt,n;
	float positive_examples=0, negative_examples=0;
	float diff_positive=0,diff_negative=0;
	float v,fitness,classification;
	float dropout_rate = 0.0f;
	gprcm_function * f = &population->individual[individual_index];
	int fields_per_example = shape_fields_per_example;

	if (mode!=0) dropout_rate=0;

	for (i = 0; i < trials; i++) {
		/* clear the state */
		gprcm_clear_state(f,
						  population->rows, population->columns,
						  population->sensors, population->actuators);

		/* randomly choose and example */
		n = i; /*rand_num(&f->program.random_seed)%trials;*/

		/* set the sensor values */
		for (j = INITIAL_FIELDS; j < fields_per_example; j++) {
			gprcm_set_sensor(f, j - INITIAL_FIELDS,
							 current_shape_data_set[n*fields_per_example+j]);
			gprcm_set_sensor(f, j - INITIAL_FIELDS + 64,
							 current_margin_data_set[n*fields_per_example+j]);
			gprcm_set_sensor(f, j - INITIAL_FIELDS + 128,
							 current_texture_data_set[n*fields_per_example+j]);
		}

		for (itt = 0; itt < RUN_STEPS; itt++) {
			/* run the program */
			gprcm_run(f, population, dropout_rate, 0, 0);
		}

		/* how close is the output to the actual crime level? */
		classification =
			(int)current_shape_data_set[n*fields_per_example];

		v = gprcm_get_actuator(f, 0, population->rows,
							   population->columns,
							   population->sensors);

		if (species_index != classification) {
			if (v > -0.1f) {
				diff_negative+=1+(v*v);
			}
			negative_examples+=1+(v*v);
		}
		else {
			if (v < 0.1f) {
				diff_positive+=1+(v*v);
			}
			positive_examples+=1+(v*v);
		}
	}
	fitness = 100;
	if (positive_examples > 0) {
		fitness -= (diff_positive*50/(float)positive_examples);	
	}
	else {
		fitness -= 50;
	}
	if (negative_examples > 0) {
		fitness -= (diff_negative*50/(float)negative_examples);
	}
	else {
		fitness -= 50;
	}

	if (fitness < 0) fitness=0;
	return fitness;
}

static void leaf_classification()
{
	int islands = 2;
	int migration_interval = 50;
	int population_per_island = 128;
	int rows = 6, columns = 12;
	int i, gen=0;
	int connections_per_gene = 8;
	int chromosomes = 1;
	gprcm_system sys;
	float min_value = -100;
	float max_value = 100;
	float elitism = 0.4f;
	float mutation_prob = 0.3f;
	float max_fitness = 98;
	int trials = 100;
	int use_crossover = 1;
	int ADF_modules = 0;
	unsigned int random_seed = (unsigned int)time(NULL);
	int sensors=4, actuators=1;
	int integers_only = 0;
	int no_of_shape_test_examples;
	int no_of_margin_test_examples;
	int no_of_texture_test_examples;
	float test_performance;
	FILE * fp;	
	char compile_command[256];
	int instruction_set[64], no_of_instructions=0;
	int data_size=0, data_fields=0;
	char * sensor_names[64*3];
	char * actuator_names[] = {
		"Classification"
	};

	for (i = 0; i < 64*3; i++) {
		sensor_names[i] = (char*)malloc(12);
		sprintf(sensor_names[i],"Feature %d",i);
	}

	/* load the data */
	no_of_shape_examples =
		load_data("data_Sha_64.txt",
				  shape_feature_data, MAX_EXAMPLES,
				  &shape_fields_per_example);
	no_of_margin_examples =
		load_data("data_Mar_64.txt",
				  margin_feature_data, MAX_EXAMPLES,
				  &margin_fields_per_example);
	no_of_texture_examples =
		load_data("data_Tex_64.txt",
				  texture_feature_data, MAX_EXAMPLES,
				  &texture_fields_per_example);

	/* create a test data set */
	no_of_shape_test_examples =
		create_test_data(shape_feature_data,
						 &no_of_shape_examples,
						 shape_fields_per_example,
						 shape_test_data);
	no_of_margin_test_examples =
		create_test_data(margin_feature_data,
						 &no_of_margin_examples,
						 margin_fields_per_example,
						 margin_test_data);
	no_of_texture_test_examples =
		create_test_data(texture_feature_data,
						 &no_of_texture_examples,
						 texture_fields_per_example,
						 texture_test_data);

	sensors = 64 * 3;;
	trials = no_of_shape_examples;

	printf("Number of training examples: %d\n",no_of_shape_examples);
	printf("Number of test examples: %d %d %d\n",
		   no_of_shape_test_examples,
		   no_of_margin_test_examples,
		   no_of_texture_test_examples);
	printf("Number of fields: %d %d %d\n",
		   shape_fields_per_example,
		   margin_fields_per_example,
		   texture_fields_per_example);

	/* create an instruction set */
	no_of_instructions =
		gprcm_default_instruction_set((int*)instruction_set);

	/* create a population */
	gprcm_init_system(&sys, islands,
					  population_per_island,
					  rows, columns,
					  sensors, actuators,
					  connections_per_gene,
					  ADF_modules,
					  chromosomes,
					  min_value, max_value,
					  integers_only,
					  data_size, data_fields,
					  &random_seed,
					  instruction_set, no_of_instructions);

	gpr_xmlrpc_server("server.rb","crime",3573,
					  "./agent",
					  sensors, actuators);
	
	test_performance = 0;
	while (test_performance < max_fitness) {
		/* use the training data */
		current_shape_data_set = shape_feature_data;
		current_margin_data_set = margin_feature_data;
		current_texture_data_set = texture_feature_data;

		/* evaluate each individual */
		gprcm_evaluate_system(&sys, trials,0,
							  (*evaluate_features));

		/* produce the next generation */
		gprcm_generation_system(&sys, migration_interval,
								elitism, mutation_prob,
								use_crossover, &random_seed,
								instruction_set,
								no_of_instructions);

		/* evaluate the test data set */
		current_shape_data_set = shape_test_data;
		current_margin_data_set = margin_test_data;
		current_texture_data_set = texture_test_data;
		test_performance =
			evaluate_features(no_of_shape_test_examples,
							  &sys.island[0], 0, 1);

		/* show the best fitness value */
		printf("Generation %05d  Fitness %.3f/%.3f%% ",
			   gen, gprcm_best_fitness(&sys.island[0]),test_performance);
		for (i = 0; i < islands; i++) {
			printf("  %.5f",gprcm_average_fitness(&sys.island[i]));
		}
		printf("\n");

		test_performance =
			(test_performance+gprcm_best_fitness(&sys.island[0]))*0.5f;

		if (((gen % 50 == 0) && (gen>0)) ||
			(test_performance > max_fitness)) {
			gprcm_draw_population("population.png",
								  640, 640, &sys.island[0]);

			gprcm_plot_history_system(&sys,
									  GPR_HISTORY_FITNESS,
									  "fitness.png",
									  "Leaf Shape Classification " \
									  "Performance",
									  640, 480);

			gprcm_plot_history_system(&sys,
									  GPR_HISTORY_AVERAGE,
									  "fitness_average.png",
									  "Leaf Shape Classification " \
									  "Average Performance",
									  640, 480);

			gprcm_plot_history_system(&sys,
									  GPR_HISTORY_DIVERSITY,
									  "diversity.png",
									  "Leaf Shape Classification " \
									  "Diversity",
									  640, 480);

			gprcm_plot_fitness(&sys.island[0],
							   "fitness_histogram.png",
							   "Leaf Shape Classification " \
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
						"gcc -Wall -std=c99 -pedantic " \
						"-o species%02d agent.c -lm", species_index);
				assert(system(compile_command)==0);
			}

			fp = fopen("fittest.dot","w");
			if (fp) {
				gprcm_dot(gprcm_best_individual_system(&sys),
						  &sys.island[0],
						  sensor_names,  actuator_names,
						  fp);
				fclose(fp);
			}
		}

		if (test_performance > max_fitness) break;
		gen++;
	}

	for (i = 0; i < MAX_FIELDS; i++) {
		free(sensor_names[i]);
	}

	/* free memory */
	gprcm_free_system(&sys);
}

int main(int argc, char* argv[])
{	
	if (argc>1) {
		/* specify the species as a command line argument */
		species_index = atoi(argv[1]);
		/* limit within range */
		if (species_index < 0) species_index = 0;
		if (species_index >= SPECIES) species_index = SPECIES-1;
		printf("Species index %d\n", species_index);
	}
	leaf_classification();
	return 1;
}

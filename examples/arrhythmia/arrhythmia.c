/*
 Classifying heart arrhythmias into 16 categories
 Copyright (C) 2012  Bob Mottram <bob@sluggish.dyndns.org>

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
#include "libgpr/gprc.h"

#define MAX_EXAMPLES 500
#define MAX_TEST_EXAMPLES 100
#define MAX_FIELDS   280

#define CATEGORIES 16

#define RUN_STEPS  1

float arrhythmia_data[MAX_EXAMPLES*MAX_FIELDS];
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
							value = atof(valuestr);
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
				assert(value>=0);
				assert(value<=16);
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
							   gprc_population * population,
							   int individual_index,
							   int custom_command)
{
	int i,j,classification,itt,c;
	float diff=0,fitness,v,max;
	float dropout_rate = 0.2f;
	gprc_function * f = &population->individual[individual_index];

	/* include all genes when testing */
	if (custom_command != 0) dropout_rate=0;

	for (i = 0; i < trials; i++) {
		/* clear the state */
		gprc_clear_state(f,
						 population->rows, population->columns,
						 population->sensors, population->actuators);

		for (j = 0; j < fields_per_example - 1; j++) {
			gprc_set_sensor(f, j,
							current_data_set[i*fields_per_example+j]);
		}

		for (itt = 0; itt < RUN_STEPS; itt++) {
			/* run the program */
			gprc_run(f, population, dropout_rate, 0, 0);
		}

		classification =
			(int)current_data_set[i*fields_per_example+
								  fields_per_example-1];

		max = -9999;
		c = -1;
		for (j = 1; j <= CATEGORIES; j++) {
			v = gprc_get_actuator(f, j,
								  population->rows,
								  population->columns,
								  population->sensors);			
			if ((v > max) || (j == 1)) {
				max = v;
				c = j;
			}
		}
		if (c != classification) diff += 1;
	}
	diff = (diff / trials)*100;
	fitness = 100 - diff;
	return fitness;
}

static void arrhythmia_classification()
{
	int islands = 4;
	int migration_interval = 200;
	int population_per_island = 64;
	int rows = 9, columns = 16;
	int i, gen=0;
	int connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int modules = 0;
	int chromosomes=3;
	gprc_system sys;
	float min_value = -100;
	float max_value = 100;
	float elitism = 0.2f;
	float mutation_prob = 0.3f;
	int trials = 100;
	int use_crossover = 1;
	unsigned int random_seed = (unsigned int)time(NULL);
	int sensors=2, actuators=CATEGORIES;
	int integers_only = 0;
	int no_of_test_examples;
	float test_performance;
	FILE *fp;
	int data_size=0, data_fields=0;
	char compile_command[256];
	int instruction_set[64], no_of_instructions=0;
	char * sensor_names[] = {
		"Age",
		"Sex",
		"Height",
		"Weight",
		"QRS Duration",
		"P-R Interval",
		"QRS",
		"T",
		"P",
		"QRST",
		"J",
		"Heart Rate",
		"Q Wave",
		"R Wave",
		"S Wave",
		"R' Wave",
		"S' Wave",
		"Intrinsic Deflections",
		"Ragged R Wave",
		"Diphasic R Wave",
		"Ragged P Wave",
		"Diphasic P Wave",
		"Ragged T Wave",
		"Diphasic T Wave",
		"DII Q Wave",
		"DII R Wave",
		"DII S Wave",
		"DII R' Wave",
		"DII S' Wave",
		"DII Intrinsic Deflections",
		"DII Ragged R Wave",
		"DII Diphasic R Wave",
		"DII Ragged P Wave",
		"DII Diphasic P Wave",
		"DII Ragged T Wave",
		"DII Diphasic T Wave",
		"DIII Q Wave",
		"DIII R Wave",
		"DIII S Wave",
		"DIII R' Wave",
		"DIII S' Wave",
		"DIII Intrinsic Deflections",
		"DIII Ragged R Wave",
		"DIII Diphasic R Wave",
		"DIII Ragged P Wave",
		"DIII Diphasic P Wave",
		"DIII Ragged T Wave",
		"DIII Diphasic T Wave",
		"AVR Q Wave",
		"AVR R Wave",
		"AVR S Wave",
		"AVR R' Wave",
		"AVR S' Wave",
		"AVR Intrinsic Deflections",
		"AVR Ragged R Wave",
		"AVR Diphasic R Wave",
		"AVR Ragged P Wave",
		"AVR Diphasic P Wave",
		"AVR Ragged T Wave",
		"AVR Diphasic T Wave",
		"AVL Q Wave",
		"AVL R Wave",
		"AVL S Wave",
		"AVL R' Wave",
		"AVL S' Wave",
		"AVL Intrinsic Deflections",
		"AVL Ragged R Wave",
		"AVL Diphasic R Wave",
		"AVL Ragged P Wave",
		"AVL Diphasic P Wave",
		"AVL Ragged T Wave",
		"AVL Diphasic T Wave",
		"AVF Q Wave",
		"AVF R Wave",
		"AVF S Wave",
		"AVF R' Wave",
		"AVF S' Wave",
		"AVF Intrinsic Deflections",
		"AVF Ragged R Wave",
		"AVF Diphasic R Wave",
		"AVF Ragged P Wave",
		"AVF Diphasic P Wave",
		"AVF Ragged T Wave",
		"AVF Diphasic T Wave",
		"V1 Q Wave",
		"V1 R Wave",
		"V1 S Wave",
		"V1 R' Wave",
		"V1 S' Wave",
		"V1 Intrinsic Deflections",
		"V1 Ragged R Wave",
		"V1 Diphasic R Wave",
		"V1 Ragged P Wave",
		"V1 Diphasic P Wave",
		"V1 Ragged T Wave",
		"V1 Diphasic T Wave",
		"V2 Q Wave",
		"V2 R Wave",
		"V2 S Wave",
		"V2 R' Wave",
		"V2 S' Wave",
		"V2 Intrinsic Deflections",
		"V2 Ragged R Wave",
		"V2 Diphasic R Wave",
		"V2 Ragged P Wave",
		"V2 Diphasic P Wave",
		"V2 Ragged T Wave",
		"V2 Diphasic T Wave",
		"V3 Q Wave",
		"V3 R Wave",
		"V3 S Wave",
		"V3 R' Wave",
		"V3 S' Wave",
		"V3 Intrinsic Deflections",
		"V3 Ragged R Wave",
		"V3 Diphasic R Wave",
		"V3 Ragged P Wave",
		"V3 Diphasic P Wave",
		"V3 Ragged T Wave",
		"V3 Diphasic T Wave",
		"V4 Q Wave",
		"V4 R Wave",
		"V4 S Wave",
		"V4 R' Wave",
		"V4 S' Wave",
		"V4 Intrinsic Deflections",
		"V4 Ragged R Wave",
		"V4 Diphasic R Wave",
		"V4 Ragged P Wave",
		"V4 Diphasic P Wave",
		"V4 Ragged T Wave",
		"V4 Diphasic T Wave",
		"V5 Q Wave",
		"V5 R Wave",
		"V5 S Wave",
		"V5 R' Wave",
		"V5 S' Wave",
		"V5 Intrinsic Deflections",
		"V5 Ragged R Wave",
		"V5 Diphasic R Wave",
		"V5 Ragged P Wave",
		"V5 Diphasic P Wave",
		"V5 Ragged T Wave",
		"V5 Diphasic T Wave",
		"V6 Q Wave",
		"V6 R Wave",
		"V6 S Wave",
		"V6 R' Wave",
		"V6 S' Wave",
		"V6 Intrinsic Deflections",
		"V6 Ragged R Wave",
		"V6 Diphasic R Wave",
		"V6 Ragged P Wave",
		"V6 Diphasic P Wave",
		"V6 Ragged T Wave",
		"V6 Diphasic T Wave",
		"DI JJ Wave",
		"DI Q Wave",
		"DI R Wave",
		"DI S Wave",
		"DI R' Wave",
		"DI S' Wave",
		"DI P Wave",
		"DI T Wave",
		"DI QRSA",
		"DI QRSTA",
		"DII JJ Wave",
		"DII Q Wave",
		"DII R Wave",
		"DII S Wave",
		"DII R' Wave",
		"DII S' Wave",
		"DII P Wave",
		"DII T Wave",
		"DII QRSA",
		"DII QRSTA",
		"DIII JJ Wave",
		"DIII Q Wave",
		"DIII R Wave",
		"DIII S Wave",
		"DIII R' Wave",
		"DIII S' Wave",
		"DIII P Wave",
		"DIII T Wave",
		"DIII QRSA",
		"DIII QRSTA",
		"AVR JJ Wave",
		"AVR Q Wave",
		"AVR R Wave",
		"AVR S Wave",
		"AVR R' Wave",
		"AVR S' Wave",
		"AVR P Wave",
		"AVR T Wave",
		"AVR QRSA",
		"AVR QRSTA",
		"AVL JJ Wave",
		"AVL Q Wave",
		"AVL R Wave",
		"AVL S Wave",
		"AVL R' Wave",
		"AVL S' Wave",
		"AVL P Wave",
		"AVL T Wave",
		"AVL QRSA",
		"AVL QRSTA",
		"AVF JJ Wave",
		"AVF Q Wave",
		"AVF R Wave",
		"AVF S Wave",
		"AVF R' Wave",
		"AVF S' Wave",
		"AVF P Wave",
		"AVF T Wave",
		"AVF QRSA",
		"AVF QRSTA",
		"V1 JJ Wave",
		"V1 Q Wave",
		"V1 R Wave",
		"V1 S Wave",
		"V1 R' Wave",
		"V1 S' Wave",
		"V1 P Wave",
		"V1 T Wave",
		"V1 QRSA",
		"V1 QRSTA",
		"V2 JJ Wave",
		"V2 Q Wave",
		"V2 R Wave",
		"V2 S Wave",
		"V2 R' Wave",
		"V2 S' Wave",
		"V2 P Wave",
		"V2 T Wave",
		"V2 QRSA",
		"V2 QRSTA",
		"V3 JJ Wave",
		"V3 Q Wave",
		"V3 R Wave",
		"V3 S Wave",
		"V3 R' Wave",
		"V3 S' Wave",
		"V3 P Wave",
		"V3 T Wave",
		"V3 QRSA",
		"V3 QRSTA",
		"V4 JJ Wave",
		"V4 Q Wave",
		"V4 R Wave",
		"V4 S Wave",
		"V4 R' Wave",
		"V4 S' Wave",
		"V4 P Wave",
		"V4 T Wave",
		"V4 QRSA",
		"V4 QRSTA",
		"V5 JJ Wave",
		"V5 Q Wave",
		"V5 R Wave",
		"V5 S Wave",
		"V5 R' Wave",
		"V5 S' Wave",
		"V5 P Wave",
		"V5 T Wave",
		"V5 QRSA",
		"V5 QRSTA",
		"V6 JJ Wave",
		"V6 Q Wave",
		"V6 R Wave",
		"V6 S Wave",
		"V6 R' Wave",
		"V6 S' Wave",
		"V6 P Wave",
		"V6 T Wave",
		"V6 QRSA",
		"V6 QRSTA",
		"Sensor 0",
		"Sensor 1",
		"Sensor 2"
	};
	
	char * actuator_names[] = {
		"Normal",
		"Ischemic changes",
		"Old Anterior Myocardial Infarction",
		"Old Interior Myocardial Infarction",
		"Sinus tachycardy",
		"Sinus bradycardy",
		"Ventricular Premature Contraction",
		"Supraventricular Premature Contraction",
		"Left bundle branch block",
		"Right bundle branch block",
		"1. degree AtrioVentricular block",
		"2. degree AV block",
		"3. degree AV block",
		"Left ventricule hypertrophy",
		"Atrial Fibrillation or Flutter",
		"Other"
	};	

	/* load the data */
	no_of_examples = load_data("arrhythmia.data",
							   arrhythmia_data, MAX_EXAMPLES,
							   &fields_per_example);

	/* create a test data set */
	no_of_test_examples = create_test_data(arrhythmia_data,
										   &no_of_examples,
										   fields_per_example,
										   test_data);

	sensors = fields_per_example-1;
	trials = no_of_examples;

	printf("Number of training examples: %d\n",no_of_examples);
	printf("Number of test examples: %d\n",no_of_test_examples);
	printf("Number of fields: %d\n",fields_per_example);

	/* create an instruction set */
	no_of_instructions =
		gprc_dynamic_instruction_set((int*)instruction_set);

	/* create a population */
	gprc_init_system(&sys, islands,
					 population_per_island,
					 rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 modules,
					 chromosomes,
					 min_value, max_value,
					 integers_only,
					 data_size, data_fields,
					 &random_seed,
					 instruction_set, no_of_instructions);

	gpr_xmlrpc_server("server.rb","arrhythmia",3573,
					  "./agent",
					  sensors, actuators);

	test_performance = 0;
	while (test_performance < 99) {
		/* use the training data */
		current_data_set = arrhythmia_data;

		/* evaluate each individual */
		gprc_evaluate_system(&sys,
							 trials,0,
							 (*evaluate_features));
		/* produce the next generation */
		gprc_generation_system(&sys,
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

		/* show the best fitness value */
		printf("Generation %05d  Fitness %.2f/%.2f%% ",
			   gen, gprc_best_fitness(&sys.island[0]),test_performance);
		for (i = 0; i < islands; i++) {
			printf("  %.3f",gprc_average_fitness(&sys.island[i]));
		}
		printf("\n");

		if (((gen % 100 == 0) && (gen>0)) || (test_performance > 99)) {
			gprc_draw_population("population.png",
								 640, 640, &sys.island[0]);

			gprc_plot_history_system(&sys,
									 GPR_HISTORY_FITNESS,
									 "fitness.png",
									 "Arrhythmia Classification Performance",
									 640, 480);

			gprc_plot_history_system(&sys,
									 GPR_HISTORY_AVERAGE,
									 "fitness_average.png",
									 "Arrhythmia Classification Average Performance",
									 640, 480);

			gprc_plot_history_system(&sys,
									 GPR_HISTORY_DIVERSITY,
									 "diversity.png", 
									 "Arrhythmia Classification Diversity",
									 640, 480);

			gprc_plot_fitness(&sys.island[0],
							  "fitness_histogram.png",
							  "Fitness Histogram",
							  640, 480);

			fp = fopen("agent.c","w");
			if (fp) {
				/* save the best program */
				gprc_c_program(&sys,
							   gprc_best_individual_system(&sys),
							   RUN_STEPS, 0, fp);
				fclose(fp);

				/* compile the program */
				sprintf(compile_command,
						"gcc -Wall -std=c99 -pedantic -o agent agent.c -lm");
				assert(system(compile_command)==0);
			}

			fp = fopen("fittest.dot","w");
			if (fp) {
				gprc_dot(gprc_best_individual_system(&sys),
						 &sys.island[0],
						 sensor_names,  actuator_names,
						 fp);
				fclose(fp);
			}
		}

		if (test_performance > 99) break;
		gen++;
	}

	/* free memory */
	gprc_free_system(&sys);
}

int main(int argc, char* argv[])
{	
	arrhythmia_classification();
	return 1;
}


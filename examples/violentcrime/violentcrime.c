/*
 Predicting levels of violent crime in America
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
#include "libgpr/som.h"

#define MAX_EXAMPLES 2000
#define MAX_TEST_EXAMPLES 200
#define MAX_FIELDS   130

/* the initial number of fields which can be ignored */
#define INITIAL_FIELDS 5

#define RUN_STEPS 2

/* dimension of the topological map */
#define SOM_DIMENSION 32

/* the number of learning itterations for the SOM */
#define SOM_LEARNING_ITTERATIONS 20

/* the number of crime levels */
#define CATEGORIES 10

float crime_data[MAX_EXAMPLES*MAX_FIELDS];
float test_data[MAX_TEST_EXAMPLES*MAX_FIELDS];
float * current_data_set;
float * current_som_outputs;

int no_of_examples = 0;
int fields_per_example = 0;

float max_crime_rate = 0;

gpr_som som;

/* create a test data set from the original data.
   The test data can be used to calculate a final fitness
   value, because it was not seen during training and so
   provides an indication of how well the system has generalised */
static int create_test_data(float * training_data, int * no_of_training_examples,
							int fields_per_example,
							float * test_data)
{
	int i,j,index;
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
			training_data[(j-1)*fields_per_example + j] = 
				training_data[j*fields_per_example + j];
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
							if ((valuestr[0]>='0') && (valuestr[0]<='9')) {
								value = atof(valuestr);
							}
						}
						else {
							value = GPR_MISSING_VALUE;
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
							   gprc_population * population,
							   int individual_index,
							   int mode)
{
	int i,j,itt,index,winner;
	float v,crime_rate,fitness=0,max,v1,v2;
	float dropout_rate=0.0f;
	gprc_function * f = &population->individual[individual_index];

	if (mode!=0) dropout_rate=0;

	for (i=0;i<trials;i++) {
		/* clear the state */
		gprc_clear_state(f,
						 population->rows, population->columns,
						 population->sensors, population->actuators);

		gprc_set_sensor(f,0,current_som_outputs[i*2]);
		gprc_set_sensor(f,1,current_som_outputs[i*2 + 1]);

		for (itt = 0; itt < RUN_STEPS; itt++) {
			/* run the program */
			gprc_run(f, population, dropout_rate, 0, 0);
		}
		/* how close is the output to the actual crime level? */
		crime_rate = current_data_set[i*fields_per_example+fields_per_example-1];
		index = (int)(crime_rate * CATEGORIES / max_crime_rate);
		if (index >= CATEGORIES) index = CATEGORIES-1;

		max = -1;
		winner=-1;
		v1 = 0;
		v2 = 0;
		for (j = 0; j < CATEGORIES; j++) {
			v = fabs(gprc_get_actuator(f,j,population->rows,
									   population->columns,
									   population->sensors));

			if (v > max) {
				v1 = v;
				max=v;
				winner=j;				
			}
			v2 += v;
		}

		if (mode==0) {
			v2 /= CATEGORIES;
			fitness += v1 - v2;
		}
		else {
			if (winner==index) {
				fitness += 100;
			}
		}
	}
	fitness /= trials;
	if (fitness < 0) fitness=0;
	return fitness;
}

static void violent_crimes_prediction()
{
	int islands = 4;
	int migration_interval = 200;
	int population_per_island = 64;
	int rows = 9, columns = 10;
	int i, gen=0;
	int connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int modules = 2;
	int chromosomes = 3;
	gprc_system sys;
	float min_value = -100;
	float max_value = 100;
	float elitism = 0.2f;
	float mutation_prob = 0.3f;
	int trials = 100;
	int use_crossover = 1;
	unsigned int random_seed = (unsigned int)time(NULL);
	unsigned int som_random_seed = 123;
	int sensors=2, som_sensors, actuators=CATEGORIES;
	int integers_only = 0;
	int no_of_test_examples;
	float test_performance;
	FILE * fp;
	char compile_command[256];
	int instruction_set[64], no_of_instructions=0;
	int * training_data_field_index;
	int inhibit_radius = 4;
	int excite_radius = 2;
	float som_learning_rate = 0.2f;
	float * som_training_outputs;
	float * som_test_outputs;
	/*
	char * sensor_names[] = {
		"Population",
		"household size",
		"Race percent black",
		"Race percent White",
		"Race percent Asian",
		"Race percent Hisp",
		"Age percent 12-21",
		"Age percent 12-29",
		"Age percent 16-24",
		"Age percent 65+",
		"Numb Urban",
		"Percent Urban",
		"Median Income",
		"Percent Wage",
		"Percent Farm Self",
		"Percent Inv Inc",
		"Percent Soc Sec",
		"Percent Pub Asst",
		"Percent Retire",
		"Median Family Income",
		"Per Capita Income",
		"White Per Cap",
		"Black Per Cap",
		"Indian Per Cap",
		"Asian Per Cap",
		"Other Per Cap",
		"Hisp Per Cap",
		"Num Under Pov",
		"Percent Pop Under Pov",
		"Percent Less 9th Grade",
		"Percent Not HS Grad",
		"Percent BS or More",
		"Percent Unemployed",
		"Percent Employ",
		"Percent Empl Manu",
		"Percent Empl Prof Serv",
		"Percent Occup Manu",
		"Percent Occup Mgmt Prof",
		"Male Percent Divorce",
		"Male Percent Nev Marr",
		"Female Pct Div",
		"Total Percent Div",
		"Percent Per Fam",
		"Percent Fam 2 Par",
		"Pct Kids 2 Par",
		"Percent Young Kids 2 Par",
		"Pct Teen 2 Par",
		"Percent Work Mom Young Kids",
		"Percent Work Mom",
		"Num Illeg",
		"Percent Illeg",
		"Num Immig",
		"Percent Immig Recent",
		"Percent Immig Rec5",
		"Percent Immig Rec8",
		"Percent Immig Rec10",
		"Percent Recent Immig",
		"Percent Rec Immig5",
		"Percent Rec Immig8",
		"Percent Rec Immig10",
		"Percent Speak Engl Only",
		"Percent Not Speak Engl Well",
		"Percent Larg House Fam",
		"Percent Larg House Occup",
		"Percent Per Occup Hous",
		"Percent Per Own Occ Hous",
		"Percent Per Rent Occ Hous",
		"Percent Pers Own Occup",
		"Percent Pers Dense Hous",
		"Percent Hous Less 3 BR",
		"Median Num BR",
		"Hous Vacant",
		"Percent Hous Occup",
		"Percent Hous Own Occ",
		"Percent Vacant Boarded",
		"Percent Vac More 6 Mos",
		"Med Yr Hous Built",
		"Percent Hous No Phone",
		"Percent WOFullPlumb",
		"Own Occ Low Quart",
		"Own Occ MedVal",
		"Own Occ Hi Quart",
		"Rent Low Q",
		"Rent Median",
		"Rent High Q",
		"Med Rent",
		"Med Rent Pct Hous Inc",
		"Med Own Cost Pct Inc",
		"Med Own Cost Pct Inc No Mtg",
		"Num In Shelters",
		"Num Street",
		"Percent Foreign Born",
		"Percent Born Same State",
		"Percent Same House 85",
		"Percent Same City 85",
		"Percent Same State 85",
		"Lemas Sworn FT",
		"Lemas SwFT Per Pop",
		"Lemas SwFT Field Ops",
		"Lemas SwFT Field Per Pop",
		"Lemas Total Req",
		"Lemas Tot Req Per Pop",
		"Polic Req Per Offic",
		"Polic Per Pop",
		"Racial Match Comm Pol",
		"Percent Polic White",
		"Percent Polic Black",
		"Percent Polic Hisp",
		"Percent Polic Asian",
		"Percent Polic Minor",
		"Offic Assgn Drug Units",
		"Num Kinds Drugs Seiz",
		"Polic Ave OT Worked",
		"Land Area",
		"Pop Dens",
		"Percent Use Pub Trans",
		"Polic Cars",
		"Polic Oper Budg",
		"Lemas Pct Polic On Patr",
		"Lemas Gang Unit Deploy",
		"Lemas Pct Offic Drug Un",
		"Polic Budg Per Pop"
	};
	*/
	char * som_sensor_names[] = { "x","y" };
	char * actuator_names[] = {
		"Level 1","Level 2","Level 3","Level 4","Level 5","Level 6","Level 7","Level 8","Level 9","Level 10",
		"Level 11","Level 12","Level 13","Level 14","Level 15","Level 16","Level 17","Level 18","Level 19","Level 20",
		"Level 21","Level 22","Level 23","Level 24","Level 25","Level 26","Level 27","Level 28","Level 29","Level 30",
		"Level 31","Level 32","Level 33","Level 34","Level 35","Level 36","Level 37","Level 38","Level 39","Level 40",
		"Level 41","Level 42","Level 43","Level 44","Level 45","Level 46","Level 47","Level 48","Level 49","Level 50",
	};

	/* load the data */
	no_of_examples = load_data("communities.data", crime_data,
							   MAX_EXAMPLES,
							   &fields_per_example);

	/* create a test data set */
	no_of_test_examples = create_test_data(crime_data,
										   &no_of_examples,
										   fields_per_example,
										   test_data);

	som_sensors = fields_per_example-1-INITIAL_FIELDS;
	trials = no_of_examples;

	printf("Number of training examples: %d\n",no_of_examples);
	printf("Number of test examples: %d\n",no_of_test_examples);
	printf("Number of fields: %d\n",fields_per_example);

	/* create a SOM */
	gpr_som_init(SOM_DIMENSION, som_sensors, &som);
	
	training_data_field_index = (int*)malloc(som_sensors*sizeof(int));

	/* initialise the SOM */
	for (i = 0; i < som_sensors; i++) {
		assert(gpr_som_init_sensor_from_data(&som,
											 i, INITIAL_FIELDS + i,
											 crime_data,
											 fields_per_example,
											 no_of_examples,
											 &random_seed) != -1);
		training_data_field_index[i] = INITIAL_FIELDS + i;
	}

	printf("Training SOM...");
	fflush(stdout);
	gpr_som_learn_from_data(&som,
							training_data_field_index,
							crime_data,
							fields_per_example,
							no_of_examples,
							SOM_LEARNING_ITTERATIONS,
							inhibit_radius, excite_radius,
							som_learning_rate,
							&som_random_seed,1);
	printf("Done\n");
	printf("Updating training outputs...");
	som_training_outputs = (float*)malloc(no_of_examples*2*sizeof(float));
	som_test_outputs = (float*)malloc(no_of_test_examples*2*sizeof(float));
	gpr_som_outputs_from_data(&som,
							  training_data_field_index,
							  crime_data,
							  fields_per_example,
							  no_of_examples,
							  som_training_outputs);
	printf("Done\n");
	printf("Updating test outputs...");
	gpr_som_outputs_from_data(&som,
							  training_data_field_index,
							  test_data,
							  fields_per_example,
							  no_of_test_examples,
							  som_test_outputs);
	printf("Done\n");
	free(training_data_field_index);

	/* find the maximum crime rate */
	max_crime_rate = 0;
	for (i = 0; i < no_of_examples; i++) {
		if (crime_data[(i*fields_per_example) +
					   fields_per_example - 1]) {
			if (crime_data[(i*fields_per_example) +
						   fields_per_example - 1] > max_crime_rate) {
				max_crime_rate =
					crime_data[(i*fields_per_example) +
							   fields_per_example - 1];
			}
		}
	}

	/*
	for (i=0;i<no_of_examples;i++) {
		printf("%.2f  %.2f\n",
			   som_training_outputs[i*2],
			   som_training_outputs[i*2 + 1]);
	}
	*/

	/* create an instruction set */
	no_of_instructions = gprc_equation_instruction_set((int*)instruction_set);

	/* create a population */
	gprc_init_system(&sys, islands,
					 population_per_island,
					 rows, columns,
					 sensors, actuators,
					 connections_per_gene,
					 modules,
					 chromosomes,
					 min_value, max_value,
					 integers_only, &random_seed,
					 instruction_set, no_of_instructions);

	gpr_xmlrpc_server("server.rb","crime",3573,
					  "./agent",
					  sensors, actuators);
	
	test_performance = 0;
	while (test_performance < 99) {
		/* use the training data */
		current_data_set = crime_data;
		current_som_outputs = som_training_outputs;

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
		current_som_outputs = som_training_outputs;
		test_performance = evaluate_features(no_of_test_examples,
											 &sys.island[0],
											 0,1);

		/* show the best fitness value */
		printf("Generation %05d  Fitness %.2f/%.2f%% ",gen, gprc_best_fitness(&sys.island[0]),test_performance);
		for (i = 0; i < islands; i++) {
			printf("  %.5f",gprc_average_fitness(&sys.island[i]));
		}
		printf("\n");

		if (((gen % 50 == 0) && (gen>0)) || (test_performance > 99)) {
			gprc_plot_history_system(&sys,
									 GPR_HISTORY_FITNESS,
									 "fitness.png", "Violent Crimes Prediction Performance",
									 640, 480);

			gprc_plot_history_system(&sys,
									 GPR_HISTORY_AVERAGE,
									 "fitness_average.png", "Violent Crimes Prediction Average Performance",
									 640, 480);

			gprc_plot_history_system(&sys,
									 GPR_HISTORY_DIVERSITY,
									 "diversity.png", "Violent Crimes Prediction Diversity",
									 640, 480);

			gprc_plot_fitness(&sys.island[0],
							  "fitness_histogram.png", "Violent Crimes Prediction Fitness Histogram",
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
						 som_sensor_names,  actuator_names,
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
	violent_crimes_prediction();
	return 1;
}


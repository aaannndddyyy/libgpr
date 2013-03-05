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

#ifndef GPRCM_H
#define GPRCM_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "globals.h"
#include "gpr.h"
#include "gprc.h"

#define GPRCM_MORPHOLOGY_ROWS                  6
#define GPRCM_MORPHOLOGY_COLUMNS               8
#define GPRCM_MORPHOLOGY_SENSORS               3
#define GPRCM_MORPHOLOGY_ACTUATORS             8
#define GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE  8

struct gprcm_func {
	/* instruction set for the morphology generator */
	int morphology_no_of_instructions;
	int morphology_instruction_set[64];

	/* morphology generator */
	gprc_function morphology;

	/* the main program */
	gprc_function program;
};
typedef struct gprcm_func gprcm_function;

/* represents a population */
struct gprcm_pop {
	/* the number of individuals in the population */
	int size;
	/* the number of rows and columns in the grid for each individual */
	int rows, columns;
	/* the number of sensors and actuators */
	int sensors, actuators;
	/* the maximum number of connections for each gene in the grid */
	int connections_per_gene;
	/* the number of ADF modules */
	int ADF_modules;
	/* the minimum and maximum constant values */
	float min_value, max_value;
	/* the number of chromosomes for each individual */
	int chromosomes;
	/* whether to only use integer maths */
	int integers_only;
	/* array containing individual programs */
	struct gprcm_func * individual;
	float * fitness;
	/* the fitness history for the population */
	struct gpr_hist history;
};
typedef struct gprcm_pop gprcm_population;

/* system containing a number of populations */
struct gprcm_sys {
	/* the number of sub-populations or islands */
	int size;
	/* the number of time steps after which
	   migrations between islands will occur */
	int migration_tick;
	/* population for each island */
	gprcm_population * island;
	/* the best fitness for each island */
	float * fitness;
	/* the fitness history for the system */
	struct gpr_hist history;
};
typedef struct gprcm_sys gprcm_system;

struct gprcm_env {
	/* the number of individuals in the population */
	int population_size;
	/* the maximum population size */
	int max_population_size;
	/* the number of rows and columns in the grid for each individual */
	int rows, columns;
	/* the number of sensors and actuators */
	int sensors, actuators;
	/* the maximum number of connections for each gene in the grid */
	int connections_per_gene;
	/* the number of ADF modules */
	int ADF_modules;
	/* the minimum and maximum constant values */
	float min_value, max_value;
	/* the number of chromosomes for each individual */
	int chromosomes;
	/* whether to only use integer maths */
	int integers_only;
	/* array containing individual programs */
	struct gprcm_func * individual;
	/* the number of matings */
	int matings;
	/* index numbers of mating parents */
	int * mating;
};
typedef struct gprcm_env gprcm_environment;

int gprcm_morphology_instruction_set(int * instruction_set);
void gprcm_init(gprcm_function * f,
				int rows, int columns, int sensors, int actuators,
				int connections_per_gene, int ADF_modules,
				unsigned int * random_seed);
void gprcm_init_sensor_sources(gprcm_system * system,
							   int no_of_sensor_sources,
							   unsigned int * random_seed);
void gprcm_init_actuator_destinations(gprcm_system * system,
									  int no_of_actuator_destinations,
									  unsigned int * random_seed);
void gprcm_clear_state(gprcm_function * f,
					   int rows, int columns,
					   int sensors, int actuators);
void gprcm_free(gprcm_function * f);
void gprcm_random(gprcm_function * f,
				  int rows, int columns,
				  int sensors, int actuators,
				  int connections_per_gene,
				  float min_value, float max_value,
				  int integers_only, unsigned int * random_seed,
				  int * instruction_set, int no_of_instructions);
int gprcm_validate(gprcm_function * f,
				   int rows, int columns,
				   int sensors, int actuators,
				   int connections_per_gene,
				   int integers_only,
				   int * instruction_set,
				   int no_of_instructions);
void gprcm_run(gprcm_function * f, gprcm_population * population,
			   float dropout_prob, int dynamic,
			   float (*custom_function)(float,float,float));
void gprcm_run_environment(gprcm_function * f,
						   gprcm_environment * population,
						   float dropout_prob, int dynamic,
						   float (*custom_function)(float,float,float));
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
						   int no_of_instructions);
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
							int no_of_instructions);
void gprcm_free_population(gprcm_population * population);
void gprcm_free_environment(gprcm_environment * population);
void gprcm_copy(gprcm_function * source, gprcm_function * dest,
				int rows, int columns, int connections_per_gene,
				int sensors, int actuators);
void gprcm_evaluate(gprcm_population * population,
					int time_steps, int reevaluate,
					float (*evaluate_program)
					(int,gprcm_population*,int,int));
float gprcm_best_fitness(gprcm_population * population);
float gprcm_worst_fitness(gprcm_population * population);
float gprcm_average_fitness(gprcm_population * population);
gprcm_function * gprcm_best_individual(gprcm_population * population);
void gprcm_set_sensor(gprcm_function * f, int index, float value);
float gprcm_get_sensor(gprcm_function * f, int index);
int gprcm_get_sensor_source(gprcm_function * f, int index);
float gprcm_get_actuator(gprcm_function * f, int index,
						 int rows, int columns, int sensors);
int gprcm_get_actuator_destination(gprcm_function * f, int index);
void gprcm_sort(gprcm_population * population);
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
				gprcm_function *child);
void gprcm_generation(gprcm_population * population,
					  float elitism,
					  float mutation_prob,
					  int use_crossover, unsigned int * random_seed,
					  int * instruction_set, int no_of_instructions);
int gprcm_save(gprcm_function * f,
			   int rows, int columns,
			   int connections_per_gene,
			   int sensors, int actuators,
			   FILE * fp);
int gprcm_load(gprcm_function * f,
			   int rows, int columns,
			   int connections_per_gene,
			   int sensors, int actuators,
			   FILE * fp);
void gprcm_save_population(gprcm_population * population,
						   FILE * fp);
void gprcm_save_environment(gprcm_environment * population,
							FILE * fp);
void gprcm_load_population(gprcm_population * population,
						   FILE * fp,
						   int * instruction_set,
						   int no_of_instructions);
void gprcm_load_environment(gprcm_environment * population,
							FILE * fp,
							int * instruction_set,
							int no_of_instructions);
void print_gprcm(gprcm_function * f,
				 int ADF_module,
				 int rows, int columns,
				 int sensors, int actuators,
				 int connections_per_gene,
				 int integers_only);
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
				   FILE * fp);
void gprcm_c_program(gprcm_system * system,
					 gprcm_function * f,
					 int itterations,
					 int dynamic,
					 FILE * fp);
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
					   int * instruction_set, int no_of_instructions);
void gprcm_free_system(gprcm_system * system);
void gprcm_evaluate_system(gprcm_system * system,
						   int time_steps, int reevaluate,
						   float (*evaluate_program)
						   (int,gprcm_population*,int,int));
void gprcm_generation_system(gprcm_system * system,
							 int migration_interval,
							 float elitism,
							 float mutation_prob,
							 int use_crossover,
							 unsigned int * random_seed,
							 int * instruction_set,
							 int no_of_instructions);
void gprcm_sort_system(gprcm_system * system);
float gprcm_best_fitness_system(gprcm_system * system);
gprcm_function * gprcm_best_individual_system(gprcm_system * system);
void gprcm_load_system(gprcm_system * system,
					   FILE * fp,
					   int * instruction_set, int no_of_instructions);
void gprcm_save_system(gprcm_system *system, FILE * fp);
int gprcm_default_instruction_set(int * instruction_set);
int gprcm_equation_instruction_set(int * instruction_set);
int gprcm_equation_dynamic_instruction_set(int * instruction_set);
int gprcm_advanced_instruction_set(int * instruction_set);
int gprcm_dynamic_instruction_set(int * instruction_set);
int gprcm_associative_instruction_set(int * instruction_set);
int gprcm_plot_history(gprcm_population * population,
					   int history_type,
					   char * filename, char * title,
					   int image_width, int image_height);
int gprcm_plot_history_system(gprcm_system * sys,
							  int history_type,
							  char * filename, char * title,
							  int image_width, int image_height);
float gprcm_median_fitness(gprcm_population * population);
int gprcm_plot_fitness(gprcm_population * population,
					   char * filename, char * title,
					   int image_width, int image_height);
int gprcm_plot_fitness_system(gprcm_system * sys,
							  char * filename, char * title,
							  int image_width, int image_height);
void gprcm_dot(gprcm_function * f, gprcm_population * population,
			   char * sensor_names[],
			   char * actuator_names[],
			   FILE * fp);
int gprcm_get_sensors(int ADF_module, int sensors);
int gprcm_get_actuators(int ADF_module, int actuators);
int gprcm_compress_ADF(gprcm_function * f,
					   int ADF_module,
					   int start_index,
					   int rows, int columns,
					   int connections_per_gene,
					   int sensors, int actuators,
					   float min_value, float max_value,
					   int max_depth,
					   int random_termination);
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
					   int random_termination);
int gprcm_mate_environment(gprcm_environment * population,
						   int parent1_index,
						   int parent2_index,
						   float mutation_prob,
						   int use_crossover,
						   int * instruction_set,
						   int no_of_instructions);
void gprcm_death(gprcm_environment * population,
				 int victim_index);
int gprcm_functions_are_equal(gprcm_function * f1,
							  gprcm_function * f2,
							  int rows, int columns,
							  int connections_per_gene,
							  int modules, int sensors);
void gprcm_show_population(unsigned char * img,
						   int img_width, int img_height, int bpp,
						   gprcm_population * population);
void gprcm_show_environment(unsigned char * img,
							int img_width, int img_height, int bpp,
							gprcm_environment * population);
void gprcm_draw_population(char * filename,
						   int img_width, int img_height,
						   gprcm_population * population);

#endif

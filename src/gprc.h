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

#ifndef GPRC_H
#define GPRC_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "globals.h"
#include "gpr.h"

enum {
	GPRC_GENE_FUNCTION_TYPE = 0,
	GPRC_GENE_CONSTANT,
	GPRC_GENE_IMAGINARY
};

/* the range of weight values */
#define GPRC_MIN_WEIGHT -1.0
#define GPRC_MAX_WEIGHT  1.0

/* the number of initial values before the list of
   connections and weights */
#define GPRC_INITIAL     3

/* the number of values per connection,
 with the first value being the connection location itself */
#define GPRC_WEIGHTS_PER_CONNECTION 2

/* the size of each gene within the cartesian grid */
#define GPRC_GENE_SIZE(connections) ((GPRC_INITIAL) + \
									 (connections* \
									  GPRC_WEIGHTS_PER_CONNECTION))

/* this structure contains the cartesian grid */
struct gprc_mod {
	/* defines the grid functions, known as genes */
	float * gene;
	/* the state value for each grid location */
	float * state;
	/* whether each gene is currently being used
	   as part of the inputs -> outputs transform */
	unsigned char * used;
};
typedef struct gprc_mod gprc_ADF_module;

/* represents a function */
struct gprc_func {
	/* the number of ADF modules */
	int ADF_modules;
	/* stores the cartesian grids for the main
	   program and any ADF modules */
	gprc_ADF_module genome[GPRC_MAX_ADF_MODULES+1];

	/* sensors can be redirected to various other
	   sources (eg a larger input set).  This defines
	   what inputs are picked from a larger set. */
	int no_of_sensor_sources;
	int * sensor_source;

	/* actuators can be redirected to various other
	   outputs (eg a larger output set).  This defines
	   what outputs are picked from a larger set. */
	int no_of_actuator_destinations;
	int * actuator_destination;

	/* random number seed */
	unsigned int random_seed;

	/* the age of this individual in generations */
	int age;

	/* temporary array */
	int * temp_genes;
};
typedef struct gprc_func gprc_function;


/* represents a population */
struct gprc_pop {
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
	struct gprc_func * individual;
	float * fitness;
	/* the fitness history for the population */
	struct gpr_hist history;
};
typedef struct gprc_pop gprc_population;

/* system containing a number of populations */
struct gprc_sys {
	/* the number of sub-populations or islands */
	int size;
	/* the number of time steps after which
	   migrations between islands will occur */
	int migration_tick;
	/* population for each island */
	gprc_population * island;
	/* the best fitness for each island */
	float * fitness;
	/* the fitness history for the system */
	struct gpr_hist history;
};
typedef struct gprc_sys gprc_system;

struct gprc_env {
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
	struct gprc_func * individual;
	/* the number of matings */
	int matings;
	/* index numbers of mating parents */
	int * mating;
};
typedef struct gprc_env gprc_environment;

int get_ADF_args(gprc_function * f, int ADF_module);
void gprc_tidy(gprc_function * f,
			   int rows, int columns,
			   int sensors, int actuators,
			   int connections_per_gene,
			   float min_value, float max_value);
void gprc_unique_outputs(gprc_function * f,
						 int rows, int columns,
						 int connections_per_gene,
						 int sensors, int actuators,
						 unsigned int * random_seed);
void gprc_valid_logical_operators(gprc_function * f,
								  int rows, int columns,
								  int connections_per_gene,
								  int sensors,
								  unsigned int * random_seed);
void gprc_dot_label(gprc_function * f,
					int ADF_module,
					int rows, int columns,
					int connections_per_gene,
					int sensors, int actuators,
					char * sensor_names[],
					char * actuator_names[],
					FILE * fp);
void gprc_dot_links(gprc_function * f,
					int ADF_module,
					int rows, int columns,
					int connections_per_gene,
					int sensors, int actuators,
					FILE * fp);
int gprc_contains_ADFs(gprc_function * f,
					   int ADF_module, int call_ADF_module,
					   int rows, int columns,
					   int connections_per_gene,
					   int sensors);
void gprc_used_functions(gprc_function * f,
						 int rows, int columns,
						 int connections_per_gene,
						 int sensors, int actuators);
void gprc_valid_ADFs(gprc_function * f,
					 int rows, int columns,
					 int connections_per_gene,
					 int sensors,
					 float min_value, float max_value);
void gprc_remove_ADFs(gprc_function * f,
					  int rows, int columns,
					  int connections_per_gene);
void gprc_init(gprc_function * f,
			   int rows, int columns, int sensors, int actuators,
			   int connections_per_gene, int ADF_modules,
			   unsigned int * random_seed);
void gprc_init_sensor_sources(gprc_system * system,
							  int no_of_sensor_sources,
							  unsigned int * random_seed);
void gprc_init_actuator_destinations(gprc_system * system,
									 int no_of_actuator_destinations,
									 unsigned int * random_seed);
void gprc_clear_state(gprc_function * f,
					  int rows, int columns,
					  int sensors, int actuators);
void gprc_free(gprc_function * f);
void gprc_random(gprc_function * f,
				 int rows, int columns,
				 int sensors, int actuators,
				 int connections_per_gene,
				 float min_value, float max_value,
				 int integers_only, unsigned int * random_seed,
				 int * instruction_set, int no_of_instructions);
void gprc_mutate(gprc_function * f,
				 int rows, int columns,
				 int sensors, int actuators,
				 int connections_per_gene,
				 int chromosomes,
				 float prob,
				 float chromosomes_prob,
				 float min_value, float max_value,
				 int integers_only,
				 int * instruction_set, int no_of_instructions);
int gprc_validate(gprc_function * f,
				  int rows, int columns,
				  int sensors, int actuators,
				  int connections_per_gene,
				  int integers_only,
				  int * instruction_set,
				  int no_of_instructions);
void gprc_run_float(gprc_function * f,
					int ADF_module,
					int rows, int columns,
					int connections_per_gene,
					int sensors, int actuators,
					float dropout_prob,
					int dynamic,
					float (*custom_function)(float,float,float));
void gprc_run_int(gprc_function * f,
				  int ADF_module,
				  int rows, int columns,
				  int connections_per_gene,
				  int sensors, int actuators,
				  float dropout_prob,
				  int dynamic,
				  float (*custom_function)(float,float,float));
void gprc_run(gprc_function * f, gprc_population * population,
			  float dropout_prob, int dynamic,
			  float (*custom_function)(float,float,float));
void gprc_run_environment(gprc_function * f,
						  gprc_environment * population,
						  float dropout_prob, int dynamic,
						  float (*custom_function)(float,float,float));
void gprc_init_population(gprc_population * population,
						  int size,
						  int rows, int columns,
						  int sensors, int actuators,
						  int connections_per_gene,
						  int ADF_modules,
						  int chromosomes,
						  float min_value, float max_value,
						  int integers_only, unsigned int * random_seed,
						  int * instruction_set, int no_of_instructions);
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
						   int no_of_instructions);
void gprc_free_population(gprc_population * population);
void gprc_free_environment(gprc_environment * population);
void gprc_copy(gprc_function * source, gprc_function * dest,
			   int rows, int columns, int connections_per_gene,
			   int sensors, int actuators);
void gprc_crossover(gprc_function *parent1, gprc_function *parent2,
					int rows, int columns,
					int sensors, int actuators,
					int connections_per_gene,
					int chromosomes,
					int allocate_memory,
					gprc_function *child);
void gprc_evaluate(gprc_population * population,
				   int time_steps, int reevaluate,
				   float (*evaluate_program)
				   (int,gprc_population*,int,int));
float gprc_best_fitness(gprc_population * population);
float gprc_worst_fitness(gprc_population * population);
float gprc_average_fitness(gprc_population * population);
gprc_function * gprc_best_individual(gprc_population * population);
void gprc_set_sensor(gprc_function * f, int index, float value);
void gprc_set_sensor_complex(gprc_function * f, int index,
							 float real, float imaginary,
							 int no_of_sensors, int no_of_actuators,
							 int rows, int columns);
void gprc_set_sensor_colour(gprc_function * f, int index,
							unsigned char R, unsigned char G, unsigned char B,
							int no_of_sensors, int no_of_actuators,
							int rows, int columns);
float gprc_get_sensor(gprc_function * f, int index);
int gprc_get_sensor_source(gprc_function * f, int index);
float gprc_get_actuator(gprc_function * f, int index,
						int rows, int columns, int sensors);
void gprc_get_actuator_complex(gprc_function * f, int index,
							   int rows, int columns,
							   int sensors, int actuators,
							   float * real, float * imaginary);
void gprc_get_actuator_colour(gprc_function * f, int index,
							  int rows, int columns,
							  int sensors, int actuators,
							  unsigned char * R,
							  unsigned char * G,
							  unsigned char * B);
int gprc_get_actuator_destination(gprc_function * f, int index);
void gprc_sort(gprc_population * population);
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
			   gprc_function *child);
void gprc_generation(gprc_population * population,
					 float elitism,
					 float mutation_prob,
					 int use_crossover, unsigned int * random_seed,
					 int * instruction_set, int no_of_instructions);
int gprc_save(gprc_function * f,
			  int rows, int columns,
			  int connections_per_gene,
			  int sensors, int actuators,
			  FILE * fp);
int gprc_load(gprc_function * f,
			  int rows, int columns,
			  int connections_per_gene,
			  int sensors, int actuators,
			  FILE * fp);
void gprc_save_population(gprc_population * population,
						  FILE * fp);
void gprc_save_environment(gprc_environment * population,
						  FILE * fp);
void gprc_load_population(gprc_population * population,
						  FILE * fp,
						  int * instruction_set,
						  int no_of_instructions);
void gprc_load_environment(gprc_environment * population,
						   FILE * fp,
						   int * instruction_set,
						   int no_of_instructions);
void print_gprc(gprc_function * f,
				int ADF_module,
				int rows, int columns,
				int sensors, int actuators,
				int connections_per_gene,
				int integers_only);

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
					   FILE * fp);
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
				  FILE * fp);
void gprc_c_program_base(int rows, int columns,
						 int connections_per_gene,
						 int sensors, int actuators,
						 int ADF_modules, int integers_only,
						 gprc_function * f,
						 int itterations,
						 int dynamic,
						 FILE * fp);
void gprc_c_program(gprc_system * system,
					gprc_function * f,
					int itterations,
					int dynamic,
					FILE * fp);
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
					  int * instruction_set, int no_of_instructions);
void gprc_free_system(gprc_system * system);
void gprc_evaluate_system(gprc_system * system,
						  int time_steps, int reevaluate,
						  float (*evaluate_program)
						  (int,gprc_population*,int,int));
void gprc_generation_system(gprc_system * system,
							int migration_interval,
							float elitism,
							float mutation_prob,
							int use_crossover,
							unsigned int * random_seed,
							int * instruction_set,
							int no_of_instructions);
void gprc_sort_system(gprc_system * system);
float gprc_best_fitness_system(gprc_system * system);
gprc_function * gprc_best_individual_system(gprc_system * system);
void gprc_load_system(gprc_system * system,
					  FILE * fp,
					  int * instruction_set, int no_of_instructions);
void gprc_save_system(gprc_system *system, FILE * fp);
int gprc_default_instruction_set(int * instruction_set);
int gprc_equation_instruction_set(int * instruction_set);
int gprc_equation_dynamic_instruction_set(int * instruction_set);
int gprc_advanced_instruction_set(int * instruction_set);
int gprc_dynamic_instruction_set(int * instruction_set);
int gprc_associative_instruction_set(int * instruction_set);
int gprc_plot_history(gprc_population * population,
					  int history_type,
					  char * filename, char * title,
					  int image_width, int image_height);
int gprc_plot_history_system(gprc_system * sys,
							 int history_type,
							 char * filename, char * title,
							 int image_width, int image_height);
float gprc_median_fitness(gprc_population * population);
int gprc_plot_fitness(gprc_population * population,
					  char * filename, char * title,
					  int image_width, int image_height);
int gprc_plot_fitness_system(gprc_system * sys,
							 char * filename, char * title,
							 int image_width, int image_height);
void gprc_dot(gprc_function * f, gprc_population * population,
			  char * sensor_names[],
			  char * actuator_names[],
			  FILE * fp);
int gprc_get_sensors(int ADF_module, int sensors);
int gprc_get_actuators(int ADF_module, int actuators);
int gprc_compress_ADF(gprc_function * f,
					  int ADF_module,
					  int start_index,
					  int rows, int columns,
					  int connections_per_gene,
					  int sensors, int actuators,
					  float min_value, float max_value,
					  int max_depth,
					  int random_termination);
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
					  int random_termination);
int gprc_mate_environment(gprc_environment * population,
						  int parent1_index,
						  int parent2_index,
						  float mutation_prob,
						  int use_crossover,
						  int * instruction_set, int no_of_instructions);
void gprc_death(gprc_environment * population,
				int victim_index);
int gprc_functions_are_equal(gprc_function * f1,
							 gprc_function * f2,
							 int rows, int columns,
							 int connections_per_gene,
							 int modules,
							 int sensors);
int gprc_no_of_dynamic_functions(gprc_function * f,
								 int rows, int columns,
								 int sensors, int actuators,
								 int connections_per_gene);
void gprc_show_genome(unsigned char * img,
					  int img_width, int img_height, int bpp,
					  int tx, int ty, int bx, int by,
					  gprc_function * f,
					  int rows, int columns,
					  int sensors, int actuators,
					  int connections_per_gene,
					  int module);
void gprc_show_population(unsigned char * img,
						  int img_width, int img_height, int bpp,
						  gprc_population * population);
void gprc_draw_population(char * filename,
						  int img_width, int img_height,
						  gprc_population * population);

#endif

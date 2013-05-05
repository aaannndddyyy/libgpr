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

#ifndef GPR_H
#define GPR_H

#define PNG_DEBUG 3

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <omp.h>
#include "globals.h"
#include <zlib.h>
#include "pnglite.h"
#include "gpr_data.h"

/* types of function */
enum {
	/* terminals */
	GPR_FUNCTION_NONE = 0,
	GPR_FUNCTION_VALUE,

	/* simple instruction set */
	GPR_FUNCTION_ADD,
	GPR_FUNCTION_SUBTRACT,
	GPR_FUNCTION_NEGATE,
	GPR_FUNCTION_MULTIPLY,
	GPR_FUNCTION_WEIGHT,
	GPR_FUNCTION_DIVIDE,
	GPR_FUNCTION_MODULUS,
	GPR_FUNCTION_FLOOR,
	GPR_FUNCTION_AVERAGE,
	GPR_FUNCTION_NOOP1,
	GPR_FUNCTION_NOOP2,
	GPR_FUNCTION_NOOP3,
	GPR_FUNCTION_NOOP4,
	GPR_FUNCTION_GREATER_THAN,
	GPR_FUNCTION_LESS_THAN,
	GPR_FUNCTION_EQUALS,
	GPR_FUNCTION_AND,
	GPR_FUNCTION_OR,
	GPR_FUNCTION_XOR,
	GPR_FUNCTION_NOT,

	GPR_FUNCTION_DATA_PUSH,
	GPR_FUNCTION_DATA_POP,
	GPR_FUNCTION_DATA_GET,
	GPR_FUNCTION_DATA_SET,

	GPR_FUNCTION_SET,
	GPR_FUNCTION_GET,

	/* advanced instruction set with more functions */
	GPR_FUNCTION_EXP,
	GPR_FUNCTION_SQUARE_ROOT,
	GPR_FUNCTION_ABS,
	GPR_FUNCTION_SINE,
	GPR_FUNCTION_ARCSINE,
	GPR_FUNCTION_COSINE,
	GPR_FUNCTION_ARCCOSINE,
	GPR_FUNCTION_POW,
	GPR_FUNCTION_SIGMOID,
	GPR_FUNCTION_MIN,
	GPR_FUNCTION_MAX,

	/* dynamic programming functions */
	GPR_FUNCTION_COPY_FUNCTION,
	GPR_FUNCTION_COPY_CONSTANT,
	GPR_FUNCTION_COPY_STATE,
	GPR_FUNCTION_COPY_BLOCK,
	GPR_FUNCTION_COPY_CONNECTION1,
	GPR_FUNCTION_COPY_CONNECTION2,
	GPR_FUNCTION_COPY_CONNECTION3,
	GPR_FUNCTION_COPY_CONNECTION4,

	GPR_FUNCTION_HEBBIAN,
	GPR_FUNCTION_DEFUN,
	GPR_FUNCTION_MAIN,
	GPR_FUNCTION_ADF,
	GPR_FUNCTION_ARG,
	GPR_FUNCTION_PROGRAM,
	GPR_FUNCTION_CUSTOM
};

#define GPR_FUNCTION_TYPES_ADVANCED (GPR_FUNCTION_MAX+1)
#define GPR_FUNCTION_TYPES_DYNAMIC (GPR_FUNCTION_COPY_CONNECTION4+1)
#define GPR_FUNCTION_TYPES_SIMPLE (GPR_FUNCTION_GET+1)
#define GPR_FUNCTION_TYPES_CARTESIAN (GPR_FUNCTION_NOT+1)

/* the function at the top of the tree */
#define GPR_TOP_LEVEL_FUNCTION GPR_FUNCTION_PROGRAM

enum {
	GPR_ORACLE_REGISTER=0,
	GPR_ORACLE_SENSOR,
	GPR_ORACLE_ACTUATOR,
	GPR_ORACLES
};

/* complex number value */
/*
struct gpr_val {
	float value;
	float imaginary;
};
typedef struct gpr_val gpr_value;*/

/* represents a function */
struct gpr_func {
	/* the type of function */
	unsigned short function_type;
	/* if this is a terminal then this is the value */
	float value;
	/* the number of function arguments */
	int argc;
	/* sub-functions */
	struct gpr_func * argv[GPR_MAX_ARGUMENTS];
};
typedef struct gpr_func gpr_function;

/* machine state */
struct gpr_st {
	/* the number of registers available to the program */
	int no_of_registers;
	float * registers;

	/* the number of sensors */
	int no_of_sensors;
	float * sensors;

	int no_of_sensor_sources;
	int * sensor_source;

	/* the number of actuators */
	int no_of_actuators;
	float * actuators;

	/* actuators can be redirected to various other
	   outputs (eg a larger output set).  This defines
	   what outputs are picked from a larger set. */
	int no_of_actuator_destinations;
	int * actuator_destination;

	/* random number seed */
	unsigned int random_seed;

	/* the age of the program in generations */
	unsigned int age;

	/* pointers to automatically defined functions */
	struct gpr_func * ADF[GPR_MAX_ARGUMENTS];
	
	/* number of arguments for each ADF */
	int ADF_argc[GPR_MAX_ARGUMENTS];

	/* temporary arguments to be passed to an ADF */
	float temp_ADF_arg[GPR_MAX_CALL_DEPTH][GPR_MAX_ARGUMENTS];

	/* a data store */
	gpr_data data;
};
typedef struct gpr_st gpr_state;

struct gpr_hist {
	int index;
	int interval;
	int tick;
	/* fitness */
	float log[GPR_MAX_HISTORY];
	/* diversity */
	float diversity[GPR_MAX_HISTORY];
	/* average fitness */
	float average[GPR_MAX_HISTORY];
};
typedef struct gpr_hist gpr_history;

/* represents a population */
struct gpr_pop {
	/* the number of individuals in the population */
	int size;
	/* array containing individual programs */
	struct gpr_func * individual;
	/* the states for each program */
	struct gpr_st * state;
	/* the fitness of each program */
	float * fitness;
	/* data store parameters */
	int data_size, data_fields;
	/* the fitness history for the population */
	struct gpr_hist history;
};
typedef struct gpr_pop gpr_population;

struct gpr_env {
	/* the maximum population size */
	int max_population_size;
	/* the number of individuals in the population */
	int population_size;
	/* array containing individual programs */
	struct gpr_func * individual;
	/* the states for each program */
	struct gpr_st * state;
	/* the number of matings */
	int matings;
	/* index numbers of mating parents */
	int * mating;
	/* data store parameters */
	int data_size, data_fields;
};
typedef struct gpr_env gpr_environment;

/* system containing a number of populations */
struct gpr_sys {
	/* the number of sub-populations or islands */
	int size;
	/* the number of time steps after which
	   migrations between islands will occur */
	int migration_tick;
	/* population for each island */
	gpr_population * island;
	/* the best fitness for each island */
	float * fitness;
	/* the fitness history for the system */
	struct gpr_hist history;
};
typedef struct gpr_sys gpr_system;

/*void gpr_clear_value(gpr_value * v);*/
float gpr_mutate_value(float value,
					   float percent,
					   unsigned int * random_seed);
int is_nan(float v);
int rand_num(unsigned int * seed);
void gpr_validate(gpr_function * f, int depth,
				  int min_depth, int max_depth,
				  int ADFs,
				  int tree_index, int max_tree_index,
				  int * result);
void gpr_max_depth(gpr_function * f, int depth, int * max_depth);
void gpr_prune(gpr_function * f, int depth, int max_depth,
			   float min_value, float max_value,
			   unsigned int * random_seed);
float gpr_run(gpr_function * f, gpr_state * state,
			  float (*custom_function)(float,float,float));
void gpr_nodes(gpr_function * f, int * ctr);
void gpr_init(gpr_function * f);
void gpr_init_state(gpr_state * state,
					int registers,
					int sensors,
					int actuators,
					int data_size, int data_fields,
					unsigned int * random_seed);
void gpr_init_sensor_sources(gpr_system * system,
							 int no_of_sensors,
							 int no_of_sensor_sources,
							 unsigned int * random_seed);
void gpr_init_actuator_destinations(gpr_system * system,
									int no_of_actuators,
									int no_of_actuator_destinations,
									unsigned int * random_seed);
void gpr_init_population(gpr_population * population,
						 int size,
						 int registers,
						 int sensors,
						 int actuators,
						 int max_tree_depth,
						 float min_value, float max_value,
						 int integers_only,
						 int ADFs,
						 int data_size, int data_fields,
						 unsigned int * random_seed,
						 int * instruction_set, int no_of_instructions);
void gpr_init_environment(gpr_environment * population,
						  int max_population_size,
						  int initial_population_size,
						  int registers,
						  int sensors, int actuators,
						  int max_tree_depth,
						  float min_value, float max_value,
						  int integers_only,
						  int ADFs,
						  int data_size, int data_fields,
						  unsigned int * random_seed,
						  int * instruction_set,
						  int no_of_instructions);
void gpr_free(gpr_function * f);
void gpr_free_state(gpr_state * state);
void gpr_free_population(gpr_population * population);
void gpr_free_environment(gpr_environment * population);
void gpr_copy(gpr_function * source, gpr_function * dest);
int gpr_default_instruction_set(int * instruction_set);
int gpr_equation_instruction_set(int * instruction_set);
int gpr_simple_instruction_set(int * instruction_set);
void gpr_random(gpr_function * f,
				int depth, int min_depth, int max_depth,
				float branching_prob,
				float min_value, float max_value,
				int integers_only,
				unsigned int * random_seed,
				int * instruction_set, int no_of_instructions);
void gpr_mutate(gpr_function * f, int depth, int max_depth, float prob,
				float min_value, float max_value,
				int integers_only, int ADFs,
				int tree_index, int max_tree_index,
				unsigned int * random_seed,
				int * instruction_set, int no_of_instructions);
void gpr_mutate_state(gpr_state * state,
					  float mutation_prob,
					  unsigned int * random_seed);
int gpr_crossover(gpr_function * parent1, gpr_function * parent2,
				  gpr_function * child, int max_depth,
				  float min_value, float max_value,
				  unsigned int * random_seed);
int gpr_crossover_ADF(gpr_function * parent1, gpr_function * parent2,
					  gpr_function * child, int max_depth,
					  float min_value, float max_value,
					  unsigned int * random_seed);
void gpr_crossover_sources(gpr_state * parent1,gpr_state * parent2,
						   gpr_state * child,
						   unsigned int * random_seed);
void gpr_mate(gpr_function * parent1, gpr_function * parent2,
			  int max_depth,
			  float min_value, float max_value,
			  float mutation_prob,
			  float pure_mutant_prob,
			  int integers_only,
			  int ADFs,
			  unsigned int * random_seed,
			  int * instruction_set, int no_of_instructions,
			  gpr_function * child,
			  gpr_state * child_state);
void gpr_evaluate(gpr_population * population,
				  int time_steps, int reevaluate,
				  float (*evaluate_program)(int,gpr_function*,gpr_state*,int));
void gpr_sort(gpr_population * population);
void gpr_generation(gpr_population * population,
					float elitism,
					int max_tree_depth,
					float min_value, float max_value,
					float mutation_prob,
					float pure_mutant_prob,
					int integers_only,
					int ADFs,
					int * instruction_set, int no_of_instructions);
float gpr_best_fitness(gpr_population * population);
float gpr_worst_fitness(gpr_population * population);
float gpr_average_fitness(gpr_population * population);
gpr_function * gpr_best_individual(gpr_population * population);
void gpr_dot(gpr_function * f, FILE * fp);
void gpr_set_sensor(gpr_state * state, int index, float value);
float gpr_get_sensor(gpr_state * state, int index);
int gpr_get_sensor_source(gpr_state * state, int index);
float gpr_get_actuator(gpr_state * state, int index);
int gpr_get_actuator_destination(gpr_state * state, int index);
void gpr_save(gpr_function * f, FILE * fp);
void gpr_save_population(gpr_population *population, FILE * fp);
int gpr_load(gpr_function * f, FILE * fp);
int gpr_load_population(gpr_population * population, FILE * fp);
void gpr_functions_are_equal(gpr_function * f1,
							 gpr_function * f2, int * result);
void gpr_S_expression(gpr_function * f, FILE * fp);
int gpr_random_function(int * instruction_set, int no_of_instructions,
						unsigned int *random_seed);
float gpr_random_value(float min_value, float max_value,
					   unsigned int * random_seed);
void gpr_init_system(gpr_system * system,
					 int islands,
					 int population_per_island,
					 int registers,
					 int sensors,
					 int actuators,
					 int max_tree_depth,
					 float min_value, float max_value,
					 int integers_only,
					 int ADFs,
					 int data_size, int data_fields,
					 unsigned int * random_seed,
					 int * instruction_set, int no_of_instructions);
void gpr_free_system(gpr_system * system);
void gpr_evaluate_system(gpr_system * system,
						 int time_steps, int reevaluate,
						 float (*evaluate_program)(int,gpr_function*,gpr_state*,int));
void gpr_generation_system(gpr_system * system,
						   int migration_interval,
						   float elitism,
						   int max_tree_depth,
						   float min_value, float max_value,
						   float mutation_prob,
						   float pure_mutant_prob,
						   int integers_only,
						   int ADFs,
						   int * instruction_set,
						   int no_of_instructions);
void gpr_sort_system(gpr_system * system);
float gpr_best_fitness_system(gpr_system * system);
gpr_function * gpr_best_individual_system(gpr_system * system);
void gpr_load_system(gpr_system * system,
					 FILE * fp,
					 int * instruction_set, int no_of_instructions);
void gpr_save_system(gpr_system *system, FILE * fp);
void gpr_arduino(gpr_function * f,
				 int baud_rate,
				 int digital_high,
				 int no_of_sensors, int no_of_actuators,
				 int no_of_registers,
				 int * digital_inputs, int no_of_digital_inputs,
				 int * analog_inputs, int no_of_analog_inputs,
				 int * digital_outputs, int no_of_digital_outputs,
				 int * analog_outputs, int no_of_analog_outputs,
				 int ADFs, FILE * fp);
void gpr_c_program(gpr_function * f,
				   int no_of_sensors, int no_of_actuators,
				   int no_of_registers,
				   int ADFs, FILE * fp);
int gpr_plot_history(gpr_population * population,
					 int history_type,
					 char * filename, char * title,
					 int image_width, int image_height);
int gpr_plot_history_system(gpr_system * sys,
							int history_type,
							char * filename, char * title,
							int image_width, int image_height);
int gpr_plot_fitness(gpr_population * population,
					 char * filename, char * title,
					 int image_width, int image_height);
int gpr_plot_fitness_system(gpr_system * sys,
							char * filename, char * title,
							int image_width, int image_height);
float gpr_median_fitness(gpr_population * population);
void gpr_clear_state(gpr_state * state);
void gpr_enforce_ADFs(gpr_function * f, gpr_state * state);
void gpr_xmlrpc_server(char * ruby_script_filename,
					   char * service_name, int port,
					   char * c_program_filename,
					   int sensors, int actuators);
int gpr_xmlrpc_client(char * service_name,
					  int port, char * hostname,
					  int sensors, int actuators,
					  float * sensor_value,
					  float * actuator_value);
int gpr_function_in_set(int function_type,
						int * instruction_set, int no_of_instructions);
void gpr_get_function_name(int function_type, float constant_value,
						   int blank_list,
						   char * name);
int gpr_mate_environment(gpr_environment * population,
						 int parent1_index, int parent2_index,
						 int max_tree_depth,
						 float min_value, float max_value,
						 float mutation_prob,
						 float pure_mutant_prob,
						 int integers_only, int ADFs,
						 int * instruction_set, int no_of_instructions);
void gpr_death(gpr_environment * population,
			   int victim_index);
void gpr_save_environment(gpr_environment *population, FILE * fp);
int gpr_load_environment(gpr_environment * population, FILE * fp);
int write_png_file(char * filename,
				   int width, int height,
				   unsigned char * buffer);

#endif

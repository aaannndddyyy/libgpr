/*
 libgpr - a library for genetic programming
 Copyright (C) 2013  Bob Mottram <bob@robotics.uk.to>

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
#include "libgpr/globals.h"
#include "libgpr/gprc.h"

#define RUN_STEPS 2

#undef ARDUINO

/* A test evaluation function.
   This tests how close the output is to
   the equation y = 3x^2 + 2x - 5 */
static float test_evaluate_program(int time_steps,
								   gprc_population * population,
								   int individual_index,
								   int custom_command)
{
	int itt,t,x,hits=0;
	float result,fitness=0, reference;
	float dropout_rate=0.1f;
	double diff=0,error;
	gprc_function * f = &population->individual[individual_index];

	/* clear the state */
	gprc_clear_state(f,
					 population->rows, population->columns,
					 population->sensors, population->actuators);

	/* for each time step */
	for (t = 0; t < time_steps; t++) {
		x = (t+1);
		/* sensor set to the current time step */
		gprc_set_sensor(f,0,x);
		for (itt = 0; itt < RUN_STEPS; itt++) {
			/* run the program */
			gprc_run(f, population, dropout_rate, 0.1f, 0);
		}
		/* observe the output */
		result = gprc_get_actuator(f,0,
								   population->rows,
								   population->columns,
								   population->sensors);
		/* the target value */
		reference = (3*x*x) + (2*x) - 5;

		if (fabs(reference) > 0.001f) {
			/* calculate the error as a percentage */
			error = (result - reference)/reference;
			diff += error*error;
			hits++;
		}
	}
	fitness = 100 - ((float)sqrt(diff/(float)hits)*100);
	return fitness;
}


static void test()
{
	int islands = 3;
	int migration_interval = 100;
	int population_per_island = 128;
	int rows = 9, columns = 16, sensors = 1, actuators = 1;
	int i, connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int chromosomes = 3;
	float min_value = -10, max_value = 10;
	gprc_system system;
	float elitism = 0.3f;
	float mutation_prob = 0.1f;
	int gen=0, time_steps = 20;
	int integers_only = 1;
	int use_crossover = 1;
	int modules = 2;
	FILE * fp;
	unsigned int random_seed = (unsigned int)time(NULL);
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;
	char * sensor_names[] = { "Sensor 0" };
	char * actuator_names[] = { "Actuator 0" };

#ifdef ARDUINO
	int no_of_digital_inputs=0;
	int no_of_analog_inputs=1;
	int no_of_digital_outputs=0;
	int no_of_analog_outputs=1;
	int analog_inputs[] = {1};
	int analog_outputs[] = {2};
#endif

	printf("Cartesian Genetic Programming test\n\n");

	/* create an instruction set */
	no_of_instructions =
		gprc_equation_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprc_init_system(&system, islands,
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

	while (gprc_best_fitness_system(&system) < 99) {
		/* evaluate each individual */
		gprc_evaluate_system(&system,
							 time_steps,0,
							 (*test_evaluate_program));

		/* produce the next generation */
		gprc_generation_system(&system,
							   migration_interval,
							   elitism,
							   mutation_prob,
							   use_crossover, &random_seed,
							   instruction_set, no_of_instructions);

		printf("Generation %05d  Fitness %.2f%% ",
			   gen, gprc_best_fitness_system(&system));
		for (i = 0; i < islands; i++) {
			printf("  %.5f",gprc_average_fitness(&system.island[i]));
		}
		printf("\n");

		if ((gen>0) && (gen % 20 == 0)) {
			gprc_draw_population("population.png",
								 640, 640, &system.island[0]);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_FITNESS,
									 "fitness.png", "Fitness History",
									 640, 480);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_FITNESS,
									 "fitness_average.png",
									 "Average Fitness History",
									 640, 480);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_DIVERSITY,
									 "diversity.png",
									 "Population Diversity History",
									 640, 480);

			gprc_plot_fitness_system(&system,
									 "fitness_histogram.png",
									 "Fitness Histogram",
									 640, 480);

			fp = fopen("fittest.dot","w");
			if (fp) {
				gprc_dot(gprc_best_individual_system(&system),
						 &system.island[0],
						 sensor_names,  actuator_names,
						 fp);
				fclose(fp);
			}

		}

		if (gprc_best_fitness_system(&system) > 99) break;
		gen++;
	}

	/* Save the best program to a file */
	fp = fopen("fittest.dat","w");
	if (fp) {
		gprc_save(gprc_best_individual_system(&system),
				  rows, columns,
				  connections_per_gene,
				  sensors, actuators,
				  data_size, data_fields,
				  fp);
		fclose(fp);
	}

	/* Save the population */
	fp = fopen("system.dat","w");
	if (fp) {
		gprc_save_system(&system, fp);
		fclose(fp);
	}

	/* save the best individual as a C program */
	fp = fopen("agent.c","w");
	if (fp) {
		gprc_c_program(&system,
					   gprc_best_individual_system(&system),
					   RUN_STEPS, 0, fp);
		fclose(fp);
	}

#ifdef ARDUINO
	gprc_arduino(&system,
				 gprc_best_individual_system(&system),
				 9600,255,
				 NULL, no_of_digital_inputs,
				 analog_inputs, no_of_analog_inputs,
				 NULL, no_of_digital_outputs,
				 analog_outputs, no_of_analog_outputs,
				 RUN_STEPS, 0, stdout);
#endif

	gprc_plot_history_system(&system,
							 GPR_HISTORY_FITNESS,
							 "fitness.png", "Fitness History",
							 640, 480);

	gprc_plot_history_system(&system,
							 GPR_HISTORY_AVERAGE,
							 "fitness_average.png",
							 "Average Fitness History",
							 640, 480);

	gprc_plot_history_system(&system,
							 GPR_HISTORY_DIVERSITY,
							 "diversity.png",
							 "Population Diversity History",
							 640, 480);

	gprc_plot_fitness(&system.island[0],
					  "fitness_histogram.png", "Fitness Histogram",
					  640, 480);


	gpr_xmlrpc_server("server.rb","libgpr",3573,
					  "agent",
					  sensors, actuators);


	/* free memory */
	gprc_free_system(&system);
}

int main(int argc, char* argv[])
{	
	test();
	return 1;
}


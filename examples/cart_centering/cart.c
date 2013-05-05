/*
 The cart centering problem, as described in
 Genetic Programming: On the programming of computers by
 means of natural selection, chapter 7
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
#include "libgpr/globals.h"
#include "libgpr/gprc.h"

#define RUN_STEPS 2

static float evaluate_controller(int time_steps,
								 gprc_population * population,
								 int individual_index,
								 int custom_command)
{
	int t, itt;
	float fitness=0,fitness_posn;
	gprc_function * f = &population->individual[individual_index];
	float velocity=0,accn=0,position=0,force=0,mass=20;
	float target_position = 55;
	float dropout_rate = 0.2f;

	/* clear the state */
	gprc_clear_state(f,
					 population->rows, population->columns,
					 population->sensors, population->actuators);

	/* for each time step */
	for (t = 0; t < time_steps; t++) {
		/* set the sensors */
		gprc_set_sensor(f,0,velocity);
		gprc_set_sensor(f,1,accn);
		gprc_set_sensor(f,2,position);
		gprc_set_sensor(f,3,force);
		/* run the program */
		for (itt = 0; itt < RUN_STEPS; itt++) {
			gprc_run(f, population, dropout_rate, 0, 0);
		}
		/* get the force */
		force = gprc_get_actuator(f,0,
								  population->rows,
								  population->columns,
								  population->sensors);

		/* physics of the cart */
		accn += force/mass;
		velocity += accn;
		position += velocity;

	}

	fitness_posn =
		80 - (fabs(position - target_position)*80.0f/target_position);
	if (fitness_posn < 0) fitness_posn = 0;

	/* How close is the output to the target position?
	   The cart should be at zero velocity when at the target */
	fitness = 
		fitness_posn +
		(20.0f/(1.0f + (fabs(velocity)*fabs(velocity))));
	return fitness;
}


static void cart_centering()
{
	int islands = 4;
	int migration_interval = 200;
	int population_per_island = 64;
	int rows = 9, columns = 10, sensors = 4, actuators = 1;
	int i, connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int modules = 0;
	int chromosomes=3;
	float min_value = -5, max_value = 5;
	gprc_system system;
	float elitism = 0.3f;
	float mutation_prob = 0.2f;
	int gen=0, time_steps = 50;
	int integers_only = 0;
	int use_crossover = 1;
	unsigned int random_seed = (unsigned int)time(NULL);
	int instruction_set[64], no_of_instructions=0;
	int data_size=10, data_fields=2;
	FILE * fp;
	char * sensor_names[] = {
		"Velocity",
		"Acceleration",
		"Position",
		"Force"
	};
	char * actuator_names[] = {
		"Force"
	};

	printf("Cart Centering Genetic Programming test\n\n");

	/* create an instruction set */
	no_of_instructions =
		gprc_equation_dynamic_instruction_set((int*)instruction_set);
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

	while (gprc_best_fitness_system(&system) < 99.99) {
		/* evaluate each individual */
		gprc_evaluate_system(&system,
							 time_steps,0,
							 (*evaluate_controller));

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
			printf("  %.3f",gprc_average_fitness(&system.island[i]));
		}
		printf("\n");

		if (((gen % 50 == 0) &&
			 (gen>0)) || (gprc_best_fitness(&system.island[0]) > 99.99)) {
			gprc_draw_population("population.png",
								 640, 640, &system.island[0]);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_FITNESS,
									 "fitness.png",
									 "Cart Centering Performance",
									 640, 480);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_AVERAGE,
									 "fitness_average.png",
									 "Cart Centering Average Performance",
									 640, 480);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_DIVERSITY,
									 "diversity.png",
									 "Cart Centering Diversity",
									 640, 480);

			gprc_plot_fitness(&system.island[0],
							  "fitness_histogram.png",
							  "Cart Centering Fitness Histogram",
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

		if (gprc_best_fitness_system(&system) > 99.99) break;
		gen++;
	}

	/* free memory */
	gprc_free_system(&system);
}

int main(int argc, char* argv[])
{	
	cart_centering();
	return 1;
}


/*
 Artificial ant problem, as described in
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

#define MAP_DIMENSION  32

#define EMPTY  0
#define FOOD   255

#define RUN_STEPS 1

static float evaluate_controller(int time_steps,
								 gprc_population * population,
								 int individual_index,
								 int custom_command)
{
	int t,i,x,y,direction=0,ant_x=0,ant_y=0,eaten=0;
	unsigned int random_seed = 534;
	float fitness=0,smell,dist;
	float turn_left=0,turn_right=0,forward=0;
	gprc_function * f = &population->individual[individual_index];
	int food_location[256];
	int itt, hits, food_items = 89;
	float dropout_rate = 0.0f;

	/* populate the map with food */
	for (i = 0; i < food_items; i++) {
		x = rand_num(&random_seed)%MAP_DIMENSION;
		y = rand_num(&random_seed)%MAP_DIMENSION;
		food_location[i*2] = x;
		food_location[i*2+1] = y;
	}

	/* clear the state */
	gprc_clear_state(f,
					 population->rows, population->columns,
					 population->sensors, population->actuators);

	/* for each time step */
	for (t = 0; t < time_steps; t++) {
		/* intensity of food smell at this location */
		smell=0;
		hits=0;
		for (i = 0; i < food_items; i++) {
			if (food_location[i*2]>-1) {
				x = food_location[i*2];
				y = food_location[i*2 + 1];
				dist = ((x - ant_x)*(x - ant_x)) +
					((y - ant_y)*(y - ant_y));

				if (dist==0) {
					/* eat food */
					eaten++;
					food_location[i*2] = -1;
				}
				else {
					smell += 100.0f/(1.0f + dist);
					hits++;
				}
			}
		}
		if (hits>0) {
			gprc_set_sensor(f,0,smell/(float)hits);
		}
		else {
			gprc_set_sensor(f,0,0);
		}

		/* run the program */
		for (itt = 0; itt < RUN_STEPS; itt++) {
			gprc_run(f, population, dropout_rate, 0, 0);
		}

		/* get the actuators */
		turn_left = fabs(gprc_get_actuator(f,0,
										   population->rows,
										   population->columns,
										   population->sensors));
		turn_right = fabs(gprc_get_actuator(f,1,
											population->rows,
											population->columns,
											population->sensors));
		forward = fabs(gprc_get_actuator(f,2,
										 population->rows,
										 population->columns,
										 population->sensors));

		/* decide which action to do */
		if ((turn_right > forward) ||
			(turn_left > forward)) {
			/* turn left or right */
			if (turn_right>turn_left) {
				direction++;
				if (direction>3) direction=0;
			}
			else {
				direction--;
				if (direction<0) direction=3;
			}
		}
		else {
			/* move forwards in the current direction */
			switch(direction) {
			case 0: { /* north */
				ant_y--;
				if (ant_y < 0) ant_y=0;
				break;
			}
			case 1: { /* east */
				ant_x++;
				if (ant_x>=MAP_DIMENSION) ant_x=MAP_DIMENSION-1;
				break;
			}
			case 2: { /* south */
				ant_y++;
				if (ant_y>=MAP_DIMENSION) ant_y=MAP_DIMENSION-1;
				break;
			}
			case 3: { /* west */
				ant_x--;
				if (ant_x < 0) ant_x=0;
				break;
			}
			}
		}		
	}
	fitness = eaten*100.0f/food_items;
	return fitness;
}


static void artificial_ant()
{
	int islands = 2;
	int migration_interval = 100;
	int population_per_island = 64;
	int rows = 9, columns = 8, sensors = 1, actuators = 3;
	int i, connections_per_gene = GPRC_MAX_ADF_MODULE_SENSORS+1;
	int modules = 0;
	int chromosomes=3;
	float min_value = -5, max_value = 5;
	gprc_system system;
	float elitism = 0.2f;
	float mutation_prob = 0.3f;
	int gen=0, time_steps = 2000;
	int integers_only = 0;
	int use_crossover = 1;
	unsigned int random_seed = (unsigned int)time(NULL);
	FILE * fp;
	int instruction_set[64], no_of_instructions=0;

	char * actuator_names[] = {
		"Left turn",
		"Right turn",
		"Forward"
	};
	char * sensor_names[] = {
		"Food sensor"
	};

	printf("Artificial Ant Genetic Programming test\n\n");

	/* create an instruction set */
	no_of_instructions =
		gprc_dynamic_instruction_set((int*)instruction_set);
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
					 integers_only, &random_seed,
					 instruction_set, no_of_instructions);

	while (gprc_best_fitness_system(&system) < 99) {
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

		/* show the best fitness value calculated
		   from the test data set */
		printf("Generation %05d  Fitness %.2f%% ",
			   gen, gprc_best_fitness_system(&system));
		for (i = 0; i < islands; i++) {
			printf("  %.5f",gprc_average_fitness(&system.island[0]));
		}
		printf("\n");

		if (((gen % 50 == 0) && (gen>0)) ||
			(gprc_best_fitness_system(&system) > 99)) {
			gprc_plot_history_system(&system,
									 GPR_HISTORY_FITNESS,
									 "fitness.png",
									 "Artificial Ant Performance",
									 640, 480);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_AVERAGE,
									 "fitness_average.png",
									 "Artificial Ant Average Performance",
									 640, 480);

			gprc_plot_history_system(&system,
									 GPR_HISTORY_DIVERSITY,
									 "diversity.png",
									 "Artificial Ant Diversity",
									 640, 480);

			gprc_plot_fitness(&system.island[0],
							  "fitness_histogram.png",
							  "Artificial Ant Fitness Histogram",
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

	/* free memory */
	gprc_free_system(&system);
}

int main(int argc, char* argv[])
{	
	artificial_ant();
	return 1;
}


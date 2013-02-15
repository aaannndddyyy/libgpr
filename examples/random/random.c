/*
 Evolving a random number generator, as described in
 Genetic Programming: On the programming of computers by means of
 natural selection, chapter 14
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
#include "libgpr/gpr.h"

/* test for randomness */
static float evaluate_function(int time_steps,
							   gpr_function * f,
							   gpr_state * state,
							   int custom_command)
{
	int t,index,ctr=0;
	float random_seed = 0.6257f;
	float fitness=0,average,standard_deviation,error;
	unsigned int histogram[16];
	unsigned int sequence[4];

	/* clear the histogram */
	for (index = 0; index < 16; index++) {
		histogram[index]=0;
	}

	/* clear the state */
	gpr_clear_state(state);

	/* for each time step */
	for (t = 0; t < time_steps; t++) {
		/* sensor set to the random seed */
		gpr_set_sensor(state,0,random_seed);
		/* run the program */
		gpr_run(f, state, 0);
		/* output random number */
		random_seed = gpr_get_actuator(state,0);
		/* function outputs a single bit */
		if (random_seed>0) {
			sequence[ctr++]=1;
		}
		else {
			sequence[ctr++]=0;
		}
		if (ctr > 3) {
			index = 0;
			if (sequence[0]==0) {
				if (sequence[1]==0) {
					if (sequence[2]==0) {
						if (sequence[3]==0) {
							index = 0;
						}
						else {
							index = 1;
						}				
					}
					else {
						if (sequence[3]==0) {
							index = 2;
						}
						else {
							index = 3;
						}				
					}
				}
				else {
					if (sequence[2]==0) {
						if (sequence[3]==0) {
							index = 4;
						}
						else {
							index = 5;
						}				
					}
					else {
						if (sequence[3]==0) {
							index = 6;
						}
						else {
							index = 7;
						}				
					}
				}				
			}
			else {
				if (sequence[1]==0) {
					if (sequence[2]==0) {
						if (sequence[3]==0) {
							index = 8;
						}
						else {
							index = 9;
						}				
					}
					else {
						if (sequence[3]==0) {
							index = 10;
						}
						else {
							index = 11;
						}				
					}
				}
				else {
					if (sequence[2]==0) {
						if (sequence[3]==0) {
							index = 12;
						}
						else {
							index = 13;
						}				
					}
					else {
						if (sequence[3]==0) {
							index = 14;
						}
						else {
							index = 15;
						}				
					}
				}				
			}
			histogram[index]++;
			ctr=0;
		}
	}

	average=0;
	for (index = 0; index < 16; index++) {
		average += (float)histogram[index];
	}
	average /= 16.0f;

	standard_deviation = 0;
	if (average > 0) {
		for (index = 0; index < 16; index++) {
			error = ((float)histogram[index] - average)/average;
			standard_deviation += error*error;
		}
	}
	standard_deviation = (float)sqrt(standard_deviation / 16.0f)*100;

	fitness = 100.0f - standard_deviation;
	/*if (fitness<0) fitness=0;*/
	return fitness;
}

static void random_number_generator()
{
	int islands = 4;
	int migration_interval = 100;
	int population_per_island = 64;
	int i, gen=0, max_depth = 8;
	gpr_system system;
	float min_value = -1;
	float max_value = 1;
	float elitism = 0.4f;
	float mutation_prob = 0.3f;
	float pure_mutant_prob = 0.2f;
	int time_steps = 100;
	FILE * fp;
	unsigned int random_seed = 8352;
	int sensors=1, actuators=1, registers=2;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;

	/* create an instruction set */
	instruction_set[0] = GPR_FUNCTION_ADD;
	instruction_set[1] = GPR_FUNCTION_SUBTRACT;
	instruction_set[2] = GPR_FUNCTION_MULTIPLY;
	instruction_set[3] = GPR_FUNCTION_DIVIDE;
	instruction_set[4] = GPR_FUNCTION_MODULUS;
	instruction_set[5] = GPR_FUNCTION_SET;
	instruction_set[6] = GPR_FUNCTION_GET;
	no_of_instructions = 7;

	/* create a population */
	gpr_init_system(&system, islands, population_per_island, registers, sensors, actuators,
					max_depth, min_value, max_value,
					integers_only, ADFs, &random_seed,
					(int*)instruction_set,no_of_instructions);

	while (gpr_best_fitness_system(&system) < 99) {

		/* evaluate each individual */
		gpr_evaluate_system(&system,
							time_steps,0,
							(*evaluate_function));

		/* produce the next generation */
		gpr_generation_system(&system,
							  migration_interval,
							  elitism,
							  max_depth,
							  min_value, max_value,
							  mutation_prob, pure_mutant_prob,
							  integers_only, ADFs,
							  (int*)instruction_set,no_of_instructions);

		/* show the best fitness value */
		printf("Generation %05d  Fitness %.2f%% ",gen, gpr_best_fitness_system(&system));
		for (i = 0; i < islands; i++) {
			printf("  %.5f",gpr_average_fitness(&system.island[i]));
		}
		printf("\n");

		if (((gen % 100 == 0) && (gen>0)) || (gpr_best_fitness(&system.island[0]) > 99)) {
			gpr_plot_history_system(&system,
									GPR_HISTORY_FITNESS,
									"fitness.png", "Random Number Generator Performance",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_AVERAGE,
									"fitness_average.png", "Random Number Generator Average Performance",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_DIVERSITY,
									"diversity.png", "Random Number Generator Diversity",
									640, 480);

			gpr_plot_fitness(&system.island[0],
							 "fitness_histogram.png", "Random Number Generator Fitness Histogram",
							 640, 480);

			fp = fopen("fittest.dot","w");
			if (fp) {
				gpr_dot(gpr_best_individual_system(&system), fp);
				fclose(fp);
			}
		}

		if (gpr_best_fitness_system(&system) > 99) break;
		gen++;
	}

	/* free memory */
	gpr_free_system(&system);
}

int main(int argc, char* argv[])
{	
	random_number_generator();
	return 1;
}


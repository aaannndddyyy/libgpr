/*
 Symbolic regression example from
 Genetic Programming: On the programming of computers by means of
 natural selection, chapter 10
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

/* test for Cos(2x) = 1 - 2Sin^2(x) */
static float evaluate_function(int time_steps,
							   gpr_function * f,
							   gpr_state * state,
							   int custom_command)
{
	int t;
	unsigned int random_seed = 527;
	float fitness=0, x, y;

	/* clear the state */
	gpr_clear_state(state);

	/* for each time step */
	for (t = 0; t < time_steps; t++) {
		x = rand_num(&random_seed)*2*3.1415927f;
		y = (float)cos(2*x);
		/* sensor set to the current time step */
		gpr_set_sensor(state,0,x);
		/* run the program */
		gpr_run(f, state, 0);
		/* how close is the output to the target equation? */
		fitness += 100-(fabs(gpr_get_actuator(state,0) - y)*50);
	}
	fitness /= (float)time_steps;
	return fitness;
}

static void test()
{
	int islands = 4;
	int migration_interval = 250;
	int population_per_island = 256;
	int i, gen, max_depth = 5;
	gpr_system system;
	float min_value = -10;
	float max_value = 10;
	float elitism = 0.2f;
	float mutation_prob = 0.4f;
	float pure_mutant_prob = 0.1f;
	int time_steps = 10;
	FILE * fp;
	unsigned int random_seed = 8352;
	int sensors=1, actuators=1, registers=0;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;

	/* create an instruction set */
	instruction_set[0] = GPR_FUNCTION_ADD;
	instruction_set[1] = GPR_FUNCTION_SUBTRACT;
	instruction_set[2] = GPR_FUNCTION_MULTIPLY;
	instruction_set[3] = GPR_FUNCTION_DIVIDE;
	instruction_set[4] = GPR_FUNCTION_SINE;
	instruction_set[5] = GPR_FUNCTION_SET;
	instruction_set[6] = GPR_FUNCTION_GET;
	no_of_instructions = 7;

	/* create a population */
	gpr_init_system(&system, islands, population_per_island, registers, sensors, actuators,
					max_depth, min_value, max_value,
					integers_only, ADFs, &random_seed,
					(int*)instruction_set,no_of_instructions);
	
	for (gen = 0; gen < 10000; gen++) {

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
									"fitness.png", "Symbolic Regression Performance",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_AVERAGE,
									"fitness_average.png", "Symbolic Regression Average Performance",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_DIVERSITY,
									"diversity.png", "Symbolic Regression Diversity",
									640, 480);

			gpr_plot_fitness(&system.island[0],
							 "fitness_histogram.png", "Symbolic Regression Fitness Histogram",
							 640, 480);
		}

		if (gpr_best_fitness_system(&system) > 99) break;
	}

	/* Save the best program in dot format */
	fp = fopen("fittest.dot","w");
	if (fp) {
		gpr_dot(gpr_best_individual_system(&system), fp);
		fclose(fp);
	}

	/* free memory */
	gpr_free_system(&system);
}

int main(int argc, char* argv[])
{	
	test();
	return 1;
}


/*
 The pursuer-evader problem as described in
 Genetic Programming: On the programming of computers by means of
 natural selection, chapter 15
 Copyright (C) 2012  Bob Mottram <bob@robotics.uk.to>

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
	const int no_of_trials=4;
	int t,r,dx,dy,trial,itt;
	unsigned int random_seed = (unsigned int)time(NULL);
	float fitness=0;
	int pursuer_x=10, pursuer_y=10;
	int pursuer_direction=0;
	int evader_x=0, evader_y=0;
	int evader_direction=0;
	const int minimum_capture_steps = 5;
	float turn_left,turn_right,forward;

	for (trial = 0; trial < no_of_trials; trial++) {

		pursuer_direction = rand_num(&random_seed)%8;
		evader_direction = rand_num(&random_seed)%8;

		/* clear the state */
		gpr_clear_state(state);

		/* randomly position the pursuer */
		itt=0;
		pursuer_x = rand_num(&random_seed)%20;
		pursuer_y = rand_num(&random_seed)%20;
		evader_x = rand_num(&random_seed)%20;
		evader_y = rand_num(&random_seed)%20;
		while (((abs(pursuer_x-evader_x)<10) &&
				((abs(pursuer_y-evader_y)<10)))) {
			pursuer_x = rand_num(&random_seed)%20;
			pursuer_y = rand_num(&random_seed)%20;
			evader_x = rand_num(&random_seed)%20;
			evader_y = rand_num(&random_seed)%20;
			itt++;
			if (itt>1000) {
				pursuer_x = 10;
				evader_x = 0;
				break;
			}
		}

		/* for each time step */
		for (t = 0; t < time_steps; t++) {
			/* update the evader */
			r = rand_num(&random_seed)%5;
			if (r==0) {
				/* turn right */
				evader_direction++;
				if (evader_direction>7) evader_direction=0;
			}
			if (r==1) {
				/* turn left */
				evader_direction--;
				if (evader_direction<0) evader_direction=7;
			}
			switch(evader_direction) {
			case 0: { /* north */
				evader_y--;
				break;
			}
			case 1: { /* north east */
				evader_y--;
				evader_x++;
				break;
			}
			case 2: { /* east */
				evader_x++;
				break;
			}
			case 3: { /* south east */
				evader_y++;
				evader_x++;
				break;
			}
			case 4: { /* south */
				evader_y++;
				break;
			}
			case 5: { /* south west */
				evader_y++;
				evader_x--;
				break;
			}
			case 6: { /* west */
				evader_x--;
				break;
			}
			case 7: { /* north west */
				evader_x--;
				evader_y--;
				break;
			}
			}

			/* sensor set to the random seed */
			dx = evader_x - pursuer_x;
			dy = evader_y - pursuer_y;
			gpr_set_sensor(state,0,(float)dx);
			gpr_set_sensor(state,1,(float)dy);
			gpr_set_sensor(state,2,((float)(rand_num(&random_seed)%10000)/5000.0f)-1.0f);
			/* run the program */
			gpr_run(f, state, 0);
			/* output random number */
			turn_left = gpr_get_actuator(state,0);
			turn_right = gpr_get_actuator(state,1);
			forward = gpr_get_actuator(state,2);

			if ((turn_left>forward) ||
				(turn_right>forward)) {
				if (turn_right>turn_left) {
					/* turn right */
					pursuer_direction++;
					if (pursuer_direction>7) pursuer_direction=0;
				}
				else {
					/* turn left */
					pursuer_direction--;
					if (pursuer_direction<0) pursuer_direction=7;
				}
			}

			switch(pursuer_direction) {
			case 0: { /* north */
				pursuer_y--;
				break;
			}
			case 1: { /* north east */
				pursuer_y--;
				pursuer_x++;
				break;
			}
			case 2: { /* east */
				pursuer_x++;
				break;
			}
			case 3: { /* south east */
				pursuer_y++;
				pursuer_x++;
				break;
			}
			case 4: { /* south */
				pursuer_y++;
				break;
			}
			case 5: { /* south west */
				pursuer_y++;
				pursuer_x--;
				break;
			}
			case 6: { /* west */
				pursuer_x--;
				break;
			}
			case 7: { /* north west */
				pursuer_x--;
				pursuer_y--;
				break;
			}
			}

			/* if the evader was caught */
			if ((evader_x==pursuer_x) &&
				(evader_y==pursuer_y)) {
				fitness += (100.0f - ((t-minimum_capture_steps)*100.0f/
									  (time_steps-minimum_capture_steps)));
				break;
			}
		}
	}
	return fitness/(float)no_of_trials;
}

static void pursuer_evader()
{
	int islands = 2;
	int migration_interval = 100;
	int population_per_island = 64;
	int i, gen=0, max_depth = 10;
	gpr_system system;
	float min_value = -1;
	float max_value = 1;
	float elitism = 0.2f;
	float mutation_prob = 0.1f;
	float pure_mutant_prob = 0.1f;
	int time_steps = 100;
	FILE * fp;
	unsigned int random_seed = 8352;
	int sensors=3, actuators=3, registers=2;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;

	/* create an instruction set */
	instruction_set[0] = GPR_FUNCTION_ADD;
	instruction_set[1] = GPR_FUNCTION_SUBTRACT;
	instruction_set[2] = GPR_FUNCTION_MULTIPLY;
	instruction_set[3] = GPR_FUNCTION_DIVIDE;
	instruction_set[4] = GPR_FUNCTION_EXP;
	instruction_set[5] = GPR_FUNCTION_LESS_THAN;
	instruction_set[6] = GPR_FUNCTION_SET;
	instruction_set[7] = GPR_FUNCTION_GET;
	no_of_instructions = 8;

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
									"fitness.png", "Pursuer Performance",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_AVERAGE,
									"fitness_average.png", "Pursuer Average Performance",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_DIVERSITY,
									"diversity.png", "Pursuer Diversity",
									640, 480);

			gpr_plot_fitness(&system.island[0],
							 "fitness_histogram.png", "Pursuer Fitness Histogram",
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
	pursuer_evader();
	return 1;
}


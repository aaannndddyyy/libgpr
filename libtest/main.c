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
#include "libgpr/gpr.h"

#undef ARDUINO

/* A test evaluation function.
   This tests how close the first output is to
   the equation y = 3x^2 + 2x - 5
   and the second output to y = 2x^2 - cos(x) + 8 */
static float test_evaluate_program(int time_steps,
								   gpr_function * f,
								   gpr_state * state,
								   int custom_command)
{
	int t,x;
	float result,fitness=0, reference;

	/* for each time step */
	for (t = 0; t < time_steps; t++) {
		x = t+1;
		/* sensor set to the current time step */
		gpr_set_sensor(state,0,x);
		/* run the program */
		gpr_run(f, state, 0);
		/* observe the first output */
		result = gpr_get_actuator(state,0);
		/* the target value */
		reference = (3*x*x) + (2*x) - 5;
		/* how close is the output to the target equation? */
		fitness += 1000-fabs(result - reference);


		/* observe the second output */
		result = gpr_get_actuator(state,1);
		/* the target value */
		reference = (2*x*x) - (float)cos(x) + 8;
		/* how close is the output to the target equation? */
		fitness += 1000-fabs(result - reference);
	}
	fitness /= (float)(time_steps*20);
	return fitness;
}

static void test()
{
	int islands = 4;
	int migration_interval = 100;
	int population_per_island = 128;
	int i, gen=0, max_depth = 8;
	gpr_system system;
	float min_value = -10;
	float max_value = 10;
	float elitism = 0.3f;
	float mutation_prob = 0.4f;
	float pure_mutant_prob = 0.2f;
	int time_steps = 20;
	FILE * fp;
	unsigned int random_seed = (unsigned int)time(NULL);
	int sensors=1, actuators=2, registers=4;
	int integers_only = 0;
	int ADFs = 1;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 0, data_fields = 0;
#ifdef ARDUINO
	int no_of_digital_inputs=0;
	int no_of_analog_inputs=sensors;
	int no_of_digital_outputs=0;
	int no_of_analog_outputs=actuators;
	int analog_inputs[] = {1};
	int analog_outputs[] = {2,3};
#endif

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gpr_init_system(&system, islands, population_per_island,
					registers, sensors, actuators,
					max_depth, min_value, max_value,
					integers_only, ADFs,
					data_size, data_fields,
					&random_seed,
					(int*)instruction_set,no_of_instructions);

	while (gpr_best_fitness_system(&system) < 99) {

		/* evaluate each individual */
		gpr_evaluate_system(&system,
							time_steps,0,
							(*test_evaluate_program));

		/* produce the next generation */
		gpr_generation_system(&system,
							  migration_interval,
							  elitism,
							  max_depth,
							  min_value, max_value,
							  mutation_prob, pure_mutant_prob,
							  integers_only,
							  ADFs,
							  (int*)instruction_set,no_of_instructions);

		/* show the best fitness value */
		printf("Generation %05d  Fitness %.2f%% ",
			   gen, gpr_best_fitness_system(&system));
		for (i = 0; i < islands; i++) {
			printf("  %.5f",gpr_average_fitness(&system.island[i]));
		}
		printf("\n");

		if ((gen>0) && (gen % 200 == 0)) {
			gpr_plot_history_system(&system,
									GPR_HISTORY_FITNESS,
									"fitness.png", "Fitness History",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_FITNESS,
									"fitness_average.png",
									"Average Fitness History",
									640, 480);

			gpr_plot_history_system(&system,
									GPR_HISTORY_DIVERSITY,
									"diversity.png",
									"Population Diversity History",
									640, 480);

			gpr_plot_fitness(&system.island[0],
							 "fitness_histogram.png",
							 "Fitness Histogram",
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

	/* Save the best program in dot format */
	fp = fopen("fittest.dot","w");
	if (fp) {
		gpr_dot(gpr_best_individual_system(&system), fp);
		fclose(fp);
	}
	/* Save the best program as an S-expression*/
	fp = fopen("fittest_s.dat","w");
	if (fp) {
		gpr_S_expression(gpr_best_individual_system(&system), fp);
		fclose(fp);
	}

	/* Save the best program to a file */
	fp = fopen("fittest.dat","w");
	if (fp) {
		gpr_save(gpr_best_individual_system(&system), fp);
		fclose(fp);
	}

	/* Save the system */
	fp = fopen("system.dat","w");
	if (fp) {
		gpr_save_system(&system, fp);
		fclose(fp);
	}

	fp = fopen("agent.c","w");
	if (fp) {
		gpr_c_program(gpr_best_individual_system(&system),
					  sensors, actuators, registers,
					  ADFs, fp);
		fclose(fp);
	}

#ifdef ARDUINO
	fp = fopen("arduino.txt","w");
	if (fp) {
		gpr_arduino(gpr_best_individual_system(&system),
					9600,255, sensors, actuators, registers,
					NULL, no_of_digital_inputs,
					analog_inputs, no_of_analog_inputs,
					NULL, no_of_digital_outputs,
					analog_outputs, no_of_analog_outputs,
					ADFs, 4, fp);
		fclose(fp);
	}
#endif

	gpr_plot_history_system(&system,
							GPR_HISTORY_FITNESS,
							"fitness.png", "Fitness History",
							640, 480);

	gpr_plot_history_system(&system,
							GPR_HISTORY_FITNESS,
							"fitness_average.png",
							"Average Fitness History",
							640, 480);

	gpr_plot_history_system(&system,
							GPR_HISTORY_DIVERSITY,
							"diversity.png",
							"Population Diversity History",
							640, 480);

	gpr_plot_fitness_system(&system,
							"fitness_histogram.png",
							"Fitness Histogram",
							640, 480);

	gpr_xmlrpc_server("server.rb","libgpr",3573,
					  "agent",
					  sensors, actuators);

	/* free memory */
	gpr_free_system(&system);
}

int main(int argc, char* argv[])
{	
	test();
	return 1;
}


/*
  Generating art using genetic programming.

  This generates a number of examples and the user can then choose
  the number of the image from which to produce the next gebneration.
  It would be better to do this with a GUI, but for maximum portability
  it's implemented as a command line application.
  Run the program in a shell and show the file manager along side it
  so that you can see the images which are produced.

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

#include <stdio.h>
#include <time.h>
#include "libgpr/globals.h"
#include "libgpr/gprcm.h"
#include "libgpr/pnglite.h"

#define RUN_STEPS  1

/* an individual produces an artwork */
void produce_art(int index,
				 gprcm_function * f,
				 gprcm_population * pop,
				 unsigned char * img,
				 int img_width, int img_height,
				 float dropout_rate,
				 char * filename)
{
	int x,y,n=0,itt,c;

	gprcm_clear_state(f, pop->rows,
					  pop->columns,
					  pop->sensors,
					  pop->actuators);

	/* for every image pixel */
	for (y = 0; y < img_height; y++) {
		for (x = 0; x < img_width; x++, n+=3) {
			gprcm_set_sensor(f, 0,
							 (x*2/(float)img_width)-0.5f);

			gprcm_set_sensor(f, 1,
							 (y*2/(float)img_height)-0.5f);

			if (n > 0) {
				gprcm_set_sensor(f, 2, img[n - 3]);
				gprcm_set_sensor(f, 3, img[n + 1 - 3]);
				gprcm_set_sensor(f, 4, img[n + 2 - 3]);
			}

			for (itt = 0; itt < RUN_STEPS; itt++) {
				/* run the program */
				gprcm_run(f, pop, dropout_rate, 0, 0);
			}

			for (c = 0; c < 3; c++) {
				img[n + c] =
					(unsigned char)(fmod(fabs(gprcm_get_actuator(f, c,
																 pop->rows,
																 pop->columns,
																 pop->sensors)),1.0f)*255);
			}
		}
	}
	write_png_file(filename, img_width, img_height, img);
}

/* produce an artwork for each individual in the population */
void produce_population_art(gprcm_population * pop,
							unsigned char * img,
							int img_width, int img_height,
							float dropout_rate)
{
	int index;
	gprcm_function * f;
	char filename[256];

	for (index = 0; index < pop->size; index++) {
		f = &pop->individual[index];
		sprintf(filename,"pic_%d.png",index);
		produce_art(index, f, pop, img,
					img_width, img_height,
					dropout_rate,
					filename);
		printf(".");
		fflush(stdout);
	}
}


static void art()
{
	int islands = 1;
	int migration_interval = 200;
	int population_per_island = 12;
	int rows = 4, columns = 10;
	int i, gen=0;
	int connections_per_gene = 3;
	int modules = 0;
	int chromosomes=1;
	gprcm_system sys;
	float min_value = -100;
	float max_value = 100;
	float elitism = 1.0f / population_per_island;
	float mutation_prob = 0.4f;
	int use_crossover = 1;
	unsigned int random_seed = (unsigned int)time(NULL);
	int sensors=5, actuators=3;
	int integers_only = 0;
	FILE *fp;
	char compile_command[256];
	int instruction_set[64], no_of_instructions=0;
	char * sensor_names[] = {
		"X",
		"Y",
		"Red",
		"Green",
		"Blue"
	};
	char * actuator_names[] = { "Red", "Green", "Blue" };
	int picked=0;
	unsigned char * img;
	int img_width=800, img_height=800;
	float dropout_rate = 0.1f;
	char str[256];
 
	img = (unsigned char*)malloc(img_width*img_height*3);

	/* create an instruction set */
	no_of_instructions = gprcm_dynamic_instruction_set(instruction_set);

	/*
	instruction_set[0] = GPR_FUNCTION_ADD;
	instruction_set[1] = GPR_FUNCTION_SUBTRACT;
	instruction_set[2] = GPR_FUNCTION_MULTIPLY;
	instruction_set[3] = GPR_FUNCTION_SINE;
	instruction_set[4] = GPR_FUNCTION_ARCSINE;
	instruction_set[5] = GPR_FUNCTION_COSINE;
	instruction_set[6] = GPR_FUNCTION_ARCCOSINE;
	instruction_set[7] = GPR_FUNCTION_SQUARE_ROOT;
	instruction_set[8] = GPR_FUNCTION_HEBBIAN;
	no_of_instructions = 9;
	*/

	/* create a population */
	gprcm_init_system(&sys, islands,
					  population_per_island,
					  rows, columns,
					  sensors, actuators,
					  connections_per_gene,
					  modules,
					  chromosomes,
					  min_value, max_value,
					  integers_only, &random_seed,
					  instruction_set, no_of_instructions);

	gpr_xmlrpc_server("server.rb","art",3573,
					  "./agent",
					  sensors, actuators);

	while (1) {

		printf("Generating artworks");

		/* produce the art */
		produce_population_art(&sys.island[0], img,
							   img_width, img_height,
							   dropout_rate);
		printf("\n");

		/* Pick the prefered individual */
		picked = -1;
		while ((picked<0) || (picked>= population_per_island)) {
			printf("Enter the number of the best image (0-%d): ",
				   population_per_island-1);
			fflush(stdout);
			gets(str);
			fflush(stdout);
			picked = atoi(str);
			printf("Picked: %d\n",picked);
		}
		for (i = 0; i < population_per_island; i++) {
			if (i != picked) {
				(&sys.island[0])->fitness[i] =
					rand_num(&random_seed)%10000/100000.0f;
			}
			else {
				(&sys.island[0])->fitness[i] = 99;
			}
		}

		/* produce the next generation */
		gprcm_generation_system(&sys,
								migration_interval,
								elitism,
								mutation_prob,
								use_crossover, &random_seed,
								instruction_set, no_of_instructions);


		fp = fopen("agent.c","w");
		if (fp) {
			/* save the best program */
			gprcm_c_program(&sys,
							gprcm_best_individual_system(&sys),
							RUN_STEPS, 0, fp);
			fclose(fp);

			/* compile the program */
			sprintf(compile_command,
					"gcc -Wall -std=c99 -pedantic -o agent agent.c -lm");
			assert(system(compile_command)==0);
		}

		fp = fopen("best.dot","w");
		if (fp) {
			gprcm_dot(gprcm_best_individual_system(&sys),
					  &sys.island[0],
					  sensor_names,  actuator_names,
					  fp);
			fclose(fp);
		}

		gen++;
	}

	free(img);

	/* free memory */
	gprcm_free_system(&sys);
}

int main(int argc, char* argv[])
{	
	art();
	return 1;
}


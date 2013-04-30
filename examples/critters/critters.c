/*
  This is a simple artificial life demo in which "critters"
  move around on a 2D toroidal flatland with randomly distributed
  food particles. They're able to smell food and also listen to
  and emit an audio waveform comprising of three fundamental
  frequencies.

  Critters can be male or female, and are only able to identify
  each other by the magnitude and type of sound.  This encourages
  an ecology into which "singing" is a signifier of sex and
  potential mate fitness.

  When a female reproduces she gives some of her energy to the child.
  The costly nature of reproduction for the female is intended to make
  her more choosy about selecting mates (via their song).

  Males can also provision females during mating - passing some of
  their energy to the female.  This helps "pay" for the energy cost
  of reproduction, and means that there is a tradeoff between
  survival of the individual and survival of its future offsprint.

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
#include "libgpr/gprc.h"

/* constants to indicate sex */
#define SEX_FEMALE  0
#define SEX_MALE    1

/* structure containing critter attributes */
struct crit {
	unsigned int energy;
	unsigned char sex;
	unsigned int age;
	float x, y, orientation;
	int mating_inhibition;
	float waveform;
	int father;
	int mother;
};
typedef struct crit critter;

/* the age at which critters can mate in time steps */
#define AGE_OF_MATURITY          100

/* the maximum age in time steps */
#define MAX_AGE                  100000

/* the number of time steps during which mating will be inhibited */
#define MATING_INHIBITION_STEPS  500

/* initial energy level */
#define INITIAL_ENERGY           100000

/* the amount of energy in food particles */
#define FOOD_ENERGY              1000

/* maximum energy level */
#define MAX_ENERGY               200000

/* Reproduction has an energy cost for females.
   This encorages them to be choosy about mates. */
#define REPRODUCTION_ENERGY      10000

/* used to indicate not speaking */
#define NO_WAVEFORM              9999

/* sensor types */
enum {
	SENSOR_SOUND1 = 0,
	SENSOR_SOUND2,
	SENSOR_SOUND3,
	SENSOR_SOUND4,
	SENSOR_SOUND5,
	SENSOR_FOOD,
	SENSOR_SMELL1,
	SENSOR_SMELL2,
	SENSOR_ORIENTATION,
	SENSOR_SEX
};

/* actuator types */
enum {
	ACTUATOR_TURN_LEFT = 0,
	ACTUATOR_TURN_RIGHT,
	ACTUATOR_SPEED,
	ACTUATOR_TURN_RATE,
	ACTUATOR_SPEAKING,
	ACTUATOR_SPEAK_FREQ1,
	ACTUATOR_SPEAK_FREQ2,
	ACTUATOR_SPEAK_FREQ3,
	ACTUATOR_MATE,
	ACTUATOR_PROVISIONING
};

/* function for audio volume decay with distance */
#define SOUND_VOLUME(v,dist) ((v)/(1.0f + (dist)))

/* returns the audio waveform for the given critter
   at the given time */
float critter_speak(gprc_function * critter_brain,
					float t,
					int rows, int columns, int sensors)
{
	float waveform;

	waveform =
		(float)sin(t + fmod(fabs(gprc_get_actuator(critter_brain,
												   ACTUATOR_SPEAK_FREQ1,
												   rows, columns, sensors)),
							3.1415927f*2)) +
		((float)sin(t + fmod(fabs(gprc_get_actuator(critter_brain,
													ACTUATOR_SPEAK_FREQ2,
													rows, columns, sensors)),
							 3.1415927f*2))*0.5f) +
		((float)sin(t + fmod(fabs(gprc_get_actuator(critter_brain,
													ACTUATOR_SPEAK_FREQ3,
													rows, columns, sensors)),
							 3.1415927f*2))*0.25f);	
	return waveform;
}

/* updates the waveform for each critter */
void update_waveforms(gprc_environment * population,
					  critter * being,
					  float t)
{
	int i;

	for (i = 0; i < population->population_size; i++) {
		if (gprc_get_actuator(&population->individual[i],
							  ACTUATOR_SPEAKING,
							  population->rows,
							  population->columns,
							  population->sensors) > 0) {
			being[i].waveform =
				critter_speak(&population->individual[i],
							  t, population->rows,
							  population->columns,
							  population->sensors);
		}
		else {
			being[i].waveform = NO_WAVEFORM;
		}
	}
}

/* returns the audio volume at a given location */
float sound_at_location(gprc_environment * population,
						float x, float y,
						critter * being,
						float t)
{
	int i;
	float x2,y2,dx,dy,dist;
	float w,waveform = 0;

	for (i = 0; i < population->population_size; i++) {
		if (being[i].waveform != NO_WAVEFORM) {
			x2 = being[i].x;
			y2 = being[i].y;
			dx = x - x2;
			dy = y - y2;
			dist = dx*dx + dy*dy;
			w = SOUND_VOLUME(being[i].waveform, dist);
			if (w > waveform) {
				waveform = w;
			}
		}
	}
	return waveform;
}

/* returns the food smell magnitude at the given location */
float smell_at_location(int critter_index,
						critter * being,
						float x, float y,
						float * food_location,
						int food_elements,
						float map_size,
						unsigned int * random_seed)
{
	int i;
	float x2,y2,dx,dy,dist;
	float s,smell = 0;

	for (i = 0; i < food_elements; i++) {
		x2 = food_location[i*2];
		y2 = food_location[i*2+1];
		dx = x - x2;
		dy = y - y2;
		dist = dx*dx + dy*dy;
		s = SOUND_VOLUME(1.0f, dist);
		if (s > smell) {
			smell = s;
		}
		if (dist < 1) {			
			food_location[i*2] =
				((rand_num(random_seed)%10000)/10000.0f)*map_size;
			food_location[i*2+1] =
				((rand_num(random_seed)%10000)/10000.0f)*map_size;
			being[critter_index].energy += FOOD_ENERGY;
			if (being[critter_index].energy > MAX_ENERGY) {
				being[critter_index].energy = MAX_ENERGY;
			}
			printf("Eat %d\n",critter_index);
		}
	}
	return smell;
}

/* updates the food smell sensors for the given critter */
void critter_smell(gprc_environment * population,
				   int critter_index,
				   critter * being,
				   float * food_location,
				   int food_elements,
				   float map_size)
{
	float smell;

	smell =
		smell_at_location(critter_index, being,
						  being[critter_index].x,
						  being[critter_index].y,
						  food_location, food_elements,
						  map_size,
						  &(&population->individual[critter_index])->random_seed);

	gprc_set_sensor(&population->individual[critter_index],
					SENSOR_SMELL2,
					gprc_get_sensor(&population->individual[critter_index],
									SENSOR_SMELL1));

	gprc_set_sensor(&population->individual[critter_index],
					SENSOR_SMELL1, smell);
}

/* updates the sound sensors for the given critter */
void critter_listen(gprc_environment * population,
					int critter_index,
					critter * being,
					float t)
{
	gprc_function * critter_brain;
	float x, y, waveform;

	critter_brain = &population->individual[critter_index];
	x = being[critter_index].x;
	y = being[critter_index].y;

	waveform =
		sound_at_location(population, x, y, being, t);

	/* the sound sensors here represent audio
	   over a number of time steps */
	gprc_set_sensor(critter_brain,
					SENSOR_SOUND5,
					gprc_get_sensor(critter_brain,
									SENSOR_SOUND4));
	gprc_set_sensor(critter_brain,
					SENSOR_SOUND4,
					gprc_get_sensor(critter_brain,
									SENSOR_SOUND3));
	gprc_set_sensor(critter_brain,
					SENSOR_SOUND3,
					gprc_get_sensor(critter_brain,
									SENSOR_SOUND2));
	gprc_set_sensor(critter_brain,
					SENSOR_SOUND2,
					gprc_get_sensor(critter_brain,
									SENSOR_SOUND1));

	gprc_set_sensor(critter_brain, SENSOR_SOUND1,
					waveform);
}

/* update the motion of the given critter */
void critter_motion(gprc_function * critter_brain,
					float * x, float * y,
					float * orientation,
					int rows, int columns, int sensors,
					float map_size)
{
	float turn_rate, speed;

	/* the rate of turn */
	turn_rate =
		fmod(fabs(gprc_get_actuator(critter_brain,ACTUATOR_TURN_RATE,
									rows, columns, sensors)),1.0f);
	/* the speed for forward motion */
	speed =
		0.2f +
		fmod(fabs(gprc_get_actuator(critter_brain,ACTUATOR_SPEED,
									rows, columns, sensors)),1.0f);

	/* turn to the right */
	if (gprc_get_actuator(critter_brain,ACTUATOR_TURN_RIGHT,
						  rows, columns, sensors) >
		gprc_get_actuator(critter_brain,ACTUATOR_TURN_LEFT,
						  rows, columns, sensors)+0.1f) {
		*orientation = *orientation + turn_rate;
	}

	/* turn to the left */
	if (gprc_get_actuator(critter_brain,ACTUATOR_TURN_LEFT,
						  rows, columns, sensors) >
		gprc_get_actuator(critter_brain,ACTUATOR_TURN_RIGHT,
						  rows, columns, sensors)+0.1f) {
		*orientation = *orientation - turn_rate;
	}

	/* keep orientation within range */
	if (*orientation < 0) {
		*orientation = *orientation + (2*3.1415927f);
	}
	if (*orientation > 2*3.1415927f) {
		*orientation = *orientation - (2*3.1415927f);
	}

	/* move */
	*x = *x + (float)(speed*sin(*orientation));
	*y = *y + (float)(speed*cos(*orientation));

	/* keep within the map coordinaets */
	if (*x < 0) *x = *x + map_size;
	if (*x >= map_size) *x = *x - map_size;
	if (*y < 0) *y = *y + map_size;
	if (*y >= map_size) *y = *y - map_size;
}

/* mate and produce offspring */
void critter_mate(int critter_index,
				  critter * being,
				  gprc_environment * population,
				  float mutation_prob,
				  int * instruction_set, int no_of_instructions,
				  unsigned int * random_seed)
{
	int i,child_index;
	float dx,dy,dist,mate1,mate2,provisioning;

	/* is the critter old enough? */
	if (being[critter_index].age < AGE_OF_MATURITY) return;
	/* is this critter female? */
	if (being[critter_index].sex != SEX_FEMALE) return;
	/* has the critter not mated for a while? */
	if (being[critter_index].mating_inhibition > 0) return;

	/* if this actuator is greater than zero then
	   the female critter is ready to mate */
	mate1 =
		gprc_get_actuator(&population->individual[critter_index],
						  ACTUATOR_MATE,
						  population->rows, population->columns,
						  population->sensors);
	if (mate1 < 0) return;

	for (i = 0; i < population->population_size; i++) {
		if (i != critter_index) {
			/* look for males who are old enough to mate
			   and also avoid mating with close relations */
			if ((being[i].sex == SEX_MALE) &&
				(being[i].age > AGE_OF_MATURITY) &&
				(being[i].mother != critter_index) &&
				(being[critter_index].father != i)) {

				/* calculate the distance between them */
				dx = being[i].x - being[critter_index].x;
				dy = being[i].y - being[critter_index].y;
				dist = dx*dx + dy*dy;
				/* within range to mate */
				if (dist < 1) {
					/* if this actuator is greater than zero then
					   the male critter is ready to mate */
					mate2 =
						gprc_get_actuator(&population->individual[i],
										  ACTUATOR_MATE,
										  population->rows,
										  population->columns,
										  population->sensors);
					if (mate2 > 0) {

						/* generate a new critter brain */
						child_index =
							gprc_mate_environment(population,
												  critter_index, i,
												  mutation_prob, 1,
												  instruction_set,
												  no_of_instructions);

						if (child_index > -1) {
							/* male provides provisioning */
							provisioning =
								fmod(fabs(gprc_get_actuator(&population->individual[i],
															ACTUATOR_PROVISIONING,
															population->rows,
															population->columns,
															population->sensors)),1.0f)*being[i].energy*0.5f;

							being[i].energy -=
								(unsigned int)provisioning;
							being[critter_index].energy +=
								(unsigned int)provisioning;

							/* reproduction cost for the female */
							if (being[critter_index].energy >
								REPRODUCTION_ENERGY) {
								being[critter_index].energy -=
									REPRODUCTION_ENERGY;
							}
							else {
								being[critter_index].energy = 0;
							}
							/* inhibit from mating for a while */
							being[critter_index].mating_inhibition =
								MATING_INHIBITION_STEPS;
							being[i].mating_inhibition =
								MATING_INHIBITION_STEPS;

							/* create a new being entry */
							if (rand_num(random_seed)%10000>5000) {
								being[child_index].sex = SEX_FEMALE;
							}
							else {
								being[child_index].sex = SEX_MALE;
							}
							being[child_index].age = 0;
							being[child_index].x =
								being[critter_index].x;
							being[child_index].y =
								being[critter_index].y;
							being[child_index].orientation =
								being[i].orientation;
							being[child_index].energy = INITIAL_ENERGY;
							being[child_index].mating_inhibition =
								AGE_OF_MATURITY;
							being[child_index].father = i;
							being[child_index].mother = critter_index;
							printf("Birth F%d M%d C%d\n",
								   critter_index, i,
								   child_index);
						}
					}					
				}
			}
		}
	}
}

/* the main program */
void critters()
{
	critter * being;
	int i, population_size = 128;
	int max_population_size = 256;
	int rows=9, columns=10, sensors=10, actuators=10;
	int connections_per_gene=8;
	int chromosomes=3;
	int modules=0;
	float min_value=-10, max_value=10;
	int integers_only=0;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;
	gprc_environment population;
	float mutation_prob = 0.2f;
	float * food_location;
	float map_size = 1000;
	float t, dropout_rate=0;
	gprc_function * critter_brain;
	int food_elements = 200;
	int data_size=0, data_fields=0;

	/* create an array to store critter properties */
	being = (critter*)malloc(max_population_size*
							 sizeof(critter));

	/* create an array to store food locations */
	food_location =
		(float*)malloc(food_elements*2*sizeof(float));

	/* random distribution of food */
	for (i = 0; i < food_elements; i++) {
		food_location[i*2] =
			((rand_num(&random_seed)%10000)/10000.0f)*map_size;
		food_location[i*2+1] =
			((rand_num(&random_seed)%10000)/10000.0f)*map_size;
	}

	/* randomly initialise critters */
	for (i = 0; i < population_size; i++) {
		/* random sex */
		if (rand_num(&random_seed)%10000 > 5000) {
			being[i].sex = SEX_FEMALE;
		}
		else {
			being[i].sex = SEX_MALE;
		}
		/* first generation have no parents */
		being[i].father=-1;
		being[i].mother=-1;
		/* no mating inhibition */
        being[i].mating_inhibition = 0;
		/* random age */
		being[i].age = rand_num(&random_seed)%MAX_AGE;
		being[i].energy = INITIAL_ENERGY;
		/* random map location */
		being[i].x =
			((rand_num(&random_seed)%10000)/10000.0f)*map_size;
		being[i].y =
			((rand_num(&random_seed)%10000)/10000.0f)*map_size;
		/* random orientation */
		being[i].orientation =
			((rand_num(&random_seed)%10000)/10000.0f)*2*3.1415927f;
	}

	/* create an instruction set */
	no_of_instructions =
		gprc_dynamic_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gprc_init_environment(&population,
						  max_population_size,
						  population_size,
						  rows, columns,
						  sensors, actuators,
						  connections_per_gene,
						  modules,
						  chromosomes,
						  min_value, max_value,
						  integers_only,
						  data_size, data_fields,
						  &random_seed,
						  instruction_set,
						  no_of_instructions);

	t = 0;
	while (population.population_size > 0) {
		/* increment the time */
		t += 0.1f;
		if (t > 999999) t = 0;

		/* update the audio for each critter */
		update_waveforms(&population, being, t);

		/* clear the number of matings */
		population.matings = 0;
		for (i = 0; i < population.population_size; i++) {
			critter_brain = &population.individual[i];

			/* listen for speech */
			critter_listen(&population, i, being, t);

			/* smell food */
			critter_smell(&population, i,
						  being,
						  food_location,
						  food_elements,
						  map_size);

			/* the critter knows its own orientation */
			gprc_set_sensor(&population.individual[i],
							SENSOR_ORIENTATION,
							being[i].orientation);

			/* the critter knows its own sex */
			gprc_set_sensor(&population.individual[i],
							SENSOR_SEX,
							being[i].sex);

			/* run the brain */
			gprc_run_environment(critter_brain,
								 &population,
								 dropout_rate, 0, 0);

			/* update the pose */
			critter_motion(critter_brain,
						   &being[i].x, &being[i].y,
						   &being[i].orientation,
						   rows, columns, sensors,
						   map_size);

			/* mating with others */
			critter_mate(i, being, &population,
						 mutation_prob,
						 instruction_set, no_of_instructions,
						 &random_seed);

			/* increment the age */
			being[i].age++;

			/* energy depletion */
			if (being[i].energy > 0) {
				being[i].energy--;
			}

			/* reduce mating inhibition */
			if (being[i].mating_inhibition > 0) {
				being[i].mating_inhibition--;
			}

			/* death */
			if ((being[i].energy <= 0) ||
				(being[i].age > MAX_AGE)) {
				being[i] = being[population.population_size];
				gprc_death(&population, i);
				printf("Population %d\n",
					   population.population_size);
				i--;
			}
		}
		if (population.population_size==0) break;
	}

	/* free the memory */
	free(food_location);
	free(being);
	gprc_free_environment(&population);

	printf("Finished\n");
}

int main(int argc, char* argv[])
{	
	critters();
	return 1;
}


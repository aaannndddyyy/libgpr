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

#include "som.h"

/* initialises a SOM structure */
void gpr_som_init(int dimension,
				  int no_of_sensors,
				  gpr_som * som)
{
	int i;
	gpr_som_element * elem;

	som->dimension = dimension;
	som->no_of_sensors = no_of_sensors;
	som->weight = (gpr_som_element**)malloc(dimension*dimension*
											sizeof(gpr_som_element*));
	assert(som->weight);
	for (i = 0; i < dimension*dimension; i++) {
		elem = (gpr_som_element *)&som->weight[i];
		
		elem->vect = (float*)malloc(no_of_sensors*sizeof(float));
	}
	som->state = (float*)malloc(dimension*dimension*sizeof(float));
	som->sensor_scale = (float*)malloc(no_of_sensors*sizeof(float));
}

/* frees a SOM structure */
void gpr_som_free(gpr_som * som)
{
	int i;
	gpr_som_element * elem;

	for (i = 0; i < som->dimension*som->dimension; i++) {
		elem = (gpr_som_element *)&som->weight[i];
		free(elem->vect);
	}
	free(som->weight);
	free(som->state);
	free(som->sensor_scale);
}

/* randomly initialise a sensor */
int gpr_som_init_sensor(gpr_som * som,
						int sensor_index,
						float min_value, float max_value,
						unsigned int * random_seed)
{
	int i;
	gpr_som_element * elem;

	if (max_value <= min_value) return -1;

	som->sensor_scale[sensor_index] = 1.0f / (max_value - min_value);

	for (i = 0; i < som->dimension*som->dimension; i++) {
		elem = (gpr_som_element *)&som->weight[i];
		elem->vect[sensor_index] =
			min_value +
			((max_value - min_value)*
			 (rand_num(random_seed)%10000)/10000.0f);
	}
	return 0;
}

/* Initialise weights for the given sensor using the given
   training data.  This ensures that there is an even
   distribution of values within the expected range */
int gpr_som_init_sensor_from_data(gpr_som * som,
								  int sensor_index,
								  int training_data_field_index,
								  float * training_data,
								  int training_data_fields_per_example,
								  int no_of_training_examples,
								  unsigned int * random_seed)
{
	float min_value=0, max_value=0, value;
	int i;

	/* get the minimum and maximum values for this field */
	for (i = 0; i < no_of_training_examples; i++) {
		value = training_data[i*training_data_fields_per_example +
							  training_data_field_index];
		if ((i == 0) || (value < min_value)) {
			min_value = value;
		}
		if ((i == 0) || (value > max_value)) {
			max_value = value;
		}
	}

	if (max_value <= min_value) return -1;

	/* initialise within the range of values */
	return gpr_som_init_sensor(som, sensor_index,
							   min_value, max_value,
							   random_seed);	
}

/* adjusts the weights within the SOM */
void gpr_som_learn(gpr_som * som,
				   float * sensors,
				   int inhibit_radius,
				   int excite_radius,
				   float learning_rate)
{
	int i,s,x=0,y=0,xx,yy,dx,dy,r;
	int inhibit_radius2;
	float max=0,diff;
	gpr_som_element * elem;

	/* location of the best response */
	for (i = 0; i < som->dimension*som->dimension; i++) {
		if (som->state[i] > max) {
			max = som->state[i];
			y = i / som->dimension;
			x = i - (y*som->dimension);
		}
	}

	excite_radius *= excite_radius;
	inhibit_radius2 = inhibit_radius*inhibit_radius;

	/* alter weights */
	for (xx = x - inhibit_radius; xx <= x + inhibit_radius; xx++) {
		if ((xx < 0) || (xx >= som->dimension)) continue;
		dx = xx - x;
		for (yy = y - inhibit_radius; yy <= y + inhibit_radius; yy++) {
			if ((yy < 0) || (yy >= som->dimension)) continue;
			dy = yy - y;
			elem =
				(gpr_som_element *)&som->weight[(yy*som->dimension) +
												xx];
			r = (dx*dx) + (dy*dy);
			if (r <= excite_radius) {
				/* excite */
				for (s = 0; s < som->no_of_sensors; s++) {
					if ((elem->vect[s] != GPR_MISSING_VALUE) &&
						(sensors[s] != GPR_MISSING_VALUE)) {
						diff = sensors[s] - elem->vect[s];
						elem->vect[s] += diff*learning_rate;
					}
				}
			}
			else {
				if (r <= inhibit_radius2) {
					/* inhibit */
					for (s = 0; s < som->no_of_sensors; s++) {
						if ((elem->vect[s] != GPR_MISSING_VALUE) &&
							(sensors[s] != GPR_MISSING_VALUE)) {
							diff = sensors[s] - elem->vect[s];
							elem->vect[s] -= diff*learning_rate;
						}
					}
				}
			}
		}
	}	
}

/* updates the SOM */
int gpr_som_update(float * sensors,
				   gpr_som * som,
				   float * x, float * y)
{
	int i,winner=-1;
	float diff,min=0,max=0;
	gpr_som_element * elem;

	/* clear the state */
	memset((void*)som->state,'\0',
		   som->dimension*som->dimension*sizeof(float));

	/* compare sensor values against weights.
	   The sensor_scale value is used to avoid bias
	   towards any particular sensor */
#pragma omp parallel for
	for (i = 0; i < som->dimension*som->dimension; i++) {
		elem = (gpr_som_element *)&som->weight[i];
		for (int s = 0; s < som->no_of_sensors; s++) {
			if ((elem->vect[s] != GPR_MISSING_VALUE) &&
				(sensors[s] != GPR_MISSING_VALUE)) {
				som->state[i] +=
					(sensors[s] - elem->vect[s])*
					(sensors[s] - elem->vect[s])*
					som->sensor_scale[s];
			}
		}
	}

	/* find the winner */
	for (i = 0; i < som->dimension*som->dimension; i++) {
		diff = som->state[i];
		if ((i == 0) || (diff < min)) {
			min = diff;
			winner = i;
		}
		if ((i == 0) || (diff > max)) {
			max = diff;			
		}
	}
	if (winner > -1) {
		/* return the peak location in the range 0.0 - 1.0 */
		*y = (int)(winner/som->dimension) / (float)som->dimension;
		*x = (int)(winner%som->dimension) / (float)som->dimension;
	}
	else {
		*x = -1;
		*y = -1;
	}

	/* normalise the state values in the range 0.0 - 1.0 */
	if (max > min) {
#pragma omp parallel for
		for (i = 0; i < som->dimension*som->dimension; i++) {
			som->state[i] = 1.0f - ((som->state[i] - min)/(max - min));
		}
	}
	return winner;
}

/* updates the SOM without altering any state */
void gpr_som_run(float * sensors,
				 gpr_som * som,
				 float * x, float * y)
{
	int s,i,winner=0;
	float diff,min=0,ds;
	gpr_som_element * elem;

	*x=0;
	*y=0;

	/* compare sensor values against weights.
	   The sensor_scale value is used to avoid bias
	   towards any particular sensor */
	for (i = 0; i < som->dimension*som->dimension; i++) {
		elem = (gpr_som_element *)&som->weight[i];
		diff=0;
		for (s = 0; s < som->no_of_sensors; s++) {
			if ((elem->vect[s] != GPR_MISSING_VALUE) &&
				(sensors[s] != GPR_MISSING_VALUE)) {
				ds = sensors[s] - elem->vect[s];
				diff += ds * ds * som->sensor_scale[s];
			}
		}
		if ((i==0) || (diff < min)) {
			min = diff;
			winner=i;
		}
	}
	*y = (int)(winner/som->dimension) / (float)som->dimension;
	*x = (int)(winner%som->dimension) / (float)som->dimension;
}

/* Updates an array containing outputs for the given
   training or test data. */
void gpr_som_outputs_from_data(gpr_som * som,
							   int * data_field_index,
							   float * data,
							   int fields_per_sample,
							   int no_of_samples,
							   float * result)
{
	int i,s;
	float x=0,y=0;
	float * sensors;

	sensors = (float*)malloc(som->no_of_sensors*sizeof(float));

	for (i = 0; i < no_of_samples; i++) {
		for (s = 0; s < som->no_of_sensors; s++) {
			sensors[s] = data[i*fields_per_sample +
							  data_field_index[s]];
		}
		gpr_som_run(sensors, som, &x, &y);
		result[i*2] = x;
		result[i*2 + 1] = y;
	}

	free(sensors);
}

/* learn from training data */
void gpr_som_learn_from_data(gpr_som * som,
							 int * data_field_index,
							 float * training_data,
							 int fields_per_sample,
							 int no_of_samples,
							 int learning_itterations,
							 int inhibit_radius, int excite_radius,
							 float learning_rate,
							 unsigned int * random_seed,
							 int show_progress)
{
	int i,j,n,s;
	float x=0,y=0;
	float * sensors;

	sensors = (float*)malloc(som->no_of_sensors*sizeof(float));

	for (i = 0; i < learning_itterations; i++) {
		for (j = 0; j < no_of_samples; j++) {
			/* pick a training sample */
			n = rand_num(random_seed)%no_of_samples;
			for (s = 0; s < som->no_of_sensors; s++) {
				sensors[s] = training_data[n*fields_per_sample +
										   data_field_index[s]];
			}
			gpr_som_update(sensors, som, &x, &y);
			gpr_som_learn(som, sensors,
						  inhibit_radius, excite_radius,
						  learning_rate);
		}
		if (show_progress > 0) {
			printf(".");
			fflush(stdout);
		}
	}

	free(sensors);
}

void gpr_som_save(gpr_som * som,
				  FILE * fp)
{
	int i,s;
	int max = som->dimension*som->dimension;
	gpr_som_element * elem;

	fprintf(fp,"%d\n",som->dimension);
	fprintf(fp,"%d\n",som->no_of_sensors);
	for (i = 0; i < max; i++) {
		elem = (gpr_som_element*)&som->weight[i];
		for (s = 0; s < som->no_of_sensors; s++) {
			fprintf(fp,"%.5f\n",elem->vect[s]);
		}
	}
}

void gpr_som_load(gpr_som * som,
				  FILE * fp)
{
	char line[256];
	int i,s,ctr=0,dimension=1,no_of_sensors=1;
	gpr_som_element * elem;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				switch(ctr) {
				case 0: {
					dimension = atoi(line);
					break;
				}
				case 1: {
					no_of_sensors = atoi(line);
					break;
				}
				}

				if (ctr==1) break;

				ctr++;
			}
		}
	}

	gpr_som_init(dimension, no_of_sensors, som);

	ctr=0;
	for (i = 0; i < dimension*dimension; i++) {
		elem = (gpr_som_element*)&som->weight[i];
		for (s = 0; s < no_of_sensors; s++) {
			if (fgets(line , 255 , fp) != NULL ) {
				if (strlen(line)>0) {
					elem->vect[s] = atof(line);
				}
			}	
		}
	}
}

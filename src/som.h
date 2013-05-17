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

#ifndef GPR_SOM_H
#define GPR_SOM_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "globals.h"
#include "gpr.h"

struct gpr_som_elem {
    float * vect;
};
typedef struct gpr_som_elem gpr_som_element;


struct gpr_som_struct {
	int dimension;
	int no_of_sensors;
    struct gpr_som_elem ** weight;
    float * state;
	float * sensor_scale;
};
typedef struct gpr_som_struct gpr_som;

void gpr_som_init(int dimension, int no_of_sensors, gpr_som * som);
void gpr_som_free(gpr_som * som);
int gpr_som_init_sensor(gpr_som * som,
						int sensor_index,
						float min_value, float max_value,
						unsigned int * random_seed);
int gpr_som_init_sensor_from_data(gpr_som * som,
								  int sensor_index,
								  int training_data_field_index,
								  float * training_data,
								  int training_data_fields_per_example,
								  int no_of_training_examples,
								  unsigned int * random_seed);
void gpr_som_learn_from_data(gpr_som * som,
							 int * data_field_index,
							 float * training_data,
							 int fields_per_sample,
							 int no_of_samples,
							 int learning_itterations,
							 int inhibit_radius, int excite_radius,
							 float learning_rate,
							 unsigned int * random_seed,
							 int show_progress);
void gpr_som_outputs_from_data(gpr_som * som,
							   int * data_field_index,
							   float * data,
							   int fields_per_sample,
							   int no_of_samples,
							   float * result);
void gpr_som_learn(gpr_som * som,
				   float * sensors,
				   int inhibit_radius,
				   int excite_radius,
				   float learning_rate);
int gpr_som_update(float * sensors,
				   gpr_som * som,
				   float * x, float * y);
void gpr_som_run(float * sensors,
				 gpr_som * som,
				 float * x, float * y);
void gpr_som_save(gpr_som * som,
				  FILE * fp);
void gpr_som_load(gpr_som * som,
				  FILE * fp);

#endif

/*
 libgpr - a library for genetic programming
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

#include "tests_som.h"

static void test_gpr_som_init()
{
	int dimension = 128;
	int no_of_sensors = 5;
	gpr_som som;

	printf("test_gpr_som_init...");

	gpr_som_init(dimension, no_of_sensors, &som);
	assert(som.dimension == dimension);
	gpr_som_free(&som);

	printf("Ok\n");
}

static void test_gpr_som_init_sensor()
{
	int dimension = 128;
	int no_of_sensors = 5;
	int i,s;
	gpr_som som;
	float min_value = -10;
	float max_value = 10;
	unsigned int random_seed = (int)time(NULL);
	gpr_som_element * elem;

	printf("test_gpr_som_init_sensor...");

	gpr_som_init(dimension, no_of_sensors, &som);
	assert(som.dimension == dimension);

	for (s = 0; s < no_of_sensors; s++) {
		assert(gpr_som_init_sensor(&som, s,
								   min_value, max_value,
								   &random_seed)==0);
	}

	for (i = 0; i < dimension*dimension; i++) {
		elem = (gpr_som_element*)&som.weight[i];
		for (s = 0; s < no_of_sensors; s++) {
			if (elem->vect[s] < min_value) {
				printf("\n%.3f\n",elem->vect[s]);
			}
			assert(elem->vect[s] >= min_value);
			if (elem->vect[s] > max_value) {
				printf("\n%.3f\n",elem->vect[s]);
			}
			assert(elem->vect[s] <= max_value);
		}
	}

	gpr_som_free(&som);

	printf("Ok\n");
}

static void test_gpr_som_init_sensor_from_data()
{
	int dimension = 128;
	int no_of_sensors = 5;
	int i,s;
	gpr_som som;
	float min_value;
	float max_value;
	unsigned int random_seed = (int)time(NULL);
	gpr_som_element * elem;
	int no_of_training_examples=10;
	int training_data_fields_per_example=0;
	float * training_data;

	printf("test_gpr_som_init_sensor_from_data...");

	training_data_fields_per_example = no_of_sensors;

	/* create the training data */
	training_data = (float*)malloc(no_of_training_examples*
								   training_data_fields_per_example*
								   sizeof(float*));
	for (i = 0; i < no_of_training_examples; i++) {
		for (s = 0; s < no_of_sensors; s++) {
			training_data[i*training_data_fields_per_example + s] = i+1+(s*2);
		}
	}

	gpr_som_init(dimension, no_of_sensors, &som);
	assert(som.dimension == dimension);

	for (s = 0; s < no_of_sensors; s++) {
		assert(gpr_som_init_sensor_from_data(&som, s,
											 s, training_data,
											 training_data_fields_per_example,
											 no_of_training_examples,
											 &random_seed)==0);
	}

	for (i = 0; i < dimension*dimension; i++) {
		elem = (gpr_som_element*)&som.weight[i];
		for (s = 0; s < no_of_sensors; s++) {
			min_value = 1+(s*2);
			max_value = no_of_training_examples+1+(s*2);
			if (elem->vect[s] < min_value) {
				printf("\n%.3f\n",elem->vect[s]);
			}
			assert(elem->vect[s] >= min_value);
			if (elem->vect[s] > max_value) {
				printf("\n%.3f\n",elem->vect[s]);
			}
			assert(elem->vect[s] <= max_value);
		}
	}

	gpr_som_free(&som);

	/* free the training data */
	free(training_data);		

	printf("Ok\n");
}

static void test_gpr_som_update()
{
	int dimension = 128;
	int no_of_sensors = 5;
	float sensors[5];
	int i,s,winner;
	gpr_som som;
	float min_value = -10;
	float max_value = 10;
	float x=0,y=0;
	unsigned int random_seed = (int)time(NULL);
	gpr_som_element * elem;

	printf("test_gpr_som_update...");

	gpr_som_init(dimension, no_of_sensors, &som);
	assert(som.dimension == dimension);

	for (i = 0; i < no_of_sensors; i++) {
		assert(gpr_som_init_sensor(&som, i,
								   min_value, max_value,
								   &random_seed)==0);
	}

	for (i = 0; i < dimension*dimension; i++) {
		elem = (gpr_som_element*)&som.weight[i];
		for (s = 0; s < no_of_sensors; s++) {
			assert(elem->vect[s]>=min_value);
			assert(elem->vect[s]<=max_value);
		}
	}

	for (i = 0; i < no_of_sensors; i++) {
		sensors[i] = (rand_num(&random_seed)%10000)/10000.0f;
	}
	
	winner = gpr_som_update((float*)&sensors,
							&som, &x, &y);
	assert(winner>-1);
	assert(x >= 0);
	assert(x <= 1.0f);
	assert(y >= 0);
	assert(y <= 1.0f);
	gpr_som_free(&som);

	printf("Ok\n");
}

static void test_gpr_som_learn()
{
	int dimension = 128;
	int no_of_sensors = 5;
	float sensors[5];
	int i,s,winner,itt;
	gpr_som som;
	float min_value = -10;
	float max_value = 10;
	float x=0,y=0;
	unsigned int random_seed = (int)time(NULL);
	gpr_som_element * elem;
	int inhibit_radius = dimension*10/100;
	int excite_radius = dimension*5/100;
	float learning_rate = 0.2f;

	printf("test_gpr_som_learn...");

	gpr_som_init(dimension, no_of_sensors, &som);
	assert(som.dimension == dimension);

	for (i = 0; i < no_of_sensors; i++) {
		assert(gpr_som_init_sensor(&som, i,
								   min_value, max_value,
								   &random_seed)==0);
	}

	for (i = 0; i < dimension*dimension; i++) {
		elem = (gpr_som_element*)&som.weight[i];
		for (s = 0; s < no_of_sensors; s++) {
			assert(elem->vect[s]>=min_value);
			assert(elem->vect[s]<=max_value);
		}
	}

	for (itt = 0; itt < 500; itt++) {
		for (i = 0; i < no_of_sensors; i++) {
			sensors[i] = (rand_num(&random_seed)%10000)/10000.0f;
		}
	
		winner = gpr_som_update((float*)&sensors,
								&som, &x, &y);
		assert(winner>-1);
		assert(x >= 0);
		assert(x <= 1.0f);
		assert(y >= 0);
		assert(y <= 1.0f);

		gpr_som_learn(&som,
					  (float*)&sensors,
					  inhibit_radius,
					  excite_radius,
					  learning_rate);
	}

	gpr_som_free(&som);

	printf("Ok\n");
}

static void test_gpr_som_save_load()
{
	int dimension = 128;
	int no_of_sensors = 5;
	float sensors[5];
	int i,s,winner,itt;
	gpr_som som1, som2;
	float min_value = -10;
	float max_value = 10;
	float x=0,y=0;
	unsigned int random_seed = (int)time(NULL);
	gpr_som_element * elem, *elem1, *elem2;
	int inhibit_radius = dimension*10/100;
	int excite_radius = dimension*5/100;
	float learning_rate = 0.2f;
	char filename[256];
	FILE * fp;

	printf("test_gpr_som_save_load...");

	gpr_som_init(dimension, no_of_sensors, &som1);
	assert(som1.dimension == dimension);
	assert(som1.no_of_sensors == no_of_sensors);

	for (i = 0; i < no_of_sensors; i++) {
		assert(gpr_som_init_sensor(&som1, i,
								   min_value, max_value,
								   &random_seed)==0);
	}

	for (i = 0; i < dimension*dimension; i++) {
		elem = (gpr_som_element*)&som1.weight[i];
		for (s = 0; s < no_of_sensors; s++) {
			assert(elem->vect[s]>=min_value);
			assert(elem->vect[s]<=max_value);
		}
	}

	for (itt = 0; itt < 500; itt++) {
		for (i = 0; i < no_of_sensors; i++) {
			sensors[i] = (rand_num(&random_seed)%10000)/10000.0f;
		}
	
		winner = gpr_som_update((float*)&sensors,
								&som1, &x, &y);
		assert(winner>-1);
		assert(x >= 0);
		assert(x <= 1.0f);
		assert(y >= 0);
		assert(y <= 1.0f);

		gpr_som_learn(&som1,
					  (float*)&sensors,
					  inhibit_radius,
					  excite_radius,
					  learning_rate);
	}

	/* save the result */
	sprintf(filename,"%slibgpr_som.txt",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"w");
	assert(fp);
	gpr_som_save(&som1, fp);
	fclose(fp);

	/* load */
	fp = fopen(filename,"r");
	assert(fp);
	gpr_som_load(&som2,fp);
	fclose(fp);

	/* compare */
	assert(som1.dimension == som2.dimension);
	assert(som1.no_of_sensors == som2.no_of_sensors);
	for (i = 0; i < som1.dimension*som1.dimension; i++) {
		elem1 = (gpr_som_element*)&som1.weight[i];
		elem2 = (gpr_som_element*)&som2.weight[i];
		for (s = 0; s < som1.no_of_sensors; s++) {
			assert(fabs(elem1->vect[s] - elem2->vect[s])<0.01f);
		}
	}

	gpr_som_free(&som1);
	gpr_som_free(&som2);

	printf("Ok\n");
}

int run_tests_som()
{
	printf("\nRunning SOM tests\n");

	test_gpr_som_init();
	test_gpr_som_init_sensor();
	test_gpr_som_init_sensor_from_data();
	test_gpr_som_update();
	test_gpr_som_learn();
	test_gpr_som_save_load();

	printf("All SOM tests completed\n");
	return 1;
}

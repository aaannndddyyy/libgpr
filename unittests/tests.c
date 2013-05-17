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

#include "tests.h"

static void test_gpr_mutate_value()
{
	unsigned int random_seed = 123;
	int i;
	float result, value = 100;

	printf("test_gpr_mutate_value...");

	for (i = 0; i < 1000; i++) {
		result = gpr_mutate_value(value, 10, &random_seed);
		assert(result >= 90);
		assert(result <= 110.01f);
	}

	printf("Ok\n");
}

static void test_rand_num()
{
	unsigned int random_seed=0;
	int i,j,t,v,min=0,max=0;
	int test1[50];
	int test2[50];
	int test3[50];
	int * result;
	int same=0,repeats=0;

	printf("test_rand_num...");

	/* run three sequences with different seeds */
	for (t = 0; t < 3; t++) {
		switch(t) {
		case 0: { random_seed = 123; result = (int*)test1; break; }
		case 1: { random_seed = 555; result = (int*)test2; break; }
		case 2: { random_seed = 8323; result = (int*)test3; break; }
		}
		for (i = 0; i < 50; i++) {
			result[i] = rand_num(&random_seed);
		}
	}

	for (i = 0; i < 50; i++) {
		/* check that the sequences are different */
		if ((test1[i]==test2[i]) ||
			(test1[i]==test3[i]) ||
			(test2[i]==test3[i])) {
			same++;
		}

		/* check the number of repeats within each sequence */
		for (j = 0; j < 50; j++) {
			if (i!=j) {
				if ((test1[i]==test1[j]) ||
					(test2[i]==test2[j]) ||
					(test3[i]==test3[j])) {
					repeats++;
				}
			}
		}
	}		
	assert(same < 2);
	assert(repeats < 2);

	/* check that the range is not too restricted */
	for (i = 0; i < 10000; i ++) {
		v = rand_num(&random_seed);
		if ((i==0) ||
			((i>0) && (v<min))) {
			min = v;
		}
		if ((i==0) ||
			((i>0) && (v>max))) {
			max = v;
		}
	}
	assert(max > min);
	assert(min >= 0);
	assert(max - min > 60000);

	printf("Ok\n");
}

static  void test_gpr_random_value()
{
	int i, empty=0;
	float testvalues[10];
	float min_value  = -10;
	float max_value = 10;
	float v;
	unsigned int random_seed = 123;

	printf("test_gpr_random_value...");

	for (i = 0; i < 10; i++) {
		testvalues[i]=0;
	}

	/* check that values are within range */
	for (i = 0; i < 10000; i++) {
		v = gpr_random_value(min_value,max_value,&random_seed);
		assert(v>=min_value);
		assert(v<=max_value);
		testvalues[(int)((v-min_value)*9/(max_value-min_value))]=1;
	}

	/* check the distribution across the range */
	for (i = 0; i < 10; i++) {
		if (testvalues[i]==0) empty++;
	}
	
	assert(empty<2);

	printf("Ok\n");
}

static void test_gpr_init()
{
	gpr_function f;
	int ctr=0;

	printf("test_gpr_init...");	

	gpr_init(&f);
	assert(f.function_type==GPR_FUNCTION_VALUE);
	assert(f.value==0);
	assert(f.argc==GPR_DEFAULT_ARGUMENTS);

	gpr_nodes(&f,&ctr);
	assert(ctr==1);
	
	gpr_free(&f);

	printf("Ok\n");
}

static void test_gpr_random()
{
	gpr_function f;
	int i, ctr=0;
	int min_depth=2;
	int max_depth=10,max_depth2,depth=0;
	int validation_result=GPR_VALIDATE_OK;
	float branching_prob=0.8f;
	unsigned int random_seed = 123;
	int instruction_set[64], no_of_instructions=0;
	int integers_only=0;
	int ADFs=0;

	printf("test_gpr_random...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* by running repeatedly we can see if there are any issues with
	   different random permutations */
	for (i = 0; i < GPR_MAX_TESTS; i++) {

		/* create a random tree */
		depth = 0;
		ctr = 0;
		gpr_random(&f,depth,min_depth,max_depth,branching_prob,
				   100,200,integers_only,&random_seed,
				   (int*)instruction_set, no_of_instructions);
		gpr_nodes(&f,&ctr);
		assert(ctr>1);

		/* check the maximum depth wasn't exceeded */
		depth=0;
		max_depth2=0;
		gpr_max_depth(&f, depth, &max_depth2);
		assert(max_depth2<=max_depth);

		/* validate the result */
		depth=0;
		gpr_validate(&f, depth,
					 min_depth, max_depth,
					 ADFs, 0, 0, &validation_result);
		switch(validation_result) {
		case GPR_VALIDATE_TREE_NOT_TERMINATED: {
			printf("\nTree not terminated\n");
			break;
		}
		case GPR_VALIDATE_TREE_TOO_DEEP: {
			printf("\nTree too deep %d %d\n", max_depth, max_depth2);
			break;
		}
		case GPR_VALIDATE_TERMINATOR_AT_ROOT: {
			printf("\nTerminator at root\n");
			break;
		}
		case GPR_VALIDATE_ADF_PATTERN: {
			printf("\nTree does not follow ADF pattern\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_VALUE: {
			printf("\nADF definition contains value\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_ARG: {
			printf("\nADF main program contains argument\n");
			break;
		}
		case GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE: {
			printf("\nADF argument out of range\n");
			break;
		}
		default: {
			if (validation_result!=GPR_VALIDATE_OK) {
				printf("\nUnknown validation return value %d\n",
					   validation_result);
			}
		}
		}
		assert(validation_result==GPR_VALIDATE_OK);

		gpr_free(&f);

		ctr=0;
		gpr_nodes(&f,&ctr);
		assert(ctr==0);
	}

	printf("Ok\n");
}

static void test_gpr_prune()
{
	gpr_function f;
	int i, ctr=0;
	int min_depth=2;
	int max_depth,max_depth2,depth;
	float branching_prob=0.8f;
	float min_value = 100;
	float max_value = 200;
	unsigned int random_seed = 123;
	int integers_only=0;
	int instruction_set[64], no_of_instructions=0;	

	printf("test_gpr_prune...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* by running repeatedly we can see if there are any issues with
	   different random permutations */
	for (i = 0;i < GPR_MAX_TESTS; i++) {

		/* create the tree */
		max_depth=10;
		depth = 0;
		ctr = 0;
		gpr_random(&f,depth,min_depth,max_depth,branching_prob,
				   min_value,max_value,integers_only,&random_seed,
				   (int*)instruction_set, no_of_instructions);
		gpr_nodes(&f,&ctr);
		assert(ctr>1);

		/* check the depth */
		depth=0;
		max_depth2=0;
		gpr_max_depth(&f,depth,&max_depth2);
		assert(max_depth2>=min_depth);

		/* prune the tree to a lesser depth */
		max_depth=5;
		depth=0;
		gpr_prune(&f, depth, max_depth,
				  min_value, max_value,
				  &random_seed);

		/* check the depth */
		depth=0;
		max_depth2=0;
		gpr_max_depth(&f,depth,&max_depth2);
		assert(max_depth2<=max_depth);

		gpr_free(&f);

		ctr=0;
		gpr_nodes(&f,&ctr);
		assert(ctr==0);
	}

	printf("Ok\n");
}

static void test_gpr_copy()
{
	gpr_function f,f2;
	int i, ctr=0,ctr2=0;
	int min_depth=2;
	int max_depth=10,depth=0;
	float branching_prob=0.8f;
	unsigned int random_seed = 123;
	int integers_only=0;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gpr_copy...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* by running repeatedly we can see if there are any issues with
	   different random permutations */
	for (i = 0; i < GPR_MAX_TESTS; i++) {

		depth = 0;
		ctr = 0;
		gpr_random(&f,depth,min_depth,max_depth,branching_prob,
				   100, 200, integers_only, &random_seed,
				   (int*)instruction_set, no_of_instructions);
		gpr_nodes(&f,&ctr);
		assert(ctr>1);

		/* make a copy */
		ctr2 = 0;
		gpr_copy(&f,&f2);
		gpr_nodes(&f2,&ctr2);
		assert(ctr==ctr2);

		gpr_free(&f);
		gpr_free(&f2);

		ctr = 0;
		gpr_nodes(&f,&ctr);
		assert(ctr==0);

		ctr = 0;
		gpr_nodes(&f2,&ctr);
		assert(ctr==0);

	}

	printf("Ok\n");
}

static void test_gpr_mutate()
{
	gpr_function f,f2;
	int i, ctr=0,ctr2=0,diffs=0;
	int min_depth=2;
	int max_depth=10,max_depth1,max_depth2,depth=0;
	float branching_prob=0.8f;
	float mutation_prob=0.5f;
	int validation_result=GPR_VALIDATE_OK;
	unsigned int random_seed = 123;
	float min_value = -5;
	float max_value = 5;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gpr_mutate...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* by running repeatedly we can see if there are any issues with
	   different random permutations */
	for (i = 0; i < GPR_MAX_TESTS; i++) {

		depth = 0;
		ctr = 0;
		gpr_random(&f,depth,min_depth,max_depth,branching_prob,
				   min_value, max_value, integers_only, &random_seed,
				   (int*)instruction_set, no_of_instructions);
		gpr_nodes(&f,&ctr);
		assert(ctr>1);

		/* make a copy */
		ctr2 = 0;
		gpr_copy(&f,&f2);
		gpr_nodes(&f2,&ctr2);
		assert(ctr==ctr2);

		/* tree depth before copying */
		depth=0;
		max_depth1=0;
		gpr_max_depth(&f, depth, &max_depth1);
		assert(max_depth1 <= max_depth);

		/* tree depth after copying */
		depth=0;
		max_depth2=0;
		gpr_max_depth(&f2, depth, &max_depth2);
		assert(max_depth2 <= max_depth);

		/* same tree depth before and after copying */
		assert(max_depth1==max_depth2);

		/* mutate the copy */
		depth = 0;
		gpr_mutate(&f2, depth, max_depth, mutation_prob,
				   min_value, max_value,
				   integers_only, ADFs, 0, 0,
				   &random_seed,
				   (int*)instruction_set, no_of_instructions);

		ctr2=0;
		gpr_nodes(&f2,&ctr2);
		assert(ctr2>0);

		/* difference in number of nodes
		   between the original and mutated version */
		diffs += abs(ctr2-ctr);

		depth=0;
		gpr_validate(&f2,depth,0,max_depth,ADFs,0,0,&validation_result);

		switch(validation_result) {
		case GPR_VALIDATE_TREE_NOT_TERMINATED: {
			printf("\nTree not terminated\n");
			break;
		}
		case GPR_VALIDATE_TREE_TOO_DEEP: {
			printf("\nTree too deep %d %d\n", max_depth, max_depth2);
			break;
		}
		case GPR_VALIDATE_TERMINATOR_AT_ROOT: {
			printf("\nTerminator at root\n");
			break;
		}
		case GPR_VALIDATE_ADF_PATTERN: {
			printf("\nTree does not follow ADF pattern\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_VALUE: {
			printf("\nADF definition contains value\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_ARG: {
			printf("\nADF main program contains argument\n");
			break;
		}
		case GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE: {
			printf("\nADF argument out of range\n");
			break;
		}
		default: {
			if (validation_result!=GPR_VALIDATE_OK) {
				printf("\nUnknown validation return value %d\n",
					   validation_result);
			}
		}
		}

		assert (validation_result==GPR_VALIDATE_OK);

		gpr_free(&f);
		gpr_free(&f2);

		ctr=0;
		gpr_nodes(&f,&ctr);
		assert(ctr==0);
	}

	/* There should be some differences between the original
	   and the mutated version */
	assert(diffs > 0);

	printf("Ok\n");
}

static void test_gpr_crossover()
{
	gpr_function parent1,parent2,child;
	int i, ctr=0,ctr2,success;
	int min_depth=2;
	int max_depth=10,max_depth2,depth=0;
	float branching_prob=0.8f;
	float min_value = 100;
	float max_value = 200;
	int validation_result = 0;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs=0;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gpr_crossover...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* by running repeatedly we can see if there are any issues with
	   different random permutations */
	for (i = 0; i < GPR_MAX_TESTS; i++) {

		/* create parent1 */
		depth=0;
		gpr_random(&parent1,depth,min_depth,max_depth,branching_prob,
				   min_value, max_value, integers_only, &random_seed,
				   (int*)instruction_set, no_of_instructions);
		ctr=0;
		gpr_nodes(&parent1,&ctr);
		assert(ctr>1);

		depth=0;
		gpr_validate(&parent1,depth,min_depth,max_depth,ADFs,
					 0,0,&validation_result);
		assert (validation_result==GPR_VALIDATE_OK);

		/* create parent2 */
		depth=0;
		gpr_random(&parent2,depth,min_depth,max_depth,branching_prob,
				   min_value, max_value, integers_only, &random_seed,
				   (int*)instruction_set, no_of_instructions);
		ctr2=0;
		gpr_nodes(&parent2,&ctr2);
		assert(ctr2>1);

		depth=0;
		gpr_validate(&parent2,depth,min_depth,max_depth,ADFs,
					 0,0,&validation_result);
		assert (validation_result==GPR_VALIDATE_OK);

		/* perform crossover */
		success = gpr_crossover(&parent1, &parent2, &child, max_depth,
								min_value, max_value,
								&random_seed);
		assert(success==1);

		ctr=0;
		gpr_nodes(&child,&ctr);
		assert(ctr>0);

		depth=0;
		gpr_validate(&child,depth,0,max_depth,ADFs,
					 0,0,&validation_result);

		switch(validation_result) {
		case GPR_VALIDATE_TREE_NOT_TERMINATED: {
			printf("\nTree not terminated\n");
			break;
		}
		case GPR_VALIDATE_TREE_TOO_DEEP: {
			depth=0;
			max_depth2=0;
			gpr_max_depth(&child, depth, &max_depth2);
			printf("\nTree too deep %d %d\n", max_depth, max_depth2);
			break;
		}
		case GPR_VALIDATE_TERMINATOR_AT_ROOT: {
			printf("\nTerminator at root\n");
			break;
		}
		case GPR_VALIDATE_ADF_PATTERN: {
			printf("\nTree does not follow ADF pattern\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_VALUE: {
			printf("\nADF definition contains value\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_ARG: {
			printf("\nADF main program contains argument\n");
			break;
		}
		case GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE: {
			printf("\nADF argument out of range\n");
			break;
		}
		default: {
			if (validation_result!=GPR_VALIDATE_OK) {
				printf("\nUnknown validation return value %d\n",
					   validation_result);
			}
		}
		}

		assert (validation_result==GPR_VALIDATE_OK);

		/* cleanup */
		gpr_free(&parent1);
		gpr_free(&parent2);
		gpr_free(&child);

	}

	printf("Ok\n");
}

static void test_gpr_mate()
{
	gpr_function parent1,parent2,child;
	gpr_state child_state;
	int i, ctr=0,ctr2;
	int min_depth=2;
	int max_depth=10,max_depth2,depth=0;
	float branching_prob=0.8f;
	float min_value = 100;
	float max_value = 200;
	float mutation_prob = 0.5f;
	float pure_mutant_prob=0.3f;
	int validation_result = 0;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 0;
	int registers=4, sensors=3, actuators=2;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_mate...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	child.argv[0]=0;
	child.argv[1]=0;

	/* by running repeatedly we can see if there are any issues with
	   different random permutations */
	for (i = 0; i < GPR_MAX_TESTS; i++) {

		/* initialise some state */
		gpr_init_state(&child_state,
					   registers, sensors, actuators,
					   data_size, data_fields,
					   &random_seed);

		/* create parent1 */
		depth=0;
		gpr_random(&parent1,depth,min_depth,max_depth,branching_prob,
				   min_value, max_value, integers_only, &random_seed,
				   (int*)instruction_set,no_of_instructions);
		ctr=0;
		gpr_nodes(&parent1,&ctr);
		assert(ctr>1);

		depth=0;
		gpr_validate(&parent1,depth,min_depth,max_depth,ADFs,
					 0,0,&validation_result);
		assert (validation_result==GPR_VALIDATE_OK);

		/* create parent2 */
		depth=0;
		gpr_random(&parent2,depth,min_depth,max_depth,branching_prob,
				   min_value, max_value, integers_only, &random_seed,
				   (int*)instruction_set,no_of_instructions);
		ctr2=0;
		gpr_nodes(&parent2,&ctr2);
		assert(ctr2>1);

		depth=0;
		gpr_validate(&parent2,depth,min_depth,max_depth,ADFs,
					 0,0,&validation_result);
		assert (validation_result==GPR_VALIDATE_OK);

		gpr_mate(&parent1, &parent2, max_depth,
				 min_value, max_value,
				 mutation_prob, pure_mutant_prob,
				 integers_only, ADFs, &random_seed,
				 (int*)instruction_set,no_of_instructions,
				 &child,&child_state);

		ctr=0;
		gpr_nodes(&child,&ctr);
		assert(ctr>0);

		depth=0;
		gpr_validate(&child,depth,0,max_depth,ADFs,
					 0,0,&validation_result);

		switch(validation_result) {
		case GPR_VALIDATE_TREE_NOT_TERMINATED: {
			printf("\nTree not terminated\n");
			break;
		}
		case GPR_VALIDATE_TREE_TOO_DEEP: {
			depth=0;
			max_depth2=0;
			gpr_max_depth(&child, depth, &max_depth2);
			printf("\nTree too deep %d %d\n", max_depth, max_depth2);
			break;
		}
		case GPR_VALIDATE_TERMINATOR_AT_ROOT: {
			printf("\nTerminator at root\n");
			break;
		}
		case GPR_VALIDATE_ADF_PATTERN: {
			printf("\nTree does not follow ADF pattern\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_VALUE: {
			printf("\nADF definition contains value\n");
			break;
		}
		case GPR_VALIDATE_ADF_CONTAINS_ARG: {
			printf("\nADF main program contains argument\n");
			break;
		}
		case GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE: {
			printf("\nADF argument out of range\n");
			break;
		}
		default: {
			if (validation_result!=GPR_VALIDATE_OK) {
				printf("\nUnknown validation return value %d\n",
					   validation_result);
			}
		}
		}

		assert (validation_result==GPR_VALIDATE_OK);

		/* cleanup */
		gpr_free(&parent1);
		gpr_free(&parent2);
		gpr_free(&child);
		gpr_free_state(&child_state);
	}

	printf("Ok\n");
}

static void test_gpr_run()
{
	gpr_function f;
	int i, j, ctr=0;
	int min_depth=2;
	int max_depth=10,depth=0;
	float branching_prob=0.8f;
	gpr_state state;
	int no_of_registers = 4;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_run...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create machine state */
	gpr_init_state(&state, no_of_registers, 0, 0,
				   data_size, data_fields,
				   &random_seed);

	/* by running repeatedly we can see if there are any issues with
	   different random permutations */
	for (i = 0; i < 100; i++) {

		/* create a random tree */
		depth = 0;
		ctr = 0;
		gpr_random(&f,depth,min_depth,max_depth,branching_prob,
				   100,200, integers_only, &random_seed,
				   (int*)instruction_set,no_of_instructions);
		gpr_nodes(&f,&ctr);
		assert(ctr>1);

		/* run the program for some steps */
		for (j = 0; j < 10; j++) {
			gpr_run(&f, &state, 0);
		}

		/* cleanup */
		gpr_free(&f);

		ctr=0;
		gpr_nodes(&f,&ctr);
		assert(ctr==0);
	}

	/* free machine state */
	assert(state.no_of_registers==no_of_registers);

	gpr_free_state(&state);
	
	printf("Ok\n");
}

static void test_gpr_sort()
{
	int population_size = 1000;
	int i, max_depth=10;
	gpr_population population;
	float min_value = 100;
	float max_value = 200;
	float f;
	gpr_function * individual, * prev_individual;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_sort...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gpr_init_population(&population, population_size, 4, 0, 0,
						max_depth, min_value, max_value,
						integers_only, ADFs,
						data_size, data_fields,
						&random_seed,
						(int*)instruction_set, no_of_instructions);

	/* set random fitness values */
	for (i = 0; i < population_size; i++) {
		f = (float)(rand_num(&random_seed)%10000);
		(&population)->fitness[i] = f;
		individual = &(population.individual[i]);
		individual->value = f;
	}

	/* sort by fitness */
	gpr_sort(&population);

	/* check the sort order */	
	for (i = 1; i < population_size; i++) {
		individual = &(population.individual[i]);
		prev_individual = &(population.individual[i-1]);
		assert(population.fitness[i-1] >= population.fitness[i]);
		assert(prev_individual->value >= individual->value);
		assert(population.fitness[i-1] == prev_individual->value);
		assert(population.fitness[i] == individual->value);
	}	

	/* free memory */
	gpr_free_population(&population);

	printf("Ok\n");
}

static void test_gpr_sort_system()
{
	int population_per_island = 256;
	int i, j, max_tree_depth = 10;
	int registers = 4;
	int sensors = 4;
	int actuators = 2;
	gpr_system system;
	float min_value = 100;
	float max_value = 200;
	float f;
	int islands = 4;
	gpr_population * population;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_sort_system...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a system */
	gpr_init_system(&system,
					islands,
					population_per_island,
					registers,
					sensors,
					actuators,
					max_tree_depth,
					min_value, max_value,
					integers_only,
					ADFs,
					data_size, data_fields,
					&random_seed,
					instruction_set, no_of_instructions);

	/* set random fitness values */
	for (i = 0; i < islands; i++) {
		f = (float)(rand_num(&random_seed)%10000);
		(&system)->fitness[i] = f;
		population = &system.island[i];
		for (j = 0; j < population->size; j++) {
			population->fitness[j] =
				(float)(rand_num(&random_seed)%10000);
		}
	}

	/* sort by fitness */
	gpr_sort_system(&system);

	/* check the sort order */	
	for (i = 1; i < islands; i++) {
		assert(system.fitness[i-1] >= system.fitness[i]);
	}	

	/* free memory */
	gpr_free_system(&system);

	printf("Ok\n");
}

/* A test evaluation function.
   This tests how close the output is to the
   equation y = 3x^2 + 2x - 5 */
float test_evaluate_program(int time_steps,
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
		/* observe the output */
		result = gpr_get_actuator(state,0);
		/* the target value */
		reference = (3*x*x) + (2*x) - 5;
		/* how close is the output to the target equation? */
		fitness += 100-fabs(result - reference);
	}
	return fitness;
}

static void test_gpr_ADF_population()
{
	int i, gen, population_size = 512;
	int max_depth=5;
	gpr_population population;
	float min_value = -5;
	float max_value = 5;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 1;
	int instruction_set[64], no_of_instructions=0;
	int validation_result=0;
	char filename[256];
	FILE * fp;
	float elitism = 0.3f;
	float mutation_prob=0.2f;
	float pure_mutant_prob=0.1f;
	int time_steps=10;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_ADF_population...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gpr_init_population(&population, population_size, 4, 1, 1,
						max_depth, min_value, max_value,
						integers_only, ADFs,
						data_size, data_fields,
						&random_seed,
						(int*)instruction_set,no_of_instructions);

	for (gen = 0; gen < 10; gen++) {

		for (i = 0; i < population.size; i++) {

			gpr_validate(&population.individual[i], 0,
						 2, max_depth,ADFs, 0, 0,
						 &validation_result);

			switch(validation_result) {
			case GPR_VALIDATE_TREE_NOT_TERMINATED: {
				printf("\nTree not terminated\n");
				break;
			}
			case GPR_VALIDATE_TREE_TOO_DEEP: {
				printf("\nTree too deep %d\n", max_depth);
				break;
			}
			case GPR_VALIDATE_TERMINATOR_AT_ROOT: {
				printf("\nTerminator at root\n");
				break;
			}
			case GPR_VALIDATE_ADF_PATTERN: {
				printf("\nTree does not follow ADF pattern\n");
				break;
			}
			case GPR_VALIDATE_ADF_CONTAINS_VALUE: {
				printf("\nADF definition contains value\n");
				break;
			}
			case GPR_VALIDATE_ADF_CONTAINS_ARG: {
				printf("\nADF main program contains argument\n");
				break;
			}
			case GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE: {
				printf("\nADF argument out of range\n");
				break;
			}
			default: {
				if (validation_result!=GPR_VALIDATE_OK) {
					printf("\nUnknown validation return value %d\n",
						   validation_result);
				}
			}
			}

			if (validation_result != GPR_VALIDATE_OK) {
				sprintf(filename,"%s","adf.dot");
				fp = fopen(filename,"w");
				assert(fp!=0);
				gpr_dot(&population.individual[i],
						fp);
				fclose(fp);
			}

			assert(validation_result==GPR_VALIDATE_OK);
		}

		/* evaluate each individual */
		gpr_evaluate(&population,
					 time_steps,0,
					 (*test_evaluate_program));

		/* produce the next generation */
		gpr_generation(&population,
					   elitism,
					   max_depth,
					   min_value, max_value,
					   mutation_prob, pure_mutant_prob,
					   integers_only,
					   ADFs,
					   (int*)instruction_set, no_of_instructions);

		/*if (gen == 0) start_fitness = gpr_best_fitness(&population);*/
	}

	assert(population.history.index == gen);


	/* save as dot file */
	sprintf(filename,"%stest.dot",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"w");
	assert(fp!=0);
	gpr_dot(&population.individual[1],fp);
	fclose(fp);

	gpr_free_population(&population);

	printf("Ok\n");
}

static void test_gpr_generation()
{
	int population_size = 512;
	int gen, max_depth = 5;
	gpr_population population;
	float min_value = -5;
	float max_value = 5;
	float elitism = 0.2f;
	float mutation_prob=0.5f;
	float pure_mutant_prob=0.2f;
	int time_steps = 10;
	/*float start_fitness = 0;*/
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_generation...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gpr_init_population(&population, population_size, 4, 1, 1,
						max_depth, min_value, max_value,
						integers_only, ADFs,
						data_size, data_fields,
						&random_seed,
						(int*)instruction_set,no_of_instructions);

	for (gen = 0; gen < 20; gen++) {

		/* evaluate each individual */
		gpr_evaluate(&population,
					 time_steps,0,
					 (*test_evaluate_program));

		/* produce the next generation */
		gpr_generation(&population,
					   elitism,
					   max_depth,
					   min_value, max_value,
					   mutation_prob, pure_mutant_prob,
					   integers_only, ADFs,
					   (int*)instruction_set, no_of_instructions);

		/*if (gen == 0) start_fitness = gpr_best_fitness(&population);*/
	}

	assert(population.history.index == gen);

	/* fitness is increasing */
	/*assert(gpr_best_fitness(&population) > start_fitness);*/

	/* free memory */
	gpr_free_population(&population);

	printf("Ok\n");
}

static void test_gpr_environment()
{
	int result, i, n, max_population_size = 64;
	int population_size = 32;
	int max_depth=5;
	gpr_environment population, population2;
	float min_value = -5;
	float max_value = 5;
	float mutation_prob=0.5f;
	float pure_mutant_prob=0.2f;
	unsigned int random_seed = 1234;
	int integers_only = 0;
	int ADFs = 0;
	int sensors = 4, actuators = 2, registers = 3;
	int parent1_index, parent2_index, child_index, victim_index;
	int instruction_set[64], no_of_instructions=0;
	FILE * fp;
	char filename[256];
	int data_size = 8, data_fields = 2;

	printf("test_gpr_environment...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gpr_init_environment(&population,
						 max_population_size,
						 population_size,
						 registers, sensors, actuators,
						 max_depth, min_value, max_value,
						 integers_only, ADFs,
						 data_size, data_fields,
						 &random_seed,
						 (int*)instruction_set,no_of_instructions);

	/* clear the number of matings */
	population.matings = 0;
	for (i = 0; i < 33; i++) {
		/* randomly select parents */
		parent1_index =
			rand_num(&random_seed)%population.population_size;
		parent2_index =
			rand_num(&random_seed)%population.population_size;
		/* mate */
		child_index =
			gpr_mate_environment(&population,
								 parent1_index, parent2_index,
								 max_depth, min_value, max_value,
								 mutation_prob, pure_mutant_prob,
								 integers_only, ADFs,
								 instruction_set, no_of_instructions);
		/* check the child indexes */
		if (i < 32) {
			n = population.matings-1;
			assert(child_index == population.population_size-1);
			assert(population.mating[n*3] == parent1_index);
			assert(population.mating[n*3+1] == parent2_index);
			assert(population.mating[n*3+2] == child_index);			
		}
		else {
			assert(child_index == -1);
		}
	}
	assert(population.population_size == 64);

	/* some individuals die */
	for (i = 0; i < 10; i++) {
		victim_index =
			rand_num(&random_seed)%population.population_size;
		gpr_death(&population, victim_index);
	}
	/* check reduced population size */
	if (population.population_size != 54) {
		printf("\nPopulation size: %d\n",
			   population.population_size);
	}
	assert(population.population_size == 54);

	/* save to file */
	sprintf(filename,"%stemp_gpr_env.dat",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"wb");
	assert(fp);
	gpr_save_environment(&population, fp);
	fclose(fp);

	/* load the file */
	fp = fopen(filename,"rb");
	assert(fp);
	gpr_load_environment(&population2, fp);
	fclose(fp);

	/* compare the populations */
	assert(population.max_population_size ==
		   population2.max_population_size);
	assert(population.population_size ==
		   population2.population_size);
	assert(population2.matings == 0);
	
	for (i = 0; i < population2.population_size; i++) {
		result = 0;
		gpr_functions_are_equal(&population.individual[i],
								&population2.individual[i],
								&result);
		if (result != 0) {
			sprintf(filename,"%s","temp1.dot");
			fp = fopen(filename,"w");
			assert(fp!=0);
			gpr_dot(&population.individual[i], fp);
			fclose(fp);
			sprintf(filename,"%s","temp2.dot");
			fp = fopen(filename,"w");
			assert(fp!=0);
			gpr_dot(&population2.individual[i], fp);
			fclose(fp);

			printf("\nresult %d/%d: %d\n",
				   i, population2.population_size, result);
			printf("Programs saved as temp1.dot and temp2.dot\n");
				   
		}
		assert(result == 0);
	}

	/* free memory */
	gpr_free_environment(&population);
	gpr_free_environment(&population2);

	printf("Ok\n");
}

static void test_gpr_generation_system()
{
	int population_per_island = 256;
	int islands = 4;
	int migration_interval = 10;
	int gen, max_depth=5;
	gpr_system sys;
	float min_value = -5;
	float max_value = 5;
	float elitism = 0.2f;
	float mutation_prob=0.5f;
	float pure_mutant_prob=0.2f;
	int time_steps = 10;
	/*float start_fitness = 0;*/
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 0;
	int instruction_set[64], no_of_instructions=0;
	char source_filename[256],compile_command[256];
	char binary_filename[256],result_filename[256];
	char line[256];
	int i, sensors=1, actuators=3, registers=4;
	int separators, itt;
	FILE * fp;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_generation_system...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	for (itt = 0; itt < 20; itt++) {
		actuators = 1 + itt;

		/* create a population */
		gpr_init_system(&sys, islands, population_per_island,
						registers, sensors, actuators,
						max_depth, min_value, max_value,
						integers_only, ADFs,
						data_size, data_fields,
						&random_seed,
						(int*)instruction_set,no_of_instructions);

		for (gen = 0; gen < 20; gen++) {

			/* evaluate each individual */
			gpr_evaluate_system(&sys,
								time_steps,0,
								(*test_evaluate_program));

			/* produce the next generation */
			gpr_generation_system(&sys,
								  migration_interval,
								  elitism,
								  max_depth,
								  min_value, max_value,
								  mutation_prob, pure_mutant_prob,
								  integers_only,
								  ADFs,
								  (int*)instruction_set,
								  no_of_instructions);

			/*if (gen == 0) start_fitness = gpr_best_fitness_system(&system);*/
		}

		/* fitness is increasing */
		/*assert(gpr_best_fitness_system(&system) > start_fitness);*/

		/* save the best individual as a C program */
		sprintf(source_filename,"%sagent.c",GPR_TEMP_DIRECTORY);
		fp = fopen(source_filename,"w");
		assert(fp);
		gpr_c_program(gpr_best_individual_system(&sys),
					  sensors, actuators, registers,
					  ADFs, fp);
		fclose(fp);

		/* compile the result */
		sprintf(binary_filename,"%s","temp_agent");
		sprintf(compile_command,
				"gcc -Wall -std=c99 -pedantic -o %s %s -lm",
				binary_filename,source_filename);
		assert(system(compile_command)==0);

		sprintf(result_filename,"%sresult.txt",GPR_TEMP_DIRECTORY);
		sprintf(compile_command,
				"./%s %f > %s",
				binary_filename,2.0f,result_filename);
		assert(system(compile_command)==0);
		fp = fopen(result_filename,"r");
		if (!fp) printf("\nNo result was produced from " \
						"compiled C program\n");
		assert(fp);
		separators=0;
		while (!feof(fp)) {
			if (fgets(line , 255 , fp) != NULL ) {
				assert(strlen(line)>0);
				for (i = 0; i < strlen(line); i++) {
					if ((line[i]==' ') ||
						(line[i]==',') ||
						(line[i]==';')) {
						separators++;
					}
				}
			}
		}
		fclose(fp);

		assert(separators==actuators-1);

		/* free memory */
		gpr_free_system(&sys);

	}

	printf("Ok\n");
}

void test_gpr_dot()
{
	gpr_function f;
	int ctr=0;
	int min_depth=2;
	int max_depth=10,depth=0;
	float branching_prob=0.8f;
	FILE * fp;
	char filename[128],str[128];
	unsigned int random_seed = 123;
	int integers_only = 0;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gpr_dot...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a random tree */
	depth = 0;
	ctr = 0;
	gpr_random(&f,depth,min_depth,max_depth,branching_prob,
			   100,200, integers_only, &random_seed,
			   (int*)instruction_set, no_of_instructions);
	gpr_nodes(&f,&ctr);
	assert(ctr>1);

	/* save as dot file */
	sprintf(filename,"%stest.dot",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"w");
	assert(fp!=0);
	gpr_dot(&f,fp);
	fclose(fp);

	gpr_free(&f);

	ctr=0;
	gpr_nodes(&f,&ctr);
	assert(ctr==0);

	/* delete the test file */
	sprintf(str,"rm %s",filename);
	ctr = system(str);

	printf("Ok\n");
}

static void test_gpr_save_load()
{
	gpr_function f1, f2;
	int ctr=0,ctr2=0;
	int min_depth=2;
	int max_depth=10,depth=0,retval;
	float branching_prob=0.8f;
	FILE * fp;
	char filename[128],str[128];
	unsigned int random_seed = 123;
	int integers_only = 0;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gpr_save_load...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a random tree */
	depth = 0;
	ctr = 0;
	gpr_random(&f1,depth,min_depth,max_depth,branching_prob,
			   100,200, integers_only, &random_seed,
			   (int*)instruction_set, no_of_instructions);
	gpr_nodes(&f1,&ctr);
	assert(ctr>1);
	assert(ctr<GPR_MAX_NODES);

	/* save */
	sprintf(filename,"%stestnode.dat",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"w");
	assert(fp);
	gpr_save(&f1, fp);
	fclose(fp);

	/* load */
	fp = fopen(filename,"r");
	assert(fp);
	retval = gpr_load(&f2,fp);
	fclose(fp);

	switch(retval) {
	case GPR_LOAD_MAX_NODES_REACHED: {
		printf("\nMaximum nodes reached\n");
		break;
	}
	case GPR_LOAD_ARG_NOT_FOUND: {
		printf("\nArgument index not found\n");
		break;
	}
	case GPR_LOAD_CHILD_NOT_FOUND: {
		printf("\nChild node number not found\n");
		break;
	}
	case GPR_LOAD_PARENT_NOT_FOUND: {
		printf("\nParent node number not found\n");
		break;
	}
	case GPR_LOAD_VALUE_NOT_FOUND: {
		printf("\nValue not found\n");
		break;
	}
	case GPR_LOAD_FUNCTION_TYPE_NOT_FOUND: {
		printf("\nFunction type not found\n");
		break;
	}
	case GPR_LOAD_LINK_NOT_FOUND: {
		printf("\nLink not found\n");
		break;
	}
	case GPR_LOAD_NODE_NOT_FOUND: {
		printf("\nNode not found\n");
		break;
	}
	}

	/* no load errors */
	assert(retval == GPR_LOAD_OK);

	/* same number of nodes as the original */
	gpr_nodes(&f2,&ctr2);
	assert(ctr==ctr2);

	/* check that the functions are the same */
	retval = 0;
	gpr_functions_are_equal(&f1, &f2, &retval);
	assert(retval==0);

	/* free memory */
	gpr_free(&f1);
	gpr_free(&f2);

	/* delete the test file */
	sprintf(str,"rm %s",filename);
	ctr = system(str);

	printf("Ok\n");
}

static void test_gpr_save_load_population()
{
	int population_size = 1000;
	int i, retval, max_depth=5;
	gpr_population population,population2;
	gpr_state * state;
	float min_value = -5;
	float max_value = 5;
	char filename[128],str[128];
	FILE * fp;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 0;
	int sensors = 1, actuators = 1, registers = 4;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_save_load_population...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* try with and without ADFs */
	for (ADFs = 0; ADFs <= 1; ADFs++) {

		/* create a population */
		gpr_init_population(&population, population_size,
							registers, sensors, actuators,
							max_depth, min_value, max_value,
							integers_only, ADFs,
							data_size, data_fields,
							&random_seed,
							(int*)instruction_set,no_of_instructions);

		/* save to file */
		sprintf(filename,"%stestpopulation.dat",GPR_TEMP_DIRECTORY);
		fp = fopen(filename,"w");
		assert(fp!=0);
		gpr_save_population(&population,fp);
		fclose(fp);

		/* load from file */
		fp = fopen(filename,"r");
		assert(fp!=0);
		retval = gpr_load_population(&population2, fp);
		fclose(fp);

		switch(retval) {
		case GPR_LOAD_POPULATION_SIZE: {
			printf("Population size does not match the " \
				   "number of individuals loaded\n");
			break;
		}
		case GPR_LOAD_MAX_NODES_REACHED: {
			printf("\nMaximum nodes reached\n");
			break;
		}
		case GPR_LOAD_ARG_NOT_FOUND: {
			printf("\nArgument index not found\n");
			break;
		}
		case GPR_LOAD_CHILD_NOT_FOUND: {
			printf("\nChild node number not found\n");
			break;
		}
		case GPR_LOAD_PARENT_NOT_FOUND: {
			printf("\nParent node number not found\n");
			break;
		}
		case GPR_LOAD_VALUE_NOT_FOUND: {
			printf("\nValue not found\n");
			break;
		}
		case GPR_LOAD_FUNCTION_TYPE_NOT_FOUND: {
			printf("\nFunction type not found\n");
			break;
		}
		case GPR_LOAD_LINK_NOT_FOUND: {
			printf("\nLink not found\n");
			break;
		}
		case GPR_LOAD_NODE_NOT_FOUND: {
			printf("\nNode not found\n");
			break;
		}
		}

		assert(retval == GPR_LOAD_OK);
		assert(population.size==population2.size);
		assert(population.history.index==population2.history.index);
		assert(population.history.tick==population2.history.tick);
		assert(population.history.interval==
			   population2.history.interval);

		/* check that the functions are the same */
		for (i = 0; i < population.size; i++) {
			retval = 0;
			gpr_functions_are_equal(&population.individual[i],
									&population2.individual[i],
									&retval);
			assert (retval==0);

			/* check whether ADFs have been defined */
			state = &population.state[i];
			if (ADFs==0) {
				assert(state->ADF[0]==0);
			}
			else {
				assert(state->ADF[0]!=0);
			}

			state = &population2.state[i];
			if (ADFs==0) {
				assert(state->ADF[0]==0);
			}
			else {
				assert(state->ADF[0]!=0);
			}
		}

		/* free memory */
		gpr_free_population(&population);
		gpr_free_population(&population2);

		/* delete the test file */
		sprintf(str,"rm %s",filename);
		retval = system(str);

	}

	printf("Ok\n");
}

static void test_gpr_init_state()
{
	int population_size = 1000;
	int p, i, max_depth=5;
	gpr_population population;
	float min_value = -5;
	float max_value = 5;
	unsigned int random_seed = 1233;
	int instruction_set[64], no_of_instructions=0;
	gpr_state * state;
	int sensors = 10;
	int actuators = 5;
	int registers = 4;
	int integers_only = 0;
	int ADFs = 0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_init_state...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gpr_init_population(&population, population_size,
						registers, sensors, actuators,
						max_depth, min_value, max_value,
						integers_only, ADFs,
						data_size, data_fields,
						&random_seed,
						(int*)instruction_set,no_of_instructions);

	for (p = 0; p < population_size; p++) {
		state = &population.state[p];

		assert(state->no_of_sensors == sensors);
		for (i = 0; i < sensors; i++) {
			assert(state->sensors[i]==0);
		}

		assert(state->no_of_actuators == actuators);
		for (i = 0; i < actuators; i++) {
			assert(state->actuators[i]==0);
		}

		assert(state->no_of_registers == registers);
		for (i = 0; i < registers; i++) {
			assert(state->registers[i]==0);
		}
		assert(state->random_seed != 0);
	}

	/* free memory */
	gpr_free_population(&population);

	printf("Ok\n");
}

static void test_gpr_save_load_system()
{
	int islands = 4;
	int population_per_island = 256;
	int i, j, k, retval, max_depth=5;
	gpr_system system1,system2;
	gpr_population *population1,*population2;
	gpr_state * state1, * state2;
	float min_value = -5;
	float max_value = 5;
	char filename[128],str[128];
	FILE * fp;
	unsigned int random_seed = 123;
	int integers_only = 0;
	int ADFs = 0;
	int sensors=10, actuators=5, registers=4;
	int no_of_sensor_sources=160;
	int no_of_actuator_destinations=72;
	int instruction_set[64], no_of_instructions=0;
	int data_size = 8, data_fields = 2;

	printf("test_gpr_save_load_system...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a population */
	gpr_init_system(&system1, islands, population_per_island,
					registers, sensors, actuators,
					max_depth, min_value, max_value,
					integers_only, ADFs,
					data_size, data_fields,
					&random_seed,
					(int*)instruction_set,no_of_instructions);

	gpr_init_sensor_sources(&system1,
							sensors,
							no_of_sensor_sources,
							&random_seed);

	gpr_init_actuator_destinations(&system1,
								   actuators,
								   no_of_actuator_destinations,
								   &random_seed);

	system1.migration_tick=2;

	/* save to file */
	sprintf(filename,"%stestsystem.dat",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"w");
	assert(fp!=0);
	gpr_save_system(&system1,fp);
	fclose(fp);

	/* load from file */
	fp = fopen(filename,"r");
	assert(fp!=0);
	gpr_load_system(&system2,fp,instruction_set,no_of_instructions);
	fclose(fp);

	assert(system1.size==system2.size);
	assert(system1.migration_tick==system2.migration_tick);

	/* check that the functions are the same */
	for (j = 0; j < system1.size; j++) {
		population1 = &system1.island[j];
		population2 = &system2.island[j];
		for (i = 0; i < population1->size; i++) {
			retval = 0;
			/* check the functions */
			gpr_functions_are_equal(&population1->individual[i],
									&population2->individual[i],
									&retval);
			assert (retval==0);
			/* check the state */
			state1 = &population1->state[i];
			state2 = &population2->state[i];
			assert(state1->no_of_sensor_sources==no_of_sensor_sources);
			assert(state1->no_of_actuator_destinations==
				   no_of_actuator_destinations);
			assert(state2->no_of_sensor_sources==no_of_sensor_sources);
			assert(state2->no_of_actuator_destinations==
				   no_of_actuator_destinations);
			for (k = 0; k < sensors; k++) {
				assert(state1->sensor_source[k]==
					   state2->sensor_source[k]);
			}
			for (k = 0; k < actuators; k++) {
				assert(state1->actuator_destination[k]==
					   state2->actuator_destination[k]);
			}
		}
	}

	/* free memory */
	gpr_free_system(&system1);
	gpr_free_system(&system2);

	/* delete the test file */
	sprintf(str,"rm %s",filename);
	retval = system(str);

	printf("Ok\n");
}

void test_gpr_S_expression()
{
	gpr_function f;
	int ctr=0;
	int min_depth=2;
	int max_depth=10,depth=0;
	float branching_prob=0.8f;
	FILE * fp;
	char filename[128],str[128];
	unsigned int random_seed = 123;
	int integers_only = 0;
	int instruction_set[64], no_of_instructions=0;

	printf("test_gpr_S_expression...");

	/* create an instruction set */
	no_of_instructions =
		gpr_default_instruction_set((int*)instruction_set);
	assert(no_of_instructions>0);

	/* create a random tree */
	depth = 0;
	ctr = 0;
	gpr_random(&f,depth,min_depth,max_depth,branching_prob,
			   100,200, integers_only, &random_seed,
			   (int*)instruction_set, no_of_instructions);
	gpr_nodes(&f,&ctr);
	assert(ctr>1);

	/* save as dot file */
	sprintf(filename,"%stest_s.txt",GPR_TEMP_DIRECTORY);
	fp = fopen(filename,"w");
	assert(fp!=0);
	gpr_S_expression(&f,fp);
	fclose(fp);

	gpr_free(&f);

	ctr=0;
	gpr_nodes(&f,&ctr);
	assert(ctr==0);

	/* delete the test file */
	sprintf(str,"rm %s",filename);
	ctr = system(str);

	printf("Ok\n");
}

void test_gpr_data()
{
	gpr_data data;
	unsigned int size = 10;
	unsigned int fields = 3;
	unsigned int field = 1;
	unsigned int i;
	float real=0, imaginary=0;

	printf("test_gpr_data...");

	gpr_data_init(&data, size, fields);

	/* test pushes */
	for (i = 0; i < size*2; i++) {
		if (i < size) {
			assert(data.tail==0);
			assert(data.head==(unsigned short)i);
		}
		else {
			assert(data.head == (unsigned short)(i - size));
			if (i < (size*2)-1) {
				if (data.tail != (unsigned short)(i - size + 1)) {
					printf("\nhead %d  tail %d\n",
						   data.head, data.tail);
				}
				assert(data.tail == (unsigned short)(i - size + 1));
			}
			else {
				assert(data.tail == (unsigned short)0);
			}
		}
		gpr_data_set_head(&data, field, (float)i, (float)i);
		gpr_data_push(&data);
	}

	for (i = 1; i < size-1; i++) {
		gpr_data_get_elem(&data, i, field, &real, &imaginary);
		assert((int)real == size+i+1);
	}

	for (i = 1; i < size; i++) {
		gpr_data_get_tail(&data, field, &real, &imaginary);
		gpr_data_pop(&data);
		assert((int)real == size+i);
	}

    gpr_data_free(&data);

	printf("Ok\n");
}

int run_tests()
{
	printf("Running tests\n");

	test_gpr_data();
	test_rand_num();
	test_gpr_mutate_value();
	test_gpr_random_value();
	test_gpr_init();
	test_gpr_random();
	test_gpr_prune();
	test_gpr_copy();
	test_gpr_mutate();
	test_gpr_crossover();
	test_gpr_mate();
	test_gpr_run();
	test_gpr_sort();
	test_gpr_sort_system();
	test_gpr_init_state();
	test_gpr_generation();
	test_gpr_generation_system();
	test_gpr_dot();
	test_gpr_save_load();
	test_gpr_save_load_population();
	test_gpr_save_load_system();
	test_gpr_S_expression();
	test_gpr_ADF_population();
	test_gpr_environment();

	printf("All tests completed\n");
	return 1;
}

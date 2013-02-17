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

#ifndef GPRCM_H
#define GPRCM_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "globals.h"
#include "gpr.h"
#include "gprc.h"

#define GPRCM_MORPHOLOGY_ROWS                  5
#define GPRCM_MORPHOLOGY_COLUMNS               8
#define GPRCM_MORPHOLOGY_SENSORS               2
#define GPRCM_MORPHOLOGY_ACTUATORS             2
#define GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE  8

struct gprcm_func {
	/* instruction set for the morphology generator */
	int morphology_no_of_instructions;
	int * morphology_instruction_set;

	/* morphology generator */
	gprc_function * morphology;

	/* the main program */
	gprc_function * program;
};
typedef struct gprcm_func gprcm_function;

/* represents a population */
struct gprcm_pop {
	/* the number of individuals in the population */
	int size;
	/* the number of rows and columns in the grid for each individual */
	int rows, columns;
	/* the number of sensors and actuators */
	int sensors, actuators;
	/* the maximum number of connections for each gene in the grid */
	int connections_per_gene;
	/* the number of ADF modules */
	int ADF_modules;
	/* the minimum and maximum constant values */
	float min_value, max_value;
	/* the number of chromosomes for each individual */
	int chromosomes;
	/* whether to only use integer maths */
	int integers_only;
	/* array containing individual programs */
	struct gprcm_func * individual;
	float * fitness;
	/* the fitness history for the population */
	struct gpr_hist history;
};
typedef struct gprcm_pop gprcm_population;

/* system containing a number of populations */
struct gprcm_sys {
	/* the number of sub-populations or islands */
	int size;
	/* the number of time steps after which
	   migrations between islands will occur */
	int migration_tick;
	/* population for each island */
	gprcm_population * island;
	/* the best fitness for each island */
	float * fitness;
	/* the fitness history for the system */
	struct gpr_hist history;
};
typedef struct gprcm_sys gprcm_system;

struct gprcm_env {
	/* the number of individuals in the population */
	int population_size;
	/* the maximum population size */
	int max_population_size;
	/* the number of rows and columns in the grid for each individual */
	int rows, columns;
	/* the number of sensors and actuators */
	int sensors, actuators;
	/* the maximum number of connections for each gene in the grid */
	int connections_per_gene;
	/* the number of ADF modules */
	int ADF_modules;
	/* the minimum and maximum constant values */
	float min_value, max_value;
	/* the number of chromosomes for each individual */
	int chromosomes;
	/* whether to only use integer maths */
	int integers_only;
	/* array containing individual programs */
	struct gprcm_func * individual;
	/* the number of matings */
	int matings;
	/* index numbers of mating parents */
	int * mating;
};
typedef struct gprcm_env gprcm_environment;

void gprcm_init(gprcm_function * f,
				int rows, int columns, int sensors, int actuators,
				int connections_per_gene, int ADF_modules,
				unsigned int * random_seed);
void gprcm_free(gprcm_function * f);

#endif

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

#include "gprcm.h"

void gprcm_init(gprcm_function * f,
				int rows, int columns, int sensors, int actuators,
				int connections_per_gene, int ADF_modules,
				unsigned int * random_seed)
{
	int * instruction_set;

	/* create an array for the instruction set to be
	   used with the morphology generator */
	f->morphology_instruction_set =
		(int*)malloc(sizeof(int)*64);
	instruction_set = f->morphology_instruction_set;

	/* create an instruction set for the morphology generator */
	f->morphology_no_of_instructions =
		gprc_equation_instruction_set(instruction_set);

	/* create the morphology generator with a fixed architecture */
	gprc_init(f->morphology,
			  GPRCM_MORPHOLOGY_ROWS,
			  GPRCM_MORPHOLOGY_COLUMNS,
			  GPRCM_MORPHOLOGY_SENSORS,
			  GPRCM_MORPHOLOGY_ACTUATORS,
			  GPRCM_MORPHOLOGY_CONNECTIONS_PER_GENE,
			  0, random_seed);

	/* create the program */
	gprc_init(f->program,
			  rows, columns, sensors, actuators,
			  connections_per_gene, ADF_modules,
			  random_seed);
}

/* free memory */
void gprcm_free(gprcm_function * f)
{
	gprc_free(f->program);
	gprc_free(f->morphology);
	free(f->morphology_instruction_set);
}

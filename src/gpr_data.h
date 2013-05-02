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

#ifndef GPR_DATA_H
#define GPR_DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GPR_DATA_MAX_ENTRIES  20
#define GPR_DATA_MAX_FIELDS    3

#define GPR_FIELD_POS(index,field,fields) \
	((((index)*(fields))+(field))*2)

struct gpr_data_struct {
	unsigned short size;
	unsigned short fields;
	unsigned short head, tail;
	float * block;
};
typedef struct gpr_data_struct gpr_data;

void gpr_data_init(gpr_data * data,
				   unsigned int size, unsigned int fields);
void gpr_data_clear(gpr_data * data);
void gpr_data_free(gpr_data * data);
void gpr_data_get_head(gpr_data * data,
					   unsigned int field,
					   float * real, float * imaginary);
void gpr_data_get_tail(gpr_data * data,
					   unsigned int field,
					   float * real, float * imaginary);
void gpr_data_push(gpr_data * data);
void gpr_data_pop(gpr_data * data);
void gpr_data_get_elem(gpr_data * data,
					   unsigned int index, unsigned int field,
					   float * real, float * imaginary);
void gpr_data_set_elem(gpr_data * data,
					   unsigned int index, unsigned int field,
					   float real, float imaginary);
void gpr_data_set_head(gpr_data * data,
					   unsigned int field,
					   float real, float imaginary);
void gpr_data_set_tail(gpr_data * data,
					   unsigned int field,
					   float real, float imaginary);

#endif

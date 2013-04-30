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

#include "gpr_data.h"

/* create a data structure */
void gpr_data_init(gpr_data * data,
				   unsigned int size, unsigned int fields)
{
	data->size = size;
	data->fields = fields;
	data->head = 0;
	data->tail = 0;
    if (size > 0) {
		data->block = (float*)malloc(size*fields*2*sizeof(float));
		gpr_data_clear(data);
	}
}

/* clear the data */
void gpr_data_clear(gpr_data * data)
{
	if (data->size > 0) {
		memset((void*)data->block,'\0',
			   data->size*data->fields*2*sizeof(float));
	}
}


/* deallocate data structure */
void gpr_data_free(gpr_data * data)
{
	if (data->size > 0) {
		free(data->block);
	}
}

/* pushes a value to the head */
void gpr_data_push(gpr_data * data, unsigned int field,
				   float real, float imaginary)
{
	unsigned int index =
		GPR_FIELD_POS(data->head,field,data->size,data->fields);
	unsigned short next;
	if (data->size==0) return;
	data->block[index] = real;
	data->block[index+1] = imaginary;

	/* increment the position of the head */
	next = data->head+1;
	if (next >= data->size) {
		next = 0;
	}
	if (next != data->tail) {
		data->head = next;
	}
	else {
		data->tail++;
		if (data->tail >= data->size) {
			data->tail = 0;
		}
		data->head = next;
	}
}

void gpr_data_pop(gpr_data * data, unsigned int field,
				  float * real, float * imaginary)
{
	unsigned int index =
		GPR_FIELD_POS(data->tail,field,data->size,data->fields);
	unsigned short next;
	if (data->size==0) return;
	*real = data->block[index];
	*imaginary = data->block[index+1];

	/* increment the position of the tail */
	if (data->tail == data->head) return;
	next = data->tail+1;
	if (next >= data->size) {
		next = 0;
	}
	if (next != data->head) {
		data->tail = next;
	}
}

static unsigned short get_data_index(unsigned short size,
									 unsigned short head,
									 unsigned short tail,
									 unsigned short field,
									 unsigned short fields,
									 unsigned int index)
{
	unsigned short n;

	if (index == 0) return 0;
	if (head > tail) {
		n = tail + ((unsigned short)index % (head - tail));
	}
	else {
		n = tail + ((unsigned short)index % ((size - tail) + head));
		if (n >= size) {
			n -= size;
		}
	}
	return GPR_FIELD_POS(n,field,size,fields);
}

/* returns a value at the head of the list */
void gpr_data_get_head(gpr_data * data,
					   unsigned int field,
					   float * real, float * imaginary)
{
	unsigned short idx;

	if ((data->size==0) ||
		(data->head == data->tail)) {
		*real = 0;
		*imaginary = 0;
		return;
	}
	idx = data->head*data->fields*2;
	*real = data->block[idx];
	*imaginary = data->block[idx+1];
}

/* returns a value at the tail of the list */
void gpr_data_get_tail(gpr_data * data,
					   unsigned int field,
					   float * real, float * imaginary)
{
	unsigned short idx;

	if ((data->size==0) ||
		(data->head == data->tail)) {
		*real = 0;
		*imaginary = 0;
		return;
	}
	idx = data->tail*data->fields*2;
	*real = data->block[idx];
	*imaginary = data->block[idx+1];
}

/* get an element from the data set */
void gpr_data_get_elem(gpr_data * data,
					   unsigned int index, unsigned int field,
					   float * real, float * imaginary)
{
	unsigned short i;

	if (data->head == data->tail) {
		*real = 0;
		*imaginary = 0;
		return;
	}
	i = get_data_index(data->size, data->head, data->tail,
					   (unsigned short)field, data->fields, index);
	*real = data->block[i];
	*imaginary = data->block[i+1];
}

/* set an element from the data set */
void gpr_data_set_elem(gpr_data * data,
					   unsigned int index, unsigned int field,
					   float real, float imaginary)
{
	unsigned short i;

	if (data->head == data->tail) {
		return;
	}
	i = get_data_index(data->size, data->head, data->tail,
					   (unsigned short)field, data->fields, index);
	data->block[i] = real;
	data->block[i+1] = imaginary;
}

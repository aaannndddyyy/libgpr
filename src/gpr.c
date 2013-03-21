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

#include "gpr.h"

/* clear a complex value */
/*void gpr_clear_value(gpr_value * v)
{
	v->value = 0;
	v->imaginary = 0;
	}*/

float gpr_mutate_value(float value,
					   float percent,
					   unsigned int * random_seed)
{
	return value *
		(100 + (percent*(((rand_num(random_seed)%20000)/
						  10000.0f)-1.0f)))/100.0f;
}

static float gpr_run_function(gpr_function * f, gpr_state * state,
							  int call_depth,
							  float (*custom_function)
							  (float,float,float));

int is_nan(float v)
{
	volatile float d = v;
    return d != d;
}

/* Lehmer random number generator */
int rand_num(unsigned int * seed)
{
	unsigned int v =
		((unsigned long long)(*seed) * 279470273UL) % 4294967291UL;
	if (v==0) v = (int)time(NULL); /* avoid the singularity */
	*seed = v;
	return abs((int)v);
}

/* is the given function type a terminal ? */
static unsigned char is_terminal(unsigned char function_type)
{
	if ((function_type==GPR_FUNCTION_VALUE) ||
		(function_type==GPR_FUNCTION_ARG) ||
		(function_type==GPR_FUNCTION_NONE)) {
		return 1;
	}
	return 0;
}

/* search and replace a function type */
static void gpr_replace(gpr_function * f,
						gpr_state * state,
						int search_function_type,
						int replace_function_type,
						int max_argc)
{
	int i=0;

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			/* if a NONE function is encountered but its parent isn't
			   of type VALUE */
			if ((f->argv[i]->function_type==GPR_FUNCTION_NONE) &&
				(f->function_type!=GPR_FUNCTION_VALUE) &&
				(f->function_type!=GPR_FUNCTION_ARG)) {
				f->argv[i]->function_type=GPR_FUNCTION_VALUE;
			}

			/* search and replace */
			if (f->argv[i]->function_type==search_function_type) {
				f->argv[i]->function_type=replace_function_type;
			}

			/* update arguments */
			if (f->argv[i]->function_type==GPR_FUNCTION_ARG) {
				if (max_argc>0) {
					f->argv[i]->value =
						abs(((int)f->argv[i]->value)%max_argc);
				}
				else {
					f->argv[i]->value = 0;
				}
			}

			/* next level of recursion */
			gpr_replace(f->argv[i], state,
						search_function_type,
						replace_function_type,
						max_argc);
		}
	}
}

/* returns 1 if the given program contains the given function */
static void gpr_contains_func(gpr_function * f,
							  int function_type,
							  int argc,
							  int * found)
{
	int i=0;

	if ((f->function_type == function_type) &&
		(f->argc == argc)) {
		*found = 1;
	}

	if (*found == 1) return;

	for (i = 0; i < f->argc; i++) {
		if ((f->argv[i] != 0) && (*found != 1)) {
			gpr_contains_func(f->argv[i],
							  function_type,
							  argc,
							  found);
		}
	}
}

/* returns 1 if the given program contains the given function */
static int gpr_contains_function(gpr_function * f,
								 int function_type, int argc)
{
	int found = 0;

	gpr_contains_func(f, function_type, argc, &found);
	return found;
}

/* Ensure that ADF calls are consistent with the number of
   ADFs defined */
static void gpr_consistent_ADF_calls(gpr_function * f,
									 gpr_state * state,
									 int no_of_ADFs)
{
	int i=0,ADF_index;
	gpr_function *f2;

	if (no_of_ADFs==0) return;

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			f2 = f->argv[i];

			if (f2->function_type==GPR_FUNCTION_ADF) {
				/* index number of the ADF */
				ADF_index = (abs((int)f2->value))%no_of_ADFs;
				f2->value = ADF_index;
				/* set the number of arguments available
				   to this ADF */
				if ((state->ADF_argc[ADF_index]==0) ||
					((state->ADF_argc[ADF_index]>0) &&
					 (f2->argc<state->ADF_argc[ADF_index]))) {
					state->ADF_argc[ADF_index] = f2->argc;
				}
			}

			/* next level of recursion */
			gpr_consistent_ADF_calls(f2,state,no_of_ADFs);
		}
	}
}



/* Enforces ADF structire at the top level */
void gpr_enforce_ADFs(gpr_function * f, gpr_state * state)
{
	int i;
	gpr_function * f2;

	f->function_type = GPR_TOP_LEVEL_FUNCTION;
	for (i = 0; i < f->argc; i++) {
		/* remember the ADFs so that they can be called later */
		state->ADF[i] = f->argv[i];
		/* clear the number of arguments */
		state->ADF_argc[i] = 0;
		/* make the first level nodes into DEFUN and MAIN */
		if (f->argv[i]!=0) {
			f2 = f->argv[i];
			/* set the value to the ADF index */
			f2->value = i;
			if (i < f->argc-1) {
				/* function definition */
				f2->function_type = GPR_FUNCTION_DEFUN;
			}
			else {
				/* contains the main program */
				f2->function_type = GPR_FUNCTION_MAIN;
				gpr_consistent_ADF_calls(f2,state,f->argc-1);
			}
		}
	}
	/* remaining ADF places are blank */
	while (i < GPR_MAX_ARGUMENTS) {
		state->ADF[i]=0;
		i++;
	}

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i] != 0) {
			f2 = f->argv[i];
			if (i < f->argc-1) {
				/* replace any values with arguments */
				gpr_replace(f2,state,
							GPR_FUNCTION_VALUE,
							GPR_FUNCTION_ARG,state->ADF_argc[i]);
			}
			else {
				gpr_replace(f2,state,
							GPR_FUNCTION_ARG,
							GPR_FUNCTION_VALUE,state->ADF_argc[i]);
			}
		}
	}
}

/* returns a pair of numbers from the arguments of the given program */
static void gpr_number_pair(gpr_function * f, gpr_state * state,
							float * value1, float * value2,
							int call_depth,
							float (*custom_function)(float,float,float))
{
	int i=0;

	*value1 = 0;
	*value2 = 0;

	/* numerator */
	while (i < f->argc/2) {
		if (f->argv[i]!=0) {
			*value1 +=
				gpr_run_function((gpr_function*)f->argv[i],
								 state, call_depth,
								 (*custom_function));
		}
		i++;
	}

	/* denominator */
	while (i < f->argc) {
		if (f->argv[i]!=0) {
			*value2 +=
				gpr_run_function((gpr_function*)f->argv[i],
								 state, call_depth,
								 (*custom_function));
		}
		i++;
	}
}

/* returns a pair of numbers from the arguments of the given program */
static float gpr_sum(gpr_function * f, gpr_state * state,
					 int call_depth,
					 float (*custom_function)(float,float,float))
{
	int i=0;
	float sum=0;

	/* numerator */
	while (i < f->argc) {
		if (f->argv[i]!=0) {
			sum +=
				gpr_run_function((gpr_function*)f->argv[i],
								 state, call_depth,
								 (*custom_function));
		}
		i++;
	}
	return sum;
}

/* saves a C program function header */
static void gpr_function_c(FILE * fp, char * function_name,
								 int argc, char * other_params)
{
	int i;

	fprintf(fp,"static float %s%d(", function_name, argc);
	for (i = 1; i <= argc; i++) {
		fprintf(fp,"float v%d", i);
		if (i < argc) {
			fprintf(fp,"%s",", ");
		}
	}
	fprintf(fp,"%s)\n{\n",other_params);
}

/* saves a string of values, eg. v1+v2+v3 */
static void gpr_chain_c(FILE * fp, int argc, char * symbol)
{
	int i;

	for (i = 1; i <= argc; i++) {
		fprintf(fp,"v%d", i);
		if (i < argc) {
			fprintf(fp,"%s",symbol);
		}
	}
}

/* for functions which require two values */
static void gpr_value_c(FILE * fp, int argc,
							  char * symbol, int num)
{
	int i, start=0, end=argc;

	if (num==1) {
		start=0;
		end = argc/2;
	}
	if (num==2) {
		start = argc/2;
		end = argc;
	}

	fprintf(fp,"%s","(");
	for (i = start; i < end; i++) {
		fprintf(fp,"v%d", i+1);
		if (i < end-1) {
			fprintf(fp,"%s",symbol);
		}
	}
	fprintf(fp,"%s",")");
}

/* C version of the pow function */
static void gpr_pow_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_pow", argc, "");
	fprintf(fp,"%s","  float result = v1;\n");
	fprintf(fp,"%s","  float v = v1;\n");
	fprintf(fp,"%s","  int i,itt = 2+(abs((int)v2)%3);\n\n");
	fprintf(fp,"%s","  for (i=0;i<itt;i++) result *= v;\n");
	fprintf(fp,"%s","  if (!isnan(result)) {\n");
	fprintf(fp,"%s","    return result;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}


static float gpr_pow(gpr_function * f,
					 gpr_state * state, int call_depth,
					 float (*custom_function)(float,float,float))
{
	int i=0,itt;
	float v,v1=0,v2=0;
	float result;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	result = v1;
	v = result;
	itt = 2+(abs((int)v2)%3);

	for (i = 0; i < itt; i++) result *= v;
	if (is_nan(result) == 0) {
		return result;
	}
	return 0;
}

static void gpr_exp_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_exp", argc, "");
	fprintf(fp,"%s","  float v = (float)exp(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_exp(gpr_function * f, gpr_state * state,
					 int call_depth,
					 float (*custom_function)(float,float,float))
{
	float v;
	v = (float)exp(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;
}

static void gpr_min_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_min", argc, "");
	fprintf(fp,"%s","  if (");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s"," < ");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",") {\n");
	fprintf(fp,"%s","    return ");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return ");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_min(gpr_function * f, gpr_state * state,
					 int call_depth,
					 float (*custom_function)(float,float,float))
{
	int i;
	float v, min_value=0;

	for (i = 0; i < f->argc; i++) {
		v = gpr_run_function((gpr_function*)f->argv[i], state,
							 call_depth,
							 (*custom_function));
		if ((i == 0) || (v < min_value)) {
			min_value = v;
		}
	}
	return min_value;
}

static void gpr_max_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_max", argc, "");
	fprintf(fp,"%s","  if (");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s"," > ");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",") {\n");
	fprintf(fp,"%s","    return ");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return ");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_max(gpr_function * f,
					 gpr_state * state, int call_depth,
					 float (*custom_function)(float,float,float))
{
	int i;
	float v, max_value=0;

	for (i = 0; i < f->argc; i++) {
		v = gpr_run_function((gpr_function*)f->argv[i], state,
							 call_depth,
							 (*custom_function));
		if ((i == 0) || (v > max_value)) {
			max_value = v;
		}
	}
	return max_value;
}

static void gpr_sigmoid_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_sigmoid", argc, "");
	fprintf(fp,"%s","  float tot = ");
	gpr_value_c(fp, argc, "+", 0);
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","    return 1.0f / (1.0f + exp(tot));");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_sigmoid(gpr_function * f, gpr_state * state,
						 int call_depth,
						 float (*custom_function)(float,float,float))
{
	int i;
	float tot = 0;

	for (i = 0; i < f->argc; i++) {
		tot += gpr_run_function((gpr_function*)f->argv[i], state,
								call_depth,
								(*custom_function));
	}
	return 1.0f / (1.0f + exp(tot));;
}


static void gpr_negate_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_negate", argc, "");
	fprintf(fp,"%s","  return -(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n");

	fprintf(fp,"%s","}\n\n");
}

static float gpr_negate(gpr_function * f, gpr_state * state,
						int call_depth,
						float (*custom_function)(float,float,float))
{
	float v = gpr_sum(f,state,call_depth,(*custom_function));
	if (is_nan(v)==0) {
		return -v;
	}
	return 0;
}

static void gpr_custom_function_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_custom", argc, "");
	fprintf(fp,"%s","  /* insert your function code here */\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_custom_function(gpr_function * f, gpr_state * state,
								 int call_depth,
								 float (*custom_function)
								 (float,float,float))
{
	return (*custom_function)(f->value,
							  gpr_run_function((gpr_function*)f->argv[0],
											   state, call_depth,
											   (*custom_function)),
							  gpr_run_function((gpr_function*)f->argv[1],
											   state, call_depth,
											   (*custom_function)));
}

static void gpr_average_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_mean", argc, "");
	fprintf(fp,"%s","  return (");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,")/%d;\n",argc);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_average(gpr_function * f, gpr_state * state,
						 int call_depth,
						 float (*custom_function)(float,float,float))
{
	int i;
	float v,sum=0;

	for (i = 0; i < f->argc; i++) {
		sum += gpr_run_function((gpr_function*)f->argv[i], state,
								call_depth,
								(*custom_function));
	}
	v = sum / f->argc;
	if (is_nan(sum)==0) {
		return v;
	}
	return 0;
}

static void gpr_add_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_add", argc, "");
	fprintf(fp,"%s","  return ");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_add(gpr_function * f,
					 gpr_state * state, int call_depth,
					 float (*custom_function)(float,float,float))
{
	float sum=gpr_sum(f,state,call_depth,(*custom_function));

	if (is_nan(sum)==0) {
		return sum;
	}
	return 0;
}

/* TODO */
static void gpr_add_complex_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_add_complex", argc, "");
	fprintf(fp,"%s","  return ");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","}\n\n");
}

/* TODO */
static float gpr_add_complex(gpr_function * f,
							 gpr_state * state,
							 int call_depth,
							 float (*custom_function)
							 (float,float,float))
{
	float sum=gpr_sum(f,state,call_depth,(*custom_function));

	if (is_nan(sum)==0) {
		return sum;
	}
	return 0;
}

static void gpr_subtract_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_subtract", argc, "");
	fprintf(fp,"%s","  return ");
	gpr_chain_c(fp, argc, "-");
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_subtract(gpr_function * f, gpr_state * state,
						  int call_depth,
						  float (*custom_function)(float,float,float))
{
	int i;
	float sum=gpr_run_function((gpr_function*)f->argv[0], state,
							   call_depth,
							   (*custom_function));

	for (i = 1; i < f->argc; i++) {
		sum -= gpr_run_function((gpr_function*)f->argv[i], state,
								call_depth,
								(*custom_function));
	}
	if (is_nan(sum)==0) {
		return sum;
	}
	return 0;
}

/* TODO */
static void gpr_subtract_complex_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_subtract_complex", argc, "");
	fprintf(fp,"%s","  return ");
	gpr_chain_c(fp, argc, "-");
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","}\n\n");
}

/* TODO */
static float gpr_subtract_complex(gpr_function * f,
								  gpr_state * state,
								  int call_depth,
								  float (*custom_function)
								  (float,float,float))
{
	int i;
	float sum=gpr_run_function((gpr_function*)f->argv[0], state,
							   call_depth,
							   (*custom_function));

	for (i = 1; i < f->argc; i++) {
		sum -= gpr_run_function((gpr_function*)f->argv[i], state,
								call_depth,
								(*custom_function));
	}
	if (is_nan(sum)==0) {
		return sum;
	}
	return 0;
}

static void gpr_multiply_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_multiply", argc, "");
	fprintf(fp,"%s","  float v = ");
	gpr_chain_c(fp, argc, "*");
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_multiply(gpr_function * f, gpr_state * state,
						  int call_depth,
						  float (*custom_function)(float,float,float))
{
	int i;
	float product=gpr_run_function((gpr_function*)f->argv[0], state,
								   call_depth,
								   (*custom_function));

	for (i = 1; i < f->argc; i++) {
		product *= gpr_run_function((gpr_function*)f->argv[i], state,
									call_depth,
									(*custom_function));
	}
	if (is_nan(product)==0) {
		return product;
	}
	return 0;
}

/* TODO */
static void gpr_multiply_complex_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_multiply_complex", argc, "");
	fprintf(fp,"%s","  float v = ");
	gpr_chain_c(fp, argc, "*");
	fprintf(fp,"%s",";\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* TODO */
static float gpr_multiply_complex(gpr_function * f,
								  gpr_state * state,
								  int call_depth,
								  float (*custom_function)
								  (float,float,float))
{
	int i;
	float product=gpr_run_function((gpr_function*)f->argv[0], state,
								   call_depth,
								   (*custom_function));

	for (i = 1; i < f->argc; i++) {
		product *= gpr_run_function((gpr_function*)f->argv[i], state,
									call_depth,
									(*custom_function));
	}
	if (is_nan(product)==0) {
		return product;
	}
	return 0;
}

static void gpr_weight_c(FILE * fp, int argc, gpr_function * f)
{
	gpr_function_c(fp, "gpr_weight", argc, "");
	fprintf(fp,     "  float v = v1 * %.10f;\n",f->value);
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_weight(gpr_function * f, gpr_state * state,
						int call_depth,
						float (*custom_function)(float,float,float))
{
	float product = gpr_run_function((gpr_function*)f->argv[0], state,
									 call_depth,
									 (*custom_function)) * f->value;

	if (is_nan(product)==0) {
		return product;
	}
	return 0;
}

static void gpr_divide_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_divide", argc, "");
	fprintf(fp,"%s","  float v,divisor=");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",";\n\n");
	fprintf(fp,"%s","  if (fabs(divisor) > 0.01f) {\n");
	fprintf(fp,"%s","    v = ");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s"," / divisor;\n");
	fprintf(fp,"%s","    if (!isnan(v)) {\n");
	fprintf(fp,"%s","      return v;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_divide(gpr_function * f, gpr_state * state,
						int call_depth,
						float (*custom_function)(float,float,float))
{
	float v,num=0,denom=0;

	gpr_number_pair(f, state, &num, &denom, call_depth,
					(*custom_function));

	if (fabs(denom) > 0.01f) {
		v =	num / denom;
		if (is_nan(v)==0) {
			return v;
		}
	}
	return 0;
}

/* TODO */
static void gpr_divide_complex_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_divide_complex", argc, "");
	fprintf(fp,"%s","  float v,divisor=");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",";\n\n");
	fprintf(fp,"%s","  if (fabs(divisor) > 0.01f) {\n");
	fprintf(fp,"%s","    v = ");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s"," / divisor;\n");
	fprintf(fp,"%s","    if (!isnan(v)) {\n");
	fprintf(fp,"%s","      return v;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* TODO */
static float gpr_divide_complex(gpr_function * f,
								gpr_state * state,
								int call_depth,
								float (*custom_function)
								(float,float,float))
{
	float v,num=0,denom=0;

	gpr_number_pair(f, state, &num, &denom, call_depth,
					(*custom_function));

	if (fabs(denom) > 0.01f) {
		v =	num / denom;
		if (is_nan(v)==0) {
			return v;
		}
	}
	return 0;
}

static void gpr_modulus_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_modulus", argc, "");
	fprintf(fp,"%s","  float v,divisor=");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",";\n\n");
	fprintf(fp,"%s","  if (fabs(divisor) > 0.001f) {\n");
	fprintf(fp,"%s","    v = (float)fmod(");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s",", divisor);\n");
	fprintf(fp,"%s","    if (!isnan(v)) {\n");
	fprintf(fp,"%s","      return v;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_modulus(gpr_function * f, gpr_state * state,
						 int call_depth,
						 float (*custom_function)(float,float,float))
{
	float v,num=0,denom=0;

	gpr_number_pair(f, state, &num, &denom, call_depth,
					(*custom_function));

	if (fabs(denom) > 0.01f) {
		v =	fmod(num, denom);
		if (is_nan(v)==0) {
			return v;
		}
	}
	return 0;
}

static void gpr_floor_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_floor", argc, "");
	fprintf(fp,"%s","  return floor(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_floor(gpr_function * f, gpr_state * state,
					   int call_depth,
					   float (*custom_function)(float,float,float))
{
	float v;

	v = floor(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;
}

static void gpr_abs_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_abs", argc, "");
	fprintf(fp,"%s","  float v = fabs(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_abs(gpr_function * f, gpr_state * state,
					 int call_depth,
					 float (*custom_function)(float,float,float))
{
	float v;

	v = fabs(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;
}

/* TODO */
static void gpr_abs_complex_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_abs_complex", argc, "");
	fprintf(fp,"%s","  float v = fabs(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* TODO */
static float gpr_abs_complex(gpr_function * f, gpr_state * state,
							 int call_depth,
							 float (*custom_function)(float,float,float))
{
	float v;

	v = fabs(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;
}

/* call an automatically defined function */
static float gpr_ADF(gpr_function * f, gpr_state * state,
					 int call_depth,
					 float (*custom_function)(float,float,float))
{
	int i, ADF_index;
	gpr_function * adf;

	if (call_depth >= GPR_MAX_CALL_DEPTH-1) return 0;

	/* index of the ADF at the first depth */
	ADF_index = abs(((int)f->value))%GPR_MAX_ARGUMENTS;

	/* get the ADF function */	
	adf = state->ADF[ADF_index];

	if (adf == 0) return 0;

	/* increment the ADF call depth */
	call_depth++;

	/* arguments of the function are saved in a temporary array */
	for (i = 0; i < f->argc; i++) {
		state->temp_ADF_arg[call_depth][i] =
			gpr_run_function((gpr_function*)f->argv[i], state,
							 call_depth,
							 (*custom_function));
	}

	/* now run the ADF function with those arguments */
	return gpr_run_function(adf, state, call_depth,
							(*custom_function));
}

static void gpr_square_root_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_square_root", argc, "");
	fprintf(fp,"  float v = (float)sqrt(gpr_abs%d(",argc);
	gpr_chain_c(fp, argc, ",");
	fprintf(fp,"%s","));\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_square_root(gpr_function * f, gpr_state * state,
							 int call_depth,
							 float (*custom_function)(float,float,float))
{
	float v = (float)sqrt(gpr_abs(f, state, call_depth,
								  (*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;	
}

/* TODO */
static void gpr_square_root_complex_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_square_root_complex", argc, "");
	fprintf(fp,"  float v = (float)sqrt(gpr_abs%d(",argc);
	gpr_chain_c(fp, argc, ",");
	fprintf(fp,"%s","));\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* TODO */
static float gpr_square_root_complex(gpr_function * f, gpr_state * state,
									 int call_depth,
									 float (*custom_function)(float,float,float))
{
	float v = (float)sqrt(gpr_abs(f, state, call_depth,
								  (*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;	
}

static void gpr_sine_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_sine", argc, "");
	fprintf(fp,"%s","  float v = (float)sin(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_sine(gpr_function * f, gpr_state * state,
					  int call_depth,
					  float (*custom_function)(float,float,float))
{
	float v;

	v = (float)sin(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;	
}

static void gpr_arcsine_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_arcsine", argc, "");
	fprintf(fp,"%s","  float v = (float)asin(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_arcsine(gpr_function * f, gpr_state * state,
						 int call_depth,
						 float (*custom_function)(float,float,float))
{
	float v;

	v = (float)asin(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;
}

static void gpr_cosine_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_cosine", argc, "");
	fprintf(fp,"%s","  float v = (float)cos(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_cosine(gpr_function * f, gpr_state * state,
						int call_depth,
						float (*custom_function)(float,float,float))
{
	float v;

	v = (float)cos(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;
}

static void gpr_arccosine_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_arccosine", argc, "");
	fprintf(fp,"%s","  float v = (float)acos(");
	gpr_chain_c(fp, argc, "+");
	fprintf(fp,"%s",");\n\n");
	fprintf(fp,"%s","  if (!isnan(v)) {\n");
	fprintf(fp,"%s","    return v;\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

static float gpr_arccosine(gpr_function * f, gpr_state * state,
						   int call_depth,
						   float (*custom_function)(float,float,float))
{
	float v;

	v = (float)acos(gpr_sum(f,state,call_depth,(*custom_function)));
	if (is_nan(v)==0) {
		return v;
	}
	return 0;	
}

static void gpr_noop_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_noop", argc, "");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* don't do anything, just run the functions */
static float gpr_noop(gpr_function * f, gpr_state * state,
					  int call_depth,
					  float (*custom_function)(float,float,float))
{
	int i;

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			gpr_run_function((gpr_function*)f->argv[i], state,
							 call_depth,
							 (*custom_function));
		}
	}
	return 0;
}

static void gpr_greater_than_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_greater_than", argc, "");
	fprintf(fp,"%s","  if (");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s"," > ");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",") {\n");
	fprintf(fp,     "    return %d;\n",GPR_TRUE);
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  return %d;\n",GPR_FALSE);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_greater_than(gpr_function * f, gpr_state * state,
							  int call_depth,
							  float (*custom_function)
							  (float,float,float))
{
	float v1=0,v2=0;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	if (v1 > v2) {
		return GPR_TRUE;
	}
	return GPR_FALSE;
}

static void gpr_less_than_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_less_than", argc, "");
	fprintf(fp,"%s","  if (");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s"," < ");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",") {\n");
	fprintf(fp,     "    return %d;\n",GPR_TRUE);
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  return %d;\n",GPR_FALSE);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_less_than(gpr_function * f, gpr_state * state,
						   int call_depth,
						   float (*custom_function)(float,float,float))
{
	float v1=0,v2=0;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	if (v1 < v2) {
		return GPR_TRUE;
	}
	return GPR_FALSE;
}

static void gpr_equals_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_equals", argc, "");
	fprintf(fp,"%s","  if ((int)");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s"," == (int)");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",") {\n");
	fprintf(fp,     "    return %d;\n",GPR_TRUE);
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  return %d;\n",GPR_FALSE);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_equals(gpr_function * f, gpr_state * state,
						int call_depth,
						float (*custom_function)(float,float,float))
{
	float v1=0,v2=0;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	if ((int)v1 == (int)v2) {
		return GPR_TRUE;
	}
	return GPR_FALSE;
}

static void gpr_and_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_and", argc, "");
	fprintf(fp,"%s","  if ((");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s",">0) && (");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",">0)) {\n");
	fprintf(fp,     "    return %d;\n",GPR_TRUE);
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  return %d;\n",GPR_FALSE);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_and(gpr_function * f, gpr_state * state,
					 int call_depth,
					 float (*custom_function)(float,float,float))
{
	float v1=0,v2=0;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	if ((v1>0) && (v2>0)) {
		return GPR_TRUE;
	}
	return GPR_FALSE;
}

static void gpr_or_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_or", argc, "");
	fprintf(fp,"%s","  if ((");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s",">0) || (");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",">0)) {\n");
	fprintf(fp,     "    return %d;\n",GPR_TRUE);
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  return %d;\n",GPR_FALSE);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_or(gpr_function * f,
					gpr_state * state, int call_depth,
					float (*custom_function)(float,float,float))
{
	float v1=0,v2=0;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	if ((v1>0) || (v2>0)) {
		return GPR_TRUE;
	}
	return GPR_FALSE;
}

static void gpr_xor_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_xor", argc, "");
	fprintf(fp,"%s","  if ((");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s",">0) != (");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",">0)) {\n");
	fprintf(fp,     "    return %d;\n",GPR_TRUE);
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  return %d;\n",GPR_FALSE);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_xor(gpr_function * f,
					 gpr_state * state, int call_depth,
					 float (*custom_function)(float,float,float))
{
	float v1=0,v2=0;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	if ((v1>0) != (v2>0)) {
		return GPR_TRUE;
	}
	return GPR_FALSE;
}

static void gpr_not_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_not", argc, "");
	fprintf(fp,"%s","  if ((int)(");
	gpr_value_c(fp, argc, "+", 1);
	fprintf(fp,"%s",") != (int)(");
	gpr_value_c(fp, argc, "+", 2);
	fprintf(fp,"%s",")) {\n");
	fprintf(fp,     "    return %d;\n",GPR_TRUE);
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  return %d;\n",GPR_FALSE);
	fprintf(fp,"%s","}\n\n");
}

static float gpr_not(gpr_function * f,
					 gpr_state * state, int call_depth,
					 float (*custom_function)(float,float,float))
{
	float v1=0,v2=0;

	gpr_number_pair(f, state, &v1, &v2, call_depth,
					(*custom_function));

	if (((int)v1) != ((int)v2)) {
		return GPR_TRUE;
	}
	return GPR_FALSE;
}

static void gpr_store_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_store", argc,"");
	fprintf(fp,"%s","  int index, v;\n\n");
	fprintf(fp,"%s","  v = (int)v1;\n");
	fprintf(fp,"%s","  if (v < 0) {\n");
	fprintf(fp,"%s","    /* set register */\n");
	fprintf(fp,"%s","    if (registers > 0) {\n");
	fprintf(fp,"%s","      index = -v % registers;\n");
	fprintf(fp,"%s","      registerlist[index] = v2;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  else {\n");
	fprintf(fp,"%s","    /* set actuator */\n");
	fprintf(fp,"%s","    if (actuators > 0) {\n");
	fprintf(fp,"%s","      index = v % actuators;\n");
	fprintf(fp,"%s","      actuator[index] = v2;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 1;\n");
	fprintf(fp,"%s","}\n\n");
}

/* stores a value in a register or actuator */
static float store(float v1, float v2, gpr_state * state)
{
	int index, v;

	v = (int)v1;
	
	if (v < 0) {
		/* set register */
		if (state->no_of_registers > 0) {
			index = -v % state->no_of_registers;
			state->registers[index] = v2;
		}
	}
	else {
		/* set actuator */
		if (state->no_of_actuators > 0) {
			index = v % state->no_of_actuators;
			state->actuators[index] = v2;
		}
	}
	return 1;
}

static float gpr_set(gpr_function * f,
					 gpr_state * state, int call_depth,
					 float (*custom_function)(float,float,float))
{
	return store(gpr_run_function((gpr_function*)f->argv[0], state,
								  call_depth,
								  (*custom_function)),
				 gpr_run_function((gpr_function*)f->argv[1], state,
								  call_depth,
								  (*custom_function)),
				 state);
}

static void gpr_fetch_c(FILE * fp, int argc)
{
	gpr_function_c(fp, "gpr_fetch", argc, "");
	fprintf(fp,     "  int oracle_type = (int)v1 %% %d;\n",GPR_ORACLES);
	fprintf(fp,"%s","  int index, v = abs((int)v2);\n\n");
	fprintf(fp,"%s","  switch(oracle_type) {\n");
	fprintf(fp,     "  case %d: {\n",GPR_ORACLE_REGISTER);
	fprintf(fp,"%s","    if (registers > 0) {\n");
	fprintf(fp,"%s","      index = v % registers;\n");
	fprintf(fp,"%s","      return registerlist[index];\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  case %d: {\n",GPR_ORACLE_SENSOR);
	fprintf(fp,"%s","    if (sensors > 0) {\n");
	fprintf(fp,"%s","      index = v % sensors;\n");
	fprintf(fp,"%s","      return sensor[index];\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,     "  case %d: {\n",GPR_ORACLE_ACTUATOR);
	fprintf(fp,"%s","    if (actuators > 0) {\n");
	fprintf(fp,"%s","      index = v % actuators;\n");
	fprintf(fp,"%s","      return actuator[index];\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* retrieves a value from register, sensor or actuator */
static float fetch(float v1, float v2, gpr_state * state)
{
	int oracle_type = (int)v1 % GPR_ORACLES;
	int index, v =	abs((int)v2);

	switch(oracle_type) {
	case GPR_ORACLE_REGISTER: {
		if (state->no_of_registers > 0) {
			index = v % state->no_of_registers;
			return state->registers[index];
		}
	}
	case GPR_ORACLE_SENSOR: {
		if (state->no_of_sensors > 0) {
			index = v % state->no_of_sensors;
			return state->sensors[index];
		}
	}
	case GPR_ORACLE_ACTUATOR: {
		if (state->no_of_actuators > 0) {
			index = v % state->no_of_actuators;
			return state->actuators[index];
		}
	}
	}
	return 0;
}

static float gpr_get(gpr_function * f,
					 gpr_state * state, int call_depth,
					 float (*custom_function)(float,float,float))
{
	return fetch(gpr_run_function((gpr_function*)f->argv[0], state,
								  call_depth,
								  (*custom_function)),
				 gpr_run_function((gpr_function*)f->argv[1], state,
								  call_depth,
								  (*custom_function)),
				 state);
}

/* Used to check that the two functions are the same.
   This should return zero if the given functions are the same */
void gpr_functions_are_equal(gpr_function * f1, gpr_function * f2,
							 int * result)
{
	int i = 0;

	if ((f1==0) || (f2==0)) return;

	if (f1->function_type != f2->function_type) {
		printf("Function type %d %d\n",
			   f1->function_type,
			   f2->function_type);
		*result = -1;
	}

	if (abs((int)(f1->value*1000) - (int)(f2->value*1000))>1) {
		printf("Function Value %d %d\n",
			   (int)(f1->value*1000),
			   (int)(f2->value*1000));
		*result = -2;
	}

	if (f1->argc != f2->argc) {		
		if ((f1->argc > 0) && (f2->argc > 0)) {
			if ((f1->argv[0]!=0) && (f2->argv[0]!=0)) {
				*result = -3;
			}
		}
	}
	else {
		for (i = 0; i < f1->argc; i++) {
			if ((f1->argv[i]!=0) && (f2->argv[i]!=0)) {
				gpr_functions_are_equal((gpr_function*)f1->argv[i],
										(gpr_function*)f2->argv[i],
										result);
			}
		}
	}
}

/* deallocate memory */
void gpr_free(gpr_function * f)
{
	for (int i = 0; i < GPR_MAX_ARGUMENTS; i++) {
		if (f->argv[i]!=0) {
			gpr_free((gpr_function*)f->argv[i]);
			free(f->argv[i]);
			f->argv[i]=0;
		}
	}
	f->function_type = GPR_FUNCTION_NONE;
	f->value = 0;
}


/* returns the number of nodes in the tree */
void gpr_nodes(gpr_function * f, int * ctr)
{
	int i=0;

	if (f->function_type != GPR_FUNCTION_NONE) {
		*ctr = *ctr + 1;
	}

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			gpr_nodes((gpr_function*)f->argv[i], ctr);
		}
	}
}

/* returns the maximum depth of the tree */
void gpr_max_depth(gpr_function * f, int depth, int * max_depth)
{
	int i=0;

	if (depth > *max_depth) *max_depth = depth;

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			gpr_max_depth((gpr_function*)f->argv[i],
						  depth+1, max_depth);
		}
	}
}

/* returns a random value in the given range */
float gpr_random_value(float min_value, float max_value,
					   unsigned int * random_seed)
{
	return min_value + ((max_value - min_value) *
						(rand_num(random_seed)%10000) / 10000.0f);
}

/* initialise sensor sources for the given system */
void gpr_init_sensor_sources(gpr_system * system,
							 int no_of_sensors,
							 int no_of_sensor_sources,
							 unsigned int * random_seed)
{
	int i,j,k;
	gpr_state * state;
	gpr_population * population;

	if (no_of_sensor_sources<=0) return;

	for (i = 0; i < system->size; i++) {
		/* for every island in the system */
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			state = &population->state[j];
			state->no_of_sensor_sources = no_of_sensor_sources;
			/* create the array */
			state->sensor_source =
				(int*)malloc(no_of_sensors*sizeof(int));
			/* set random values */
			for (k = 0; k < no_of_sensors; k++) {
				state->sensor_source[k] =
					rand_num(random_seed)%no_of_sensor_sources;
			}
		}
	}
}							 

/* initialise actuator_detinations for the given system */
void gpr_init_actuator_destinations(gpr_system * system,
									int no_of_actuators,
									int no_of_actuator_destinations,
									unsigned int * random_seed)
{
	int i,j,k;
	gpr_state * state;
	gpr_population * population;

	if (no_of_actuator_destinations<=0) return;

	for (i = 0; i < system->size; i++) {
		/* for every island in the system */
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			state = &population->state[j];
			state->no_of_actuator_destinations =
				no_of_actuator_destinations;
			/* create the array */
			state->actuator_destination =
				(int*)malloc(no_of_actuators*sizeof(int));
			/* initial random values */
			for (k = 0; k < no_of_actuators; k++) {
				state->actuator_destination[k] =
					rand_num(random_seed)%no_of_actuator_destinations;
			}
		}
	}
}							 

/* initialise machine state */
void gpr_init_state(gpr_state * state,
					int registers,
					int sensors,
					int actuators,
					unsigned int * random_seed)
{
	int i;

	state->random_seed = rand_num(random_seed);
	state->age = 0;

	state->registers = 0;
	state->sensors = 0;
	state->actuators = 0;
	state->sensor_source = 0;
	state->actuator_destination = 0;

	/* don't do anything with sensor sources
	   or actuator destinations.  We can enable this later if it
	   is needed */
	state->no_of_registers = registers;
	state->no_of_sensors = sensors;
	state->no_of_actuators = actuators;
	state->no_of_sensor_sources = 0;
	state->no_of_actuator_destinations = 0;

	if (registers>0) {
		state->registers = (float*)malloc(sizeof(float)*registers);
#ifdef DEBUG
		assert(state->registers!=0);
#endif
		memset((void*)state->registers,'\0',sizeof(float)*registers);
	}
	if (sensors>0) {
		state->sensors = (float*)malloc(sizeof(float)*sensors);
#ifdef DEBUG
		assert(state->sensors!=0);
#endif
		memset((void*)state->sensors,'\0',sizeof(float)*sensors);
	}
	if (actuators>0) {
		state->actuators = (float*)malloc(sizeof(float)*actuators);
#ifdef DEBUG
		assert(state->actuators!=0);
#endif
		memset((void*)state->actuators,'\0',sizeof(float)*actuators);
	}
	/* clear pointers to ADFs */
	memset((void*)state->ADF,'\0',
		   sizeof(gpr_function*)*GPR_MAX_ARGUMENTS);
	/* clear temporary arguments */
	memset((void*)state->temp_ADF_arg,'\0',
		   sizeof(float*)*GPR_MAX_CALL_DEPTH);
	/* clear temporary arguments */
	for (i = 0; i < GPR_MAX_CALL_DEPTH; i++) {
		memset((void*)state->temp_ADF_arg[i],'\0',
			   sizeof(float)*GPR_MAX_ARGUMENTS);
	}
}

void gpr_free_state(gpr_state * state)
{
	if (state->no_of_registers>0) {
		free(state->registers);
	}
	if (state->no_of_sensors>0) {
		free(state->sensors);
	}
	if (state->no_of_actuators>0) {
		free(state->actuators);
	}
	if (state->no_of_sensor_sources>0) {
		free(state->sensor_source);
	}
	if (state->no_of_actuator_destinations>0) {
		free(state->actuator_destination);
	}
}

/* initialise a function */
void gpr_init(gpr_function * f)
{
	int i, j;

	f->function_type = GPR_FUNCTION_VALUE;
	f->value = 0;
	f->argc = GPR_DEFAULT_ARGUMENTS;
	for (i = 0; i < GPR_MAX_ARGUMENTS; i++) {
		f->argv[i] = (gpr_function*)malloc(sizeof(gpr_function));
#ifdef DEBUG
		assert(f->argv[i]!=0);
#endif
		f->argv[i]->function_type = GPR_FUNCTION_NONE;
		f->argv[i]->value=0;
		f->argv[i]->argc=GPR_DEFAULT_ARGUMENTS;
		for (j = 0; j < GPR_MAX_ARGUMENTS; j++) {
			f->argv[i]->argv[j]=0;
		}
	}
}

void gpr_copy(gpr_function * source, gpr_function * dest)
{
	int i;

	gpr_init(dest);

	dest->function_type = source->function_type;
	dest->value = source->value;
	dest->argc = source->argc;

	for (i = 0; i < source->argc; i++) {
		if (source->argv[i]!=0) {
			gpr_copy((gpr_function*)source->argv[i],
					 (gpr_function*)dest->argv[i]);
		}
		else {
			if (dest->argv[i]!=0) {
				free(dest->argv[i]);
				dest->argv[i]=0;
			}
		}
	}
}

/* ensure that the program tree doesn't go below a maximum depth */
void gpr_prune(gpr_function * f, int depth, int max_depth,
			   float min_value, float max_value,
			   unsigned int * random_seed)
{
	int i = 0;

	if (depth >= max_depth-1) {
		if (is_terminal(f->function_type)==0) {
			/* remove any subtree */
			gpr_free(f);
			/* make this a terminal */
			gpr_init(f);
			f->value = gpr_random_value(min_value,max_value,
										random_seed);
		}
		return;
	}

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			gpr_prune((gpr_function*)f->argv[i],depth+1,max_depth,
					  min_value, max_value, random_seed);
		}
	}
}

/* validate the given tree */
void gpr_validate(gpr_function * f, int depth,
				  int min_depth, int max_depth,
				  int ADFs,
				  int tree_index, int max_tree_index,
				  int * result)
{
	int i = 0;
	gpr_function * f2;

	/* a terminator was not found at the maximum depth */
	if (depth >= max_depth-1) {
		if (is_terminal(f->function_type)==0) {
			*result = GPR_VALIDATE_TREE_NOT_TERMINATED;
			return;
		}
	}

	/* the maximum depth was exceeded */
	if (depth > max_depth) {
		*result = GPR_VALIDATE_TREE_TOO_DEEP;
		return;
	}

	/* a terminal was found at the root of the program tree */
	if (depth < min_depth) {
		if (is_terminal(f->function_type)!=0) {
			*result = GPR_VALIDATE_TERMINATOR_AT_ROOT;
			return;
		}
	}

	/* if ADFs are being used */
	if ((depth==0) && (ADFs>0)) {
		/* is the top level node what we expect? */
		if (f->function_type!=GPR_TOP_LEVEL_FUNCTION) {
			*result = GPR_VALIDATE_ADF_PATTERN;
			return;
		}
		/* examine the first level nodes */
		for (i = 0; i < f->argc; i++) {
			if (f->argv[i] != 0) {
				f2 = f->argv[i];
				if (i<f->argc-1) {
					/* all but the last node should be of DEFUN type */
					if (f2->function_type!=GPR_FUNCTION_DEFUN) {
						*result = GPR_VALIDATE_ADF_PATTERN;
						return;
					}
				}
				else {
					/* the final node at depth 1 should
					   be of NOOP type */
					if (f2->function_type!=GPR_FUNCTION_MAIN) {
						*result = GPR_VALIDATE_ADF_PATTERN;
						return;
					}
				}
			}
		}
		max_tree_index = f->argc;
	}

	/* examine the next depth */
	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			f2 = f->argv[i];

			/* if ADFs are being used */
			if (ADFs>0) {
				/* if this is the top level then keep track
				   of which branch we are following */
				if (depth == 0) {
					tree_index = i;

					if (tree_index == max_tree_index-1) {
						if (f2->function_type != GPR_FUNCTION_MAIN) {
							*result = GPR_VALIDATE_ADF_PATTERN;
							return;
						}
					}
				}
				else {
					/* if this is a DEFUN node */
					if (tree_index < max_tree_index-1) {
						/* function definition */
						if (f2->function_type ==
							GPR_FUNCTION_VALUE) {
							*result = GPR_VALIDATE_ADF_CONTAINS_VALUE;
							return;
						}
						if (f2->function_type ==
							GPR_FUNCTION_ARG) {
							if ((f2->value>GPR_MAX_ARGUMENTS) ||
								(f2->value<0)) {
								*result =
									GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE;
								return;
							}
						}
					}
					else {
						/* if this is the main program node */
						if (f2->function_type == GPR_FUNCTION_ARG) {
							*result = GPR_VALIDATE_ADF_CONTAINS_ARG;
							return;
						}
						
					}
				}
			}

			/* next level of recursion */
			gpr_validate(f2, depth+1,
						 min_depth, max_depth, ADFs,
						 tree_index, max_tree_index,
						 result);
		}
	}
}

/* returns a positive value if the provided instruction set is valid */
static int gpr_validate_instruction_set(int * instruction_set,
										int no_of_instructions)
{
	int i;

	for (i = 0; i < no_of_instructions; i++) {
		if ((instruction_set[i]==GPR_FUNCTION_HEBBIAN) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_FUNCTION) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_CONSTANT) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_STATE) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_BLOCK) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_CONNECTION1) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_CONNECTION2) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_CONNECTION3) ||
			(instruction_set[i]==GPR_FUNCTION_COPY_CONNECTION4)) {
			return -1;
		}
	}
	return 1;
}

/* creates an instruction set containing mathematical functions
   suitable for producing an equation */
int gpr_equation_instruction_set(int * instruction_set)
{
	instruction_set[0] = GPR_FUNCTION_VALUE;
	instruction_set[1] = GPR_FUNCTION_ADD;
	instruction_set[2] = GPR_FUNCTION_SUBTRACT;
	instruction_set[3] = GPR_FUNCTION_NEGATE;
	instruction_set[4] = GPR_FUNCTION_MULTIPLY;
	instruction_set[5] = GPR_FUNCTION_WEIGHT;
	instruction_set[6] = GPR_FUNCTION_DIVIDE;
	instruction_set[7] = GPR_FUNCTION_MODULUS;
	instruction_set[8] = GPR_FUNCTION_AVERAGE;
	instruction_set[9] = GPR_FUNCTION_NOOP1;
	instruction_set[10] = GPR_FUNCTION_NOOP2;
	instruction_set[11] = GPR_FUNCTION_NOOP3;
	instruction_set[12] = GPR_FUNCTION_NOOP4;
	instruction_set[13] = GPR_FUNCTION_NOOP1;
	instruction_set[14] = GPR_FUNCTION_NOOP2;
	instruction_set[15] = GPR_FUNCTION_NOOP3;
	instruction_set[16] = GPR_FUNCTION_NOOP4;
	instruction_set[17] = GPR_FUNCTION_NOOP1;
	instruction_set[18] = GPR_FUNCTION_NOOP2;
	instruction_set[19] = GPR_FUNCTION_NOOP3;
	instruction_set[20] = GPR_FUNCTION_NOOP4;
	instruction_set[21] = GPR_FUNCTION_EXP;
	instruction_set[22] = GPR_FUNCTION_SQUARE_ROOT;
	instruction_set[23] = GPR_FUNCTION_ABS;
	instruction_set[24] = GPR_FUNCTION_SINE;
	instruction_set[25] = GPR_FUNCTION_ARCSINE;
	instruction_set[26] = GPR_FUNCTION_COSINE;
	instruction_set[27] = GPR_FUNCTION_ARCCOSINE;
	instruction_set[28] = GPR_FUNCTION_POW;
	instruction_set[29] = GPR_FUNCTION_GET;
	instruction_set[30] = GPR_FUNCTION_SET;
	return 31;
}

/* creates a default instruction set containing all possible
   instructions and returns the number of instructions*/
int gpr_default_instruction_set(int * instruction_set)
{
	int i;

	for (i = GPR_TERMINALS; i < GPR_FUNCTION_TYPES_ADVANCED; i++) {
		instruction_set[i-GPR_TERMINALS] = i;
	}
	return GPR_FUNCTION_TYPES_ADVANCED - GPR_TERMINALS;
}

/* creates a simple instruction suitable for FPGAs or GPUs
   and returns the number of instructions*/
int gpr_simple_instruction_set(int * instruction_set)
{
	int i;

	for (i = GPR_TERMINALS; i < GPR_FUNCTION_TYPES_SIMPLE; i++) {
		instruction_set[i-GPR_TERMINALS] = i;
	}
	return GPR_FUNCTION_TYPES_SIMPLE - GPR_TERMINALS;
}

/* returns a random function type from the instruction set */
int gpr_random_function(int * instruction_set, int no_of_instructions,
						unsigned int * random_seed)
{
	return instruction_set[rand_num(random_seed)%no_of_instructions];
}

/* returns the number of arguments for the given function type */
static void gpr_function_args(int function_type,
							  int * min_args, int * max_args)
{
	*min_args = 2;
	*max_args = 2;

	if ((function_type==GPR_FUNCTION_ADD) ||
		(function_type==GPR_FUNCTION_SUBTRACT) ||
		(function_type==GPR_FUNCTION_AVERAGE) ||
		(function_type==GPR_FUNCTION_DIVIDE) ||
		(function_type==GPR_FUNCTION_MULTIPLY) ||
		(function_type==GPR_FUNCTION_MODULUS) ||
		(function_type==GPR_FUNCTION_ABS) ||
		(function_type==GPR_FUNCTION_NOOP1) ||
		(function_type==GPR_FUNCTION_NOOP2) ||
		(function_type==GPR_FUNCTION_NOOP3) ||
		(function_type==GPR_FUNCTION_NOOP4) ||
		(function_type==GPR_FUNCTION_DEFUN) ||
		(function_type==GPR_FUNCTION_MAIN) ||
		(function_type==GPR_FUNCTION_PROGRAM) ||
		(function_type==GPR_FUNCTION_ADF) ||
		(function_type==GPR_FUNCTION_SIGMOID) ||
		(function_type==GPR_FUNCTION_MIN) ||
		(function_type==GPR_FUNCTION_MAX)) {
		*max_args = GPR_MAX_ARGUMENTS;
	}

	if (function_type==GPR_FUNCTION_WEIGHT) {
		*min_args = 1;
		*max_args = 1;
	}
}

/* create a random program */
void gpr_random(gpr_function * f, int depth,
				int min_depth, int max_depth,
				float branching_prob,
				float min_value, float max_value,
				int integers_only,
				unsigned int * random_seed,
				int * instruction_set, int no_of_instructions)
{
	int i,min_args=0,max_args=0;

	gpr_init(f);

	f->value = gpr_random_value(min_value,max_value,random_seed);
	if (integers_only>0) f->value = (int)f->value;

	if ((depth >= max_depth-1) || 
		((depth >= min_depth) &&
		 (rand_num(random_seed)%10000 > branching_prob*10000))) {
		return;
	}

	f->function_type =
		gpr_random_function(instruction_set, no_of_instructions,
							random_seed);
	gpr_function_args(f->function_type, &min_args, &max_args);
	if (max_args > min_args) {
		f->argc = min_args + rand_num(random_seed)%(max_args-min_args);
	}
	else {
		f->argc = min_args;
	}

	if (is_terminal(f->function_type)==0) {
		for (i = 0; i < f->argc; i++) {
			gpr_random((gpr_function*)f->argv[i],depth+1,
					   min_depth,max_depth,
					   branching_prob, min_value, max_value,
					   integers_only, random_seed,
					   instruction_set, no_of_instructions);
		}
	}
}

/* returns the node with the given index number */
static void get_node(gpr_function * f, int index, int * ctr,
					 gpr_function ** node,
					 int depth, int * returned_depth)
{
	int i=0;

	if (f->function_type != GPR_FUNCTION_NONE) {
		if (*ctr == index) {
			*node = (gpr_function*)f;
			*returned_depth = depth;
		}

		*ctr = *ctr + 1;
	}

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			get_node((gpr_function*)f->argv[i], index, ctr, node,
					 depth+1, returned_depth);
		}
	}
}

/* crossover the sensor sources and actuator destinations */
void gpr_crossover_sources(gpr_state * parent1,gpr_state * parent2,
						   gpr_state * child,
						   unsigned int * random_seed)
{
	int crossover_point,i;

	/* cross over the sensor sources */
	if ((parent1->no_of_sensor_sources > 0) &&
		(parent2->no_of_sensor_sources > 0) &&
		(child->no_of_sensor_sources > 0)) {
		crossover_point = rand_num(random_seed)%parent1->no_of_sensors;
		for (i = 0; i < parent1->no_of_sensors; i++) {
			if (i < parent1->no_of_sensors/2) {
				child->sensor_source[crossover_point] =
					parent1->sensor_source[crossover_point];
			}
			else {
				child->sensor_source[crossover_point] =
					parent2->sensor_source[crossover_point];
			}
			crossover_point++;
			if (crossover_point>=parent1->no_of_sensors) {
				crossover_point = 0;
			}
		}		
	}

	/* cross over the actuator destinations */
	if ((parent1->no_of_actuator_destinations > 0) &&
		(parent2->no_of_actuator_destinations > 0) &&
		(child->no_of_actuator_destinations > 0)) {
		crossover_point =
			rand_num(random_seed)%parent1->no_of_actuators;
		for (i = 0; i < parent1->no_of_actuators; i++) {
			if (i < parent1->no_of_actuators/2) {
				child->actuator_destination[crossover_point] =
					parent1->actuator_destination[crossover_point];
			}
			else {
				child->actuator_destination[crossover_point] =
					parent2->actuator_destination[crossover_point];
			}
			crossover_point++;
			if (crossover_point>=parent1->no_of_actuators) {
				crossover_point = 0;
			}
		}		
	}
}

/* A special version of crossover in which ADFs are used */
int gpr_crossover_ADF(gpr_function * parent1, gpr_function * parent2,
					  gpr_function * child, int max_depth,
					  float min_value, float max_value,
					  unsigned int * random_seed)
{
	gpr_function * first_parent, * second_parent;
	gpr_function * first_parent_subtree;
	gpr_function * second_parent_subtree;
	gpr_function * child_subtree;
	int first_parent_node_index, second_parent_node_index, max_argc;
	int node_index;

	/* randomly choose the order of the parents */
	if (rand_num(random_seed)%10000>5000) {
		first_parent = parent1;
		second_parent = parent2;
	}
	else {
		first_parent = parent2;
		second_parent = parent1;
	}

	/* child is initially a copy of the first parent */
	gpr_copy(first_parent,child);

	max_argc = first_parent->argc;
	if (second_parent->argc < max_argc) {
		max_argc = second_parent->argc;
	}

	/* cross over each ADF and the main program */
	for (node_index = 0; node_index < max_argc; node_index++) {

		/* select a first level node */
		if ((node_index < max_argc-1) &&
			(max_argc>1)) {
			/* ADF */
			first_parent_node_index = node_index;
			second_parent_node_index = first_parent_node_index;
		}
		else {
			/* main program */
			first_parent_node_index = first_parent->argc-1;
			second_parent_node_index = second_parent->argc-1;
			node_index = max_argc;
		}

		/* get the equivalent parent subtrees */
		first_parent_subtree =
			first_parent->argv[first_parent_node_index];
		second_parent_subtree =
			second_parent->argv[second_parent_node_index];

		if ((first_parent_subtree==0) || (second_parent_subtree==0)) {
			return 1;
		}

		/* the subtree must be of non-zero size */
		if ((is_terminal(first_parent_subtree->function_type)!=0) ||
			(is_terminal(second_parent_subtree->function_type)!=0)) {
			return 1;
		}

		/* get the child subtree */
		child_subtree = child->argv[first_parent_node_index];

		/* clear the child subtree */
		gpr_free(child_subtree);

		/* cross over the subtrees */
		gpr_crossover(first_parent_subtree, second_parent_subtree,
					  child_subtree, max_depth-1,
					  min_value, max_value,
					  random_seed);
	}
	return 1;
}

/* child is a crossed over version of the two parents */
int gpr_crossover(gpr_function * parent1, gpr_function * parent2,
				  gpr_function * child, int max_depth,
				  float min_value, float max_value,
				  unsigned int * random_seed)
{
	gpr_function * first_parent, * second_parent;
	int crossed=0;
	int first_parent_nodes=0;
	int second_parent_nodes=0;
	int first_parent_crossover_node=0;
	int second_parent_crossover_node=0;
	int ctr=0, depth, returned_depth_child, returned_depth_parent2;
	gpr_function * child_node;
	gpr_function * second_parent_node;

	/* randomly choose the order of the parents */
	if (rand_num(random_seed)%10000>5000) {
		first_parent = parent1;
		second_parent = parent2;
	}
	else {
		first_parent = parent2;
		second_parent = parent1;
	}

	/* child is a copy of the first parent */
	gpr_copy(first_parent,child);

	/* get the number of nodes */
	gpr_nodes(first_parent, &first_parent_nodes);
	gpr_nodes(second_parent, &second_parent_nodes);

	if ((first_parent_nodes>1) && (second_parent_nodes>1)) {
		/* identify crossover points */
		first_parent_crossover_node =
			1 + (rand_num(random_seed)%(first_parent_nodes-1));
		second_parent_crossover_node =
			1 + (rand_num(random_seed)%(second_parent_nodes-1));

		/* get nodes for the two points */
		ctr = 0;
		child_node = 0;
		depth = 0;
		returned_depth_child = 0;
		get_node(child, first_parent_crossover_node, &ctr, &child_node,
				 depth, &returned_depth_child);
		ctr = 0;
		second_parent_node = 0;
		depth = 0;
		returned_depth_parent2 = 0;
		get_node(second_parent, second_parent_crossover_node,
				 &ctr, &second_parent_node,
				 depth, &returned_depth_parent2);

		if ((child_node!=0) && (second_parent_node!=0)) {
			if ((returned_depth_child < max_depth-3) &&
				(returned_depth_parent2 < max_depth-3)) {
				/* don't cross over terminal nodes */
				if ((is_terminal(child_node->function_type)==0) &&
					(is_terminal(second_parent_node->function_type)==0)) {

					/* clear the parent1 subtree */
					gpr_free(child_node);
				
					/* copy the parent2 subtree */
					gpr_copy(second_parent_node,child_node);

					/* ensure that the child tree
					   doesn't exceed the maximum depth */
					depth = 0;
					gpr_prune(child, depth, max_depth,
							  min_value, max_value,
							  random_seed);
					crossed=1;
				}
			}

			/* child is a clone of one parent or the other */
			if (crossed==0) {
				if (rand_num(random_seed)%10000>5000) {
					gpr_free(child);
					gpr_copy(second_parent,child);
				}
			}

			return 1;
		}
	}
	/* child is a clone of one parent or the other */
	if (rand_num(random_seed)%10000>5000) {
		gpr_free(child);
		gpr_copy(second_parent,child);
	}
	return 0;	
}

/* mutates the given state */
void gpr_mutate_state(gpr_state * state,
					  float mutation_prob,
					  unsigned int * random_seed)
{
	int i;
	int prob = (int)(mutation_prob*10000);

	/* mutate sensor sources */
	if (state->no_of_sensor_sources>0) {
		for (i = 0; i < state->no_of_sensors; i++) {
			if (rand_num(random_seed)%10000 < prob) {
				/* point mutation */
				state->sensor_source[i] =
					rand_num(random_seed)%state->no_of_sensor_sources;
			}
		}
	}

	/* mutate actuator destinations */
	if (state->no_of_actuator_destinations>0) {
		for (i = 0; i < state->no_of_actuators; i++) {
			if (rand_num(random_seed)%10000 < prob) {
				/* point mutation */
				state->actuator_destination[i] =
					rand_num(random_seed)%
					state->no_of_actuator_destinations;
			}
		}
	}
}

/* inserts a subtree as an argument */
static void gpr_insert_argument(gpr_function * f,
								int depth,
								int min_depth, int max_depth,
								float min_value, float max_value,
								int integers_only,
								unsigned int * random_seed,
								int * instruction_set,
								int no_of_instructions,
								int index,
								gpr_function * subtree)
{
	int i;
	gpr_function * temp;

	if (f->argc>=GPR_MAX_ARGUMENTS-1) return;

	temp = f->argv[f->argc];

	for (i = f->argc; i >= index+1; i--) {
		f->argv[i] = f->argv[i-1];
	}
	f->argv[index] = temp;

	gpr_free(f->argv[index]);

	if (subtree==NULL) {
		/* if a null subtree is provided then produce
		   a random tree */
		gpr_random(f->argv[index],depth+1,min_depth,max_depth,
				   GPR_DEFAULT_BRANCHING_PROB,
				   min_value,max_value,
				   integers_only,
				   random_seed,
				   instruction_set, no_of_instructions);
	}
	else {
		/* copy the subtree */
		gpr_copy(subtree,f->argv[index]);
	}

	/* ensure that the subtree
	   doesn't exceed the maximum depth */
	gpr_prune(f->argv[index], depth+1, max_depth,
			  min_value, max_value,
			  random_seed);

	f->argc++;
}

/* removes one of the arguments of a function */
static void gpr_remove_argument(gpr_function * f, int index)
{
	gpr_function * victim = 0;
	int i;

	/* maintain a minimum number of arguments */
	if (f->argc <= 2) return;

#ifdef DEBUG
	assert(f->argc < GPR_MAX_ARGUMENTS);
	assert(index < f->argc);
#endif

	if (f->argv[index]!=0) {
		victim = f->argv[index];
		gpr_free(victim);
		victim->function_type = GPR_FUNCTION_VALUE;

		/* shuffle the arguments */
		for (i = index+1; i < f->argc; i++) {
			f->argv[i-1] = f->argv[i];
		}
		f->argc--;
		f->argv[f->argc] = victim;
	}
	else {
		for (i = index+1; i < f->argc; i++) {
			f->argv[i-1] = f->argv[i];
		}
		f->argc--;
		f->argv[f->argc] = 0;
	}
}

void gpr_mutate(gpr_function * f, int depth, int max_depth, float prob,
				float min_value, float max_value,
				int integers_only, int ADFs,
				int tree_index, int max_tree_index,
				unsigned int *random_seed,
				int * instruction_set, int no_of_instructions)
{
	gpr_function *fn,*fn2;
	int i=0,arg_index,min_args=0,max_args=0;
	int min_depth=0,mutation_type;

	if (f==0) return;

	while (i < f->argc) {
		if (f->argv[i]!=0) {
			if (depth==0) {
				/* keep track of the number of nodes at depth 1 */
				tree_index = i;
				max_tree_index = f->argc;
			}
			fn = (gpr_function *)f->argv[i];

			gpr_mutate(fn, depth+1, max_depth, prob,
					   min_value, max_value,
					   integers_only, ADFs,
					   tree_index, max_tree_index, random_seed,
					   instruction_set, no_of_instructions);

			if (rand_num(random_seed)%10000 < prob*10000) {
				
				mutation_type = rand_num(random_seed)%6;
				switch(mutation_type) {
				case 0: { /* point mutations */				   
					if (fn->function_type != GPR_FUNCTION_NONE) {
						if (is_terminal(fn->function_type) == 0) {
							/* Randomly change the function type */
							fn->function_type =
								gpr_random_function(instruction_set,
													no_of_instructions,
													random_seed);

							/* insert ADF calls */
							if ((ADFs>0) && (depth>=1) &&
								(max_tree_index > 1)) {
								if (rand_num(random_seed)%10000 <
									(int)(GPR_ADF_PROB*10000)) {
									fn->function_type = GPR_FUNCTION_ADF;
									fn->value =
										rand_num(random_seed)%
										(max_tree_index-1);
								}
							}

							gpr_function_args(fn->function_type,
											  &min_args, &max_args);
							if (max_args > min_args) {
								fn->argc = min_args +
									rand_num(random_seed)%
									(max_args-min_args);
							}
							else {
								fn->argc = min_args;
							}
						}
						else {
							/* randomly change the terminal value */
							if (integers_only<=0) {
								if (rand_num(random_seed)%10000>5000) {
									/* incremental */
									fn->value =
										gpr_mutate_value(fn->value + 0.01f,
														 GPR_MUTATE_VALUE_PERCENT,
														 random_seed);
								}
								else {
									/* random */
									fn->value =
										gpr_random_value(min_value,
														 max_value,
														 random_seed);
								}
							}
							else {
								if (rand_num(random_seed)%10000>5000) {
									/* incremental */
									if (rand_num(random_seed)%2==0) {
										fn->value = (int)fn->value + 1;
									}
									else {
										fn->value = (int)fn->value - 1;
									}
								}
								else {
									/* random */
									fn->value =
										(int)gpr_random_value(min_value,
															  max_value,
															  random_seed);
								}
							}
							
							if (fn->value > GPR_MAX_CONSTANT) {
								fn->value = GPR_MAX_CONSTANT;
							}
							if (fn->value < -GPR_MAX_CONSTANT) {
								fn->value = -GPR_MAX_CONSTANT;
							}
						}
					}
					break;
				}
				case 1: { /* new random subtree */
					if (is_terminal(fn->function_type) == 0) {
						/* add or change program tree */
						if (depth < max_depth-2) {
							/* delete the existing tree */
							gpr_free(fn);
							/* create a new random tree */
							gpr_random(fn,depth+1,min_depth,max_depth,
									   GPR_DEFAULT_BRANCHING_PROB,
									   min_value,max_value,
									   integers_only,
									   random_seed,
									   instruction_set,
									   no_of_instructions);
						}
					}
					break;
				}
				case 2: { /* permutation mutation */
					if (is_terminal(fn->function_type) == 0) {
						if (depth < max_depth-2) {
							if ((ADFs==0) ||
								((depth >= 1) && (ADFs>0))) {
								/* pick one of the other arguments */
								arg_index =
									rand_num(random_seed)%f->argc;
								if (i != arg_index) {
									/* swap the arguments */
									fn2 = f->argv[arg_index];
									f->argv[i] = fn2;
									f->argv[arg_index] = fn;
								}
							}
						}
					}
					break;
				}
				case 3: { /* tree deletion mutation */
					if (is_terminal(fn->function_type) == 0) {
						if (depth > 0) {							
							gpr_free(fn);
							fn->function_type = GPR_FUNCTION_VALUE;
						}
					}
					break;
				}
				case 4: { /* argument deletion mutation */
					if (is_terminal(fn->function_type) == 0) {
						if (f->argc > 2) {
							gpr_remove_argument(f, i);
						}
					}
					break;
				}
				case 5: { /* mutation inserts a random argument */
					if (is_terminal(fn->function_type) == 0) {
						/* get the range of arguments */
						gpr_function_args(f->function_type,
										  &min_args, &max_args);
						/* less than the max number of arguments? */
						if (f->argc < max_args) {
							arg_index = i;
							/* in the case of ADFs don't insert
							   beyond the main function */
							if ((ADFs>0) && (depth==0) &&
								(arg_index==f->argc-1)) {
								arg_index=0;
							}
							/* insert a random subtree */
							gpr_insert_argument(f,depth,min_depth,
												max_depth,
												min_value, max_value,
												integers_only,
												random_seed,
												instruction_set,
												no_of_instructions,
												arg_index,NULL);
						}
					}
					break;
				}
				}
			}
		}
		i++;
	}
}

/* two parents mate and produce a child */
void gpr_mate(gpr_function * parent1, gpr_function * parent2,
			  int max_depth,
			  float min_value, float max_value,
			  float mutation_prob,
			  float pure_mutant_prob,
			  int integers_only,
			  int ADFs,
			  unsigned int * random_seed,
			  int * instruction_set, int no_of_instructions,
			  gpr_function * child, gpr_state * child_state)
{
	int depth=0;

	if (rand_num(random_seed)%10000 >
		pure_mutant_prob*10000) {
		/* produce some combination of derivative of the parents */
		if (ADFs == 0) {
			/* standard crossover */
			gpr_crossover(parent1, parent2,
						  child, max_depth,
						  min_value, max_value,
						  random_seed);
		}
		else {
			/* ADF version of crossover */
			gpr_crossover_ADF(parent1, parent2,
							  child, max_depth,
							  min_value, max_value,
							  random_seed);
		}
		gpr_mutate(child, depth, max_depth, mutation_prob,
				   min_value, max_value, integers_only,
				   ADFs, 0, 0, random_seed,
				   instruction_set, no_of_instructions);
	}
	else {
		/* introduce a purely random child */
		/*gpr_free(child);*/
		/* create a new random tree */
		gpr_random(child,depth,2,max_depth,GPR_DEFAULT_BRANCHING_PROB,
				   min_value,max_value, integers_only,
				   random_seed, instruction_set, no_of_instructions);
	}
	/* set the top level node function */
	child->function_type = GPR_TOP_LEVEL_FUNCTION;

	/* enforce ADF structure */
	if (ADFs>0) gpr_enforce_ADFs(child, child_state);
}

/* run the given function */
static float gpr_run_function(gpr_function * f,
							  gpr_state * state,
							  int call_depth,
							  float (*custom_function)
							  (float,float,float))
{
	if (f==0) return 0;
	switch(f->function_type) {
	case GPR_FUNCTION_CUSTOM: {
		return gpr_custom_function(f, state, call_depth,
								   (*custom_function));
	}
	case GPR_FUNCTION_VALUE: {
		return f->value;
	}
	case GPR_FUNCTION_ARG: {
		return state->temp_ADF_arg[call_depth][abs((int)f->value)];
	}
	case GPR_FUNCTION_NEGATE: {
		return gpr_negate(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_AVERAGE: {
		return gpr_average(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_POW: {
		return gpr_pow(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_EXP: {
		return gpr_exp(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_SIGMOID: {
		return gpr_sigmoid(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_MIN: {
		return gpr_min(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_MAX: {
		return gpr_max(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_ADD: {
		return gpr_add(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_ADD_COMPLEX: {
		return gpr_add_complex(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_SUBTRACT: {
		return gpr_subtract(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_SUBTRACT_COMPLEX: {
		return gpr_subtract_complex(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_MULTIPLY: {
		return gpr_multiply(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_MULTIPLY_COMPLEX: {
		return gpr_multiply_complex(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_WEIGHT: {
		return gpr_weight(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_DIVIDE: {
		return gpr_divide(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_DIVIDE_COMPLEX: {
		return gpr_divide_complex(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_MODULUS: {
		return gpr_modulus(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_FLOOR: {
		return gpr_floor(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_SQUARE_ROOT: {
		return gpr_square_root(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_SQUARE_ROOT_COMPLEX: {
		return gpr_square_root_complex(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_ABS: {
		return gpr_abs(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_ABS_COMPLEX: {
		return gpr_abs_complex(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_SINE: {
		return gpr_sine(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_ARCSINE: {
		return gpr_arcsine(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_COSINE: {
		return gpr_cosine(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_ARCCOSINE: {
		return gpr_arccosine(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_GREATER_THAN: {
		return gpr_greater_than(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_LESS_THAN: {
		return gpr_less_than(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_EQUALS: {
		return gpr_equals(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_AND: {
		return gpr_and(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_OR: {
		return gpr_or(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_XOR: {
		return gpr_xor(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_NOT: {
		return gpr_not(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_SET: {
		return gpr_set(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_GET: {
		return gpr_get(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_NOOP1: {
		return gpr_noop(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_NOOP2: {
		return gpr_noop(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_NOOP3: {
		return gpr_noop(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_NOOP4: {
		return gpr_noop(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_DEFUN: {
		return gpr_add(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_MAIN: {
		return gpr_noop(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_ADF: {
		return gpr_ADF(f, state, call_depth,(*custom_function));
	}
	case GPR_FUNCTION_PROGRAM: {
		return gpr_noop(f, state, call_depth,(*custom_function));
	}
	}
	return 0;
}

float gpr_run(gpr_function * f, gpr_state * state,
			  float (*custom_function)(float,float,float))
{
	if (state->ADF[0]==0) {
		return gpr_run_function(f, state, 0,
								(*custom_function));
	}
	if (f->argv[f->argc-1]!=0) {
		/* only run the program part of the tree */
		return gpr_run_function(f->argv[f->argc-1], state, 0,
								(*custom_function));
	}
	return 0;
}

/* When the population is initialised make some random
   function calls within the main program */
static void gpr_ADF_calls(gpr_function * f, float prob,
						  unsigned int * random_seed)
{
	int i;

	if (f==0) return;

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			if (is_terminal(f->argv[i]->function_type)==0) {
				if (rand_num(random_seed)%10000 < prob*10000) {
					f->argv[i]->function_type =	GPR_FUNCTION_ADF;
					f->argv[i]->value =
						rand_num(random_seed)%(GPR_MAX_ARGUMENTS-1);
				}
			}
			gpr_ADF_calls((gpr_function*)f->argv[i],
						  prob, random_seed);
		}
	}
}

void gpr_init_population(gpr_population * population,
						 int size,
						 int registers,
						 int sensors,
						 int actuators,
						 int max_tree_depth,
						 float min_value, float max_value,
						 int integers_only,
						 int ADFs,
						 unsigned int * random_seed,
						 int * instruction_set, int no_of_instructions)
{
	int i, j, depth;
	gpr_function * f;
	gpr_state * state;

	/* check that the instruction set is valid */
	if (gpr_validate_instruction_set(instruction_set,
									 no_of_instructions)==-1) {
		printf("\nInstruction set contains invalid functions\n");
		return;
	}

	population->size = size;
	population->history.index = 0;
	population->history.interval = 1;
	population->history.tick = 0;

	/* the program for each individual */
	population->individual =
		(gpr_function*)malloc(sizeof(gpr_function)*size);

	/* the retained state for each individual */
	population->state = (gpr_state*)malloc(sizeof(gpr_state)*size);

	/* fitness values */
	population->fitness = (float*)malloc(sizeof(float)*size);

	/* initialise */
	for (i = 0; i < size; i++) {
		state = &population->state[i];
		gpr_init_state(state,registers,sensors,actuators,
					   random_seed);

		depth = 0;
		f = &population->individual[i];
		gpr_random(f, depth,2,max_tree_depth,
				   GPR_DEFAULT_BRANCHING_PROB,
				   min_value,max_value, integers_only,
				   &state->random_seed,
				   instruction_set, no_of_instructions);

		/* the top node is always of a given type */
		f->function_type = GPR_TOP_LEVEL_FUNCTION;

		/* enforce ADF structure at the top level */
		if (ADFs > 0) {
			for (j = 0; j < f->argc; j++) {
				gpr_ADF_calls(f->argv[j],0.1f,random_seed);
			}
			gpr_enforce_ADFs(f,state);
		}
		else {
			/* clear the ADFs */
			memset((void*)state->ADF,'\0',
				   sizeof(gpr_function*)*GPR_MAX_ARGUMENTS);
		}

		/* fitness has not been evaluated */
		population->fitness[i] = 0;
	}

}

/* initialise an environment with a population */
void gpr_init_environment(gpr_environment * population,
						  int max_population_size,
						  int initial_population_size,
						  int registers,
						  int sensors, int actuators,
						  int max_tree_depth,
						  float min_value, float max_value,
						  int integers_only,
						  int ADFs,
						  unsigned int * random_seed,
						  int * instruction_set,
						  int no_of_instructions)
{
	int i, j, depth;
	gpr_function * f;
	gpr_state * state;

	/* check that the instruction set is valid */
	if (gpr_validate_instruction_set(instruction_set,
									 no_of_instructions)==-1) {
		printf("\nInstruction set contains invalid functions\n");
		return;
	}

	population->max_population_size = max_population_size;
	population->population_size = initial_population_size;

	/* the program for each individual */
	population->individual =
		(gpr_function*)malloc(sizeof(gpr_function)*
							  max_population_size);

	/* the retained state for each individual */
	population->state = (gpr_state*)malloc(sizeof(gpr_state)*
										   max_population_size);

	/* initialise */
	for (i = 0; i < max_population_size; i++) {
		state = &population->state[i];
		gpr_init_state(state,registers,sensors,actuators,
					   random_seed);

		depth = 0;
		f = &population->individual[i];
		gpr_random(f, depth,2,max_tree_depth,
				   GPR_DEFAULT_BRANCHING_PROB,
				   min_value,max_value, integers_only,
				   &state->random_seed,
				   instruction_set, no_of_instructions);

		/* the top node is always of a given type */
		f->function_type = GPR_TOP_LEVEL_FUNCTION;

		/* enforce ADF structure at the top level */
		if (ADFs > 0) {
			for (j = 0; j < f->argc; j++) {
				gpr_ADF_calls(f->argv[j],0.1f,random_seed);
			}
			gpr_enforce_ADFs(f,state);
		}
		else {
			/* clear the ADFs */
			memset((void*)state->ADF,'\0',
				   sizeof(gpr_function*)*GPR_MAX_ARGUMENTS);
		}
	}
	population->matings = 0;
	population->mating =
		(int*)malloc(sizeof(int)*3*max_population_size);
}

/* initialise a system which contains multiple sub-populations */
void gpr_init_system(gpr_system * system,
					 int islands,
					 int population_per_island,
					 int registers,
					 int sensors,
					 int actuators,
					 int max_tree_depth,
					 float min_value, float max_value,
					 int integers_only,
					 int ADFs,
					 unsigned int * random_seed,
					 int * instruction_set, int no_of_instructions)
{
	int i;

	/* check that the instruction set is valid */
	if (gpr_validate_instruction_set(instruction_set,
									 no_of_instructions)==-1) {
		printf("\nInstruction set contains invalid functions\n");
		return;
	}

	system->size = islands;
	system->migration_tick=0;
	system->island =
		(gpr_population*)malloc(islands*sizeof(gpr_population));
	system->fitness = (float*)malloc(islands*sizeof(float));

	/* clear the fitness values */
	for (i = 0; i < islands; i++) {

		/* create a population for the island */
		gpr_init_population(&system->island[i],
							population_per_island,
							registers,
							sensors,
							actuators,
							max_tree_depth,
							min_value, max_value,
							integers_only,
							ADFs,
							random_seed,
							instruction_set, no_of_instructions);

		/* clear the average fitness for the population */
		system->fitness[i] = 0;
	}
}

/* frees memory for a system */
void gpr_free_system(gpr_system * system)
{
	int i;

	for (i = 0; i < system->size; i++) {
		gpr_free_population(&system->island[i]);
	}
	free(system->island);
	free(system->fitness);
}


/* frees memory for a population */
void gpr_free_population(gpr_population * population)
{
	int i;

	for (i = 0; i < population->size; i++) {
		gpr_free(&population->individual[i]);
		gpr_free_state(&population->state[i]);
	}

	free(population->individual);
	free(population->state);
	free(population->fitness);
}

/* frees memory for an environment */
void gpr_free_environment(gpr_environment * population)
{
	int i;

	for (i = 0; i < population->max_population_size; i++) {
		gpr_free(&population->individual[i]);
		gpr_free_state(&population->state[i]);
	}

	free(population->individual);
	free(population->state);
	free(population->mating);
}

/* clear the retained state */
void gpr_clear_state(gpr_state * state)
{
	int i;

	memset((void*)state->registers,'\0',
		   sizeof(float)*state->no_of_registers);
	memset((void*)state->sensors,'\0',
		   sizeof(float)*state->no_of_sensors);
	memset((void*)state->actuators,'\0',
		   sizeof(float)*state->no_of_actuators);
	memset((void*)state->temp_ADF_arg,'\0',
		   sizeof(float*)*GPR_MAX_CALL_DEPTH);
	for (i = 0; i < GPR_MAX_CALL_DEPTH; i++) {
		memset((void*)state->temp_ADF_arg[i],'\0',
			   sizeof(float)*GPR_MAX_ARGUMENTS);
	}
}

/* Evaluates the fitness of all individuals in the population.
   Here we use openmp to speed up the process, since each
   evaluation is independent */
void gpr_evaluate(gpr_population * population,
				  int time_steps, int reevaluate,
				  float (*evaluate_program)(int,
											gpr_function*,
											gpr_state*,int))
{
#pragma omp parallel for
	for (int i = 0; i < population->size; i++) {
		if ((population->fitness[i]==0) ||
			(reevaluate>0)) {
			/* clear the retained state */
			gpr_clear_state(&population->state[i]);

			/* run the evaluation function */
			population->fitness[i] =
				(*evaluate_program)(time_steps,
									&population->individual[i],
									&population->state[i], 0);
		}
		/* population gets older */
		(&population->state[i])->age++;
		if ((&population->state[i])->age > GPR_MAX_AGE) {
			population->fitness[i] = 0;
		}
	}
}

/* evaluates a system containing multiple sub-populations */
void gpr_evaluate_system(gpr_system * system,
						 int time_steps, int reevaluate,
						 float (*evaluate_program)(int,
												   gpr_function*,
												   gpr_state*,int))
{
	int i;

#pragma omp parallel for
	for (i = 0; i < system->size; i++) {
		/* evaluate the island population */
		gpr_evaluate(&system->island[i],
					 time_steps, reevaluate,
					 (*evaluate_program));
		/* set the fitness */
		system->fitness[i] = gpr_average_fitness(&system->island[i]);
	}
}

/* sorts individuals in order of fitness */
void gpr_sort(gpr_population * population)
{
	int i, j, best;
	float max, temp_fitness;
	gpr_function temp_individual;
	gpr_state temp_state;

	for (i = 0; i < population->size-1; i++) {
		max = population->fitness[i];
		best = i;
		for (j = i+1; j < population->size; j++) {
			if (population->fitness[j] > max) {
				max = population->fitness[j];
				best = j;
			}
		}
		if (best != i) {
			/* swap fitness */
			temp_fitness = population->fitness[i];
			population->fitness[i] = population->fitness[best];
			population->fitness[best] = temp_fitness;

			/* swap individual */
			temp_individual = population->individual[i];
			population->individual[i] = population->individual[best];
			population->individual[best] = temp_individual;

			/* swap state */
			temp_state = population->state[i];
			population->state[i] = population->state[best];
			population->state[best] = temp_state;
		}
	}
}

/* sorts populations in order of average fitness */
void gpr_sort_system(gpr_system * system)
{
	int i, j, best;
	float max, temp_fitness;
	gpr_population temp_island;

	for (i = 0; i < system->size-1; i++) {
		max = system->fitness[i];
		best = i;
		for (j = i+1; j < system->size; j++) {
			if (system->fitness[j] > max) {
				max = system->fitness[j];
				best = j;
			}
		}
		if (best != i) {
			/* swap fitness */
			temp_fitness = system->fitness[i];
			system->fitness[i] = system->fitness[best];
			system->fitness[best] = temp_fitness;

			/* swap island */
			temp_island = system->island[i];
			system->island[i] = system->island[best];
			system->island[best] = temp_island;
		}
	}
}

/* Returns a fitness histogram for the given population */
static void gpr_fitness_histogram(gpr_population * population,
								  int *histogram, int histogram_levels,
								  float *min_fitness,
								  float *max_fitness)
{
	int i, index;

	*min_fitness=999999;
	*max_fitness=-999999;

	/* get the minimum and maximum fitness values */
	for (i=0;i<population->size;i++) {
		if (population->fitness[i]>*max_fitness) {
			*max_fitness = population->fitness[i];
		}
		if ((population->fitness[i]<*min_fitness) &&
			(population->fitness[i]>0)) {
			*min_fitness = population->fitness[i];
		}
	}

	/* clear the histogram */
	memset((void*)histogram,'\0',histogram_levels*sizeof(int));

	if (*max_fitness <= *min_fitness) return;

	for (i = 0; i < population->size; i++) {
		if (population->fitness[i]>0) {
			index =
				(int)((population->fitness[i] - *min_fitness)*
					  (histogram_levels-1) /
					  (*max_fitness - *min_fitness));
			histogram[index]++;
		}
	}	
}

/* Returns a fitness histogram for the given population */
static void gpr_fitness_histogram_system(gpr_system * system,
										 int *histogram,
										 int histogram_levels,
										 float *min_fitness,
										 float *max_fitness)
{
	int i, j, index;
	gpr_population * population;

	*min_fitness=999999;
	*max_fitness=-999999;

	/* get the minimum and maximum fitness values */
	for (i = 0; i < system->size; i++) {
		population = &system->island[i];
		for (j=0;j<population->size;j++) {
			if (population->fitness[j]>*max_fitness) {
				*max_fitness = population->fitness[j];
			}
			if ((population->fitness[j]<*min_fitness) &&
				(population->fitness[j]>0)) {
				*min_fitness = population->fitness[j];
			}
		}
	}

	/* clear the histogram */
	memset((void*)histogram,'\0',histogram_levels*sizeof(int));

	if (*max_fitness <= *min_fitness) return;

	for (i = 0; i < system->size; i++) {
		population = &system->island[i];
		for (j = 0; j < population->size; j++) {
			if (population->fitness[j]>0) {
				index =
					(int)((population->fitness[j] - *min_fitness)*
						  (histogram_levels-1) /
						  (*max_fitness - *min_fitness));			
				histogram[index]++;
			}
		}
	}
}

/* returns a value in the range 0.0 - 1.0 indicating the diversity within the population */
static float gpr_diversity(gpr_population * population)
{
	int i, hits=0, histogram[GPR_HISTOGRAM_LEVELS];
	float min_fitness=0, max_fitness=0;
	float average=0,variance=0;
	float occupied_fraction=0;

	gpr_fitness_histogram(population,
						  histogram, GPR_HISTOGRAM_LEVELS,
						  &min_fitness, &max_fitness);

	/* calculate the average and occupied fraction */
	for (i = 0; i < GPR_HISTOGRAM_LEVELS; i++) {
		if (histogram[i]>0) {
			average += histogram[i];
			occupied_fraction++;
			hits++;
		}
	}
	if (hits==0) return 0;

	average /= hits;
	occupied_fraction /= GPR_HISTOGRAM_LEVELS;

	/* calculate the variance */
	if (average>0) {
		for (i = 0; i < GPR_HISTOGRAM_LEVELS; i++) {
			if (histogram[i]>0) {
				variance += fabs(histogram[i] - average);
			}
		}
		if (hits>0) {
			variance = (variance/(float)hits)/average;
		}
	}

	return occupied_fraction * (1.0f/(1.0f+variance));
}

/* Kills an individual with the given array index within
   the given environment.
   This really just swaps the pointers for maximum efficiency */
void gpr_death(gpr_environment * population,
			   int victim_index)
{
	if ((victim_index < 0) ||
		(victim_index >= population->population_size)) {
		return;
	}

	if (population->population_size > 1) {
		gpr_free(&population->individual[victim_index]);
		gpr_copy(&population->individual[population->population_size-1],
				 &population->individual[victim_index]);
	}

	population->population_size--;
}

/* Two parents mate and produce an offspring.
   Here the array indexes of the parents are given and the child index
   is returned.  Indexes of parents and children are stored within
   the mating array */
int gpr_mate_environment(gpr_environment * population,
						 int parent1_index, int parent2_index,
						 int max_tree_depth,
						 float min_value, float max_value,
						 float mutation_prob,
						 float pure_mutant_prob,
						 int integers_only, int ADFs,
						 int * instruction_set, int no_of_instructions)
{
	int child_index;
	gpr_state * state;
	unsigned int * random_seed;
	gpr_function * parent1, * parent2;

	/* has the maximum population size been reached? */
	if (population->population_size >=
		population->max_population_size) {
		return -1;
	}

	/* is the mating array full? */
	if (population->matings >=
		population->max_population_size) {
		printf("Too many matings.  Matings value should be cleared.\n");
		return -1;
	}	

	child_index = population->population_size;
	state = &population->state[child_index];
	random_seed = &state->random_seed;
	parent1 = &population->individual[parent1_index];
	parent2 = &population->individual[parent2_index];

	/* clear the child */
	gpr_free(&population->individual[child_index]);

	/* produce a new child */
	gpr_mate(parent1, parent2,
			 max_tree_depth,
			 min_value, max_value,
			 mutation_prob,
			 pure_mutant_prob,
			 integers_only,
			 ADFs,
			 random_seed,
			 instruction_set, no_of_instructions,
			 &population->individual[child_index],
			 &population->state[child_index]);

	/* cross over the sensor sources and actuator destinations */
	gpr_crossover_sources(&population->state[parent1_index],
						  &population->state[parent2_index],
						  &population->state[child_index],
						  random_seed);
	
	/* mutate the child's sensor sources and actuator destinations */
	gpr_mutate_state(&population->state[child_index],
					 mutation_prob,
					 random_seed);

	/* just born */
	state->age = 0;

	/* store the indexes of the parents */
	population->mating[population->matings*3] = parent1_index;
	population->mating[population->matings*3+1] = parent2_index;
	/* store the index of the child */
	population->mating[population->matings*3+2] =
		population->population_size;
	/* increment the number of matings during this cycle */
	population->matings++;
	/* increment the size of the population */
	population->population_size++;

	return child_index;
}

/* Produce the next generation.
   This assumes that fitness has already been evaluated */
void gpr_generation(gpr_population * population,
					float elitism,
					int max_tree_depth,
					float min_value, float max_value,
					float mutation_prob,
					float pure_mutant_prob,
					int integers_only,
					int ADFs,
					int * instruction_set, int no_of_instructions)
{
	int i, threshold;
	float diversity,mutation_prob_range;

	/* sort the population in order of fitness */
	gpr_sort(population);

	diversity = gpr_diversity(population);
	mutation_prob_range = (1.0f-mutation_prob)/2;
	mutation_prob +=
		mutation_prob_range -
		(mutation_prob_range*diversity);	

	/* store the fitness history */
	population->history.tick++;
	if (population->history.tick >= population->history.interval) {
		population->history.log[population->history.index] =
			population->fitness[0];
		population->history.average[population->history.index] =
			gpr_average_fitness(population);
		population->history.diversity[population->history.index] =
			gpr_diversity(population)*100;
		population->history.index++;
		population->history.tick = 0;

		if (population->history.index >= GPR_MAX_HISTORY-2) {
			for (i = 0; i < GPR_MAX_HISTORY/2; i++) {
				population->history.log[i] =
					population->history.log[i*2];
				population->history.average[i] =
					population->history.average[i*2];
				population->history.diversity[i] =
					population->history.diversity[i*2];
			}
			population->history.index /= 2;
			population->history.interval*=2;
		}
	}

	/* range checking */
	if ((elitism < 0.1f) || (elitism > 0.9f)) {
		elitism = 0.3f;
	}

	/* index setting the threshold for the fittest individuals */
	threshold = (int)((1.0f - elitism)*(population->size-1));

	/*#pragma omp parallel for*/
	for (i = 0; i < population->size - threshold; i++) {
		/* randomly choose parents from the fittest
		   section of the population */
		gpr_state * state = &population->state[threshold + i];
		unsigned int * random_seed = &state->random_seed;
		unsigned int index1 = rand_num(random_seed)%threshold;
		unsigned int index2 = rand_num(random_seed)%threshold;
		gpr_function * parent1 = &population->individual[index1];
		gpr_function * parent2 = &population->individual[index2];

		/* just born */
		state->age = 0;

		/* clear the child */
		gpr_free(&population->individual[threshold + i]);

		/* produce a new child */
		gpr_mate(parent1, parent2,
				 max_tree_depth,
				 min_value, max_value,
				 mutation_prob,
				 pure_mutant_prob,
				 integers_only,
				 ADFs,
				 &population->state[i].random_seed,
				 instruction_set, no_of_instructions,
				 &population->individual[threshold + i],
				 &population->state[threshold + i]);

		/* cross over the sensor sources and actuator destinations */
		gpr_crossover_sources(&population->state[index1],
							  &population->state[index2],
							  &population->state[threshold+i],
							  random_seed);

		/* mutate the child's sensor sources and actuator destinations */
		gpr_mutate_state(&population->state[threshold+i],
						 mutation_prob,
						 random_seed);

		/* fitness not yet evaluated */
		population->fitness[threshold + i] = 0;
	}
}

/* Produce the next generation for a system containing multiple
   sub-populations.
   This assumes that fitness has already been evaluated */
void gpr_generation_system(gpr_system * system,
						   int migration_interval,
						   float elitism,
						   int max_tree_depth,
						   float min_value, float max_value,
						   float mutation_prob,
						   float pure_mutant_prob,
						   int integers_only,
						   int ADFs,
						   int * instruction_set,
						   int no_of_instructions)
{
	int i, depth, migrant_index;
	int island1_index, island2_index;
	unsigned int * random_seed;
	gpr_population * population1, * population2;
	gpr_function * f;
	gpr_state * state;

#pragma omp parallel for
	for (i = 0; i < system->size; i++) {
		gpr_generation(&system->island[i],
					   elitism,
					   max_tree_depth,
					   min_value, max_value,
					   mutation_prob,
					   pure_mutant_prob*
					   (1+(i*pure_mutant_prob/system->size)),
					   integers_only,
					   ADFs,
					   instruction_set, no_of_instructions);		
	}

	/* sort by average fitness */
	gpr_sort_system(system);

	/* migrate individuals between islands */
	system->migration_tick--;
	if (system->migration_tick <= 0) {
		/* reset the counter */
		system->migration_tick = migration_interval;

		random_seed = &((&system->island[0])->state[0].random_seed);

		island1_index = rand_num(random_seed)%system->size; 
		island2_index = rand_num(random_seed)%system->size; 
		if ((island1_index != island2_index)) {
			/* migrate */
			population1 = &system->island[island1_index];
			population2 = &system->island[island2_index];

			/* pick a migrant */
			migrant_index = rand_num(random_seed)%population1->size;

			/* copy it to the other island */
			gpr_free(&population2->individual[population2->size-1]);

			gpr_copy(&population1->individual[migrant_index],
					 &population2->individual[population2->size-1]);

			population2->fitness[population2->size-1] =
				population1->fitness[migrant_index];

			/* free the original */
			gpr_free(&population1->individual[migrant_index]);
			population1->fitness[migrant_index] = 0;
			depth=0;
			/* replace the original with a new random individual */
			gpr_random(&population1->individual[migrant_index],
					   depth, 2, max_tree_depth,
					   GPR_DEFAULT_BRANCHING_PROB,
					   min_value, max_value, integers_only,
					   &population1->state[migrant_index].random_seed,
					   instruction_set, no_of_instructions);
			/* enforce ADF structure */
			f = &population1->individual[migrant_index];
			state = &population1->state[migrant_index];
			f->function_type = GPR_TOP_LEVEL_FUNCTION;
			if (ADFs > 0) {
				for (i = 0; i < f->argc; i++) {
					gpr_ADF_calls(f->argv[i],0.1f,random_seed);
				}
				gpr_enforce_ADFs(f,state);
			}
			else {
				/* clear the ADFs */
				memset((void*)state->ADF,'\0',
					   sizeof(gpr_function*)*GPR_MAX_ARGUMENTS);
			}
		}
	}
}

/* returns the highest fitness value */
float gpr_best_fitness(gpr_population * population)
{
	return population->fitness[0];
}

/* returns the median fitness value */
float gpr_median_fitness(gpr_population * population)
{
	return population->fitness[population->size/2];
}

/* returns the highest fitness value for the given system */
float gpr_best_fitness_system(gpr_system * system)
{
	return gpr_best_fitness(&system->island[0]);
}

/* returns the lowest fitness value */
float gpr_worst_fitness(gpr_population * population)
{
	return population->fitness[population->size-1];
}

/* returns the average fitness of the population */
float gpr_average_fitness(gpr_population * population)
{
	int i;
	float av = 0;

	for (i = 0; i < population->size; i++) {
		av += population->fitness[i];
	}
	if (population->size>0) {
		return av / population->size;
	}
	else {
		printf("Population size is zero\n");
		return 0;
	}
}


/* returns the fittest individual in the population */
gpr_function * gpr_best_individual(gpr_population * population)
{
	return &population->individual[0];
}

/* returns the fittest individual in the given system */
gpr_function * gpr_best_individual_system(gpr_system * system)
{
	return gpr_best_individual(&system->island[0]);
}

/* set a sensor to the given value */
void gpr_set_sensor(gpr_state * state, int index, float value)
{
	state->sensors[index] = value;
}

/* get a sensor value */
float gpr_get_sensor(gpr_state * state, int index)
{
	return state->sensors[index];
}

/* returns the sensor source identifier for the given sensor */
int gpr_get_sensor_source(gpr_state * state, int index)
{
	return state->sensor_source[index];
}

/* returns the value of an actuator */
float gpr_get_actuator(gpr_state * state, int index)
{
	return state->actuators[index];
}

/* returns the actuator destination identifier for the given actuator */
int gpr_get_actuator_destination(gpr_state * state, int index)
{
	return state->actuator_destination[index];
}

/* returns the name of the given function */
void gpr_get_function_name(int function_type, float constant_value,
						   int blank_list,
						   char * name)
{
	switch(function_type) {
	case GPR_FUNCTION_CUSTOM: { sprintf(name,"%s","cust"); return; }
	case GPR_FUNCTION_NEGATE: { sprintf(name,"%s","neg"); return; }
	case GPR_FUNCTION_AVERAGE: { sprintf(name,"%s","mean"); return; }
	case GPR_FUNCTION_POW: { sprintf(name,"%s","pow"); return; }
	case GPR_FUNCTION_EXP: { sprintf(name,"%s","exp"); return; }
	case GPR_FUNCTION_SIGMOID: { sprintf(name,"%s","sig"); return; }
	case GPR_FUNCTION_MIN: { sprintf(name,"%s","min"); return; }
	case GPR_FUNCTION_MAX: { sprintf(name,"%s","max"); return; }
	case GPR_FUNCTION_VALUE: {
		sprintf(name,"%.2f",constant_value);
		return;
	}
	case GPR_FUNCTION_NOOP1: {
		if (blank_list==0) {
			sprintf(name,"%s","list");
		}
		else {
			sprintf(name,"%s"," ");
		}
		return;
	}
	case GPR_FUNCTION_NOOP2: {
		if (blank_list==0) {
			sprintf(name,"%s","list");
		}
		else {
			sprintf(name,"%s"," ");
		}
		return;
	}
	case GPR_FUNCTION_NOOP3: {
		if (blank_list==0) {
			sprintf(name,"%s","list");
		}
		else {
			sprintf(name,"%s"," ");
		}
		return;
	}
	case GPR_FUNCTION_NOOP4: {
		if (blank_list==0) {
			sprintf(name,"%s","list");
		}
		else {
			sprintf(name,"%s"," ");
		}
		return;
	}
	case GPR_FUNCTION_SQUARE_ROOT: {
		sprintf(name,"%s","sqrt");
		return;
	}
	case GPR_FUNCTION_SQUARE_ROOT_COMPLEX: {
		sprintf(name,"%s","sqrti");
		return;
	}
	case GPR_FUNCTION_ABS: { sprintf(name,"%s","abs"); return; }
	case GPR_FUNCTION_ABS_COMPLEX: { sprintf(name,"%s","absi"); return; }
	case GPR_FUNCTION_SINE: { sprintf(name,"%s","sin"); return; }
	case GPR_FUNCTION_ARCSINE: { sprintf(name,"%s","asin"); return; }
	case GPR_FUNCTION_COSINE: { sprintf(name,"%s","cos"); return; }
	case GPR_FUNCTION_ARCCOSINE: { sprintf(name,"%s","acos"); return; }
	case GPR_FUNCTION_ADD: { sprintf(name,"%s","+"); return; }
	case GPR_FUNCTION_ADD_COMPLEX: { sprintf(name,"%s","+i"); return; }
	case GPR_FUNCTION_SUBTRACT: { sprintf(name,"%s","-"); return; }
	case GPR_FUNCTION_SUBTRACT_COMPLEX: { sprintf(name,"%s","-i"); return; }
	case GPR_FUNCTION_MULTIPLY: { sprintf(name,"%s","*"); return; }
	case GPR_FUNCTION_MULTIPLY_COMPLEX: { sprintf(name,"%s","*i"); return; }
	case GPR_FUNCTION_WEIGHT: { sprintf(name,"%s","w"); return; }
	case GPR_FUNCTION_DIVIDE: { sprintf(name,"%s","/"); return; }
	case GPR_FUNCTION_DIVIDE_COMPLEX: { sprintf(name,"%s","/i"); return; }
	case GPR_FUNCTION_MODULUS: { sprintf(name,"%s","mod"); return; }
	case GPR_FUNCTION_FLOOR: { sprintf(name,"%s","floor"); return; }
	case GPR_FUNCTION_GREATER_THAN: { sprintf(name,"%s",">"); return; }
	case GPR_FUNCTION_LESS_THAN: { sprintf(name,"%s","<"); return; }
	case GPR_FUNCTION_EQUALS: { sprintf(name,"%s","="); return; }
	case GPR_FUNCTION_AND: { sprintf(name,"%s","and"); return; }
	case GPR_FUNCTION_OR: { sprintf(name,"%s","or"); return; }
	case GPR_FUNCTION_XOR: { sprintf(name,"%s","xor"); return; }
	case GPR_FUNCTION_NOT: { sprintf(name,"%s","not"); return; }
	case GPR_FUNCTION_SET: { sprintf(name,"%s","store"); return; }
	case GPR_FUNCTION_GET: { sprintf(name,"%s","fetch"); return; }
	case GPR_FUNCTION_COPY_FUNCTION: {
		sprintf(name,"%s","cpfn");
		return;
	}
	case GPR_FUNCTION_COPY_CONSTANT: {
		sprintf(name,"%s","cpco");
		return;
	}
	case GPR_FUNCTION_COPY_STATE: { sprintf(name,"%s","cpst"); return; }
	case GPR_FUNCTION_COPY_BLOCK: { sprintf(name,"%s","cpbl"); return; }
	case GPR_FUNCTION_COPY_CONNECTION1: {
		sprintf(name,"%s","con1");
		return;
	}
	case GPR_FUNCTION_COPY_CONNECTION2: {
		sprintf(name,"%s","con1");
		return;
	}
	case GPR_FUNCTION_COPY_CONNECTION3: {
		sprintf(name,"%s","con1");
		return;
	}
	case GPR_FUNCTION_COPY_CONNECTION4: {
		sprintf(name,"%s","con1");
		return;
	}
	case GPR_FUNCTION_HEBBIAN: { sprintf(name,"%s","hebb"); return; }
	case GPR_FUNCTION_DEFUN: { sprintf(name,"defun %d",
									   (int)constant_value); return; }
	case GPR_FUNCTION_MAIN: { sprintf(name,"%s","main"); return; }
	case GPR_FUNCTION_ADF: { sprintf(name,"ADF %d",
									 1+(int)constant_value); return; }
	case GPR_FUNCTION_ARG: { sprintf(name,"arg %d",
									 ((int)constant_value)%
									 GPR_MAX_ARGUMENTS); return; }
	case GPR_FUNCTION_PROGRAM: { sprintf(name,"%s","progn"); return; }
	}
}

/* returns the name of the given C function used in arduino programs */
static void get_c_function_name(gpr_function * f, char * name)
{
	switch(f->function_type) {
	case GPR_FUNCTION_CUSTOM: {
		sprintf(name,"%s","gpr_custom");
		return;
	}
	case GPR_FUNCTION_NEGATE: {
		sprintf(name,"%s","gpr_negate");
		return;
	}
	case GPR_FUNCTION_AVERAGE: {
		sprintf(name,"%s","gpr_mean");
		return;
	}
	case GPR_FUNCTION_POW: { sprintf(name,"%s","gpr_pow"); return; }
	case GPR_FUNCTION_EXP: { sprintf(name,"%s","gpr_exp"); return; }
	case GPR_FUNCTION_SIGMOID: {
		sprintf(name,"%s","gpr_sigmoid");
		return;
	}
	case GPR_FUNCTION_MIN: { sprintf(name,"%s","gpr_min"); return; }
	case GPR_FUNCTION_MAX: { sprintf(name,"%s","gpr_max"); return; }
	case GPR_FUNCTION_NONE: { sprintf(name,"%s","0"); return; }
	case GPR_FUNCTION_VALUE: { sprintf(name,"%.10f",f->value); return; }
	case GPR_FUNCTION_NOOP1: { sprintf(name,"%s","gpr_noop"); return; }
	case GPR_FUNCTION_NOOP2: { sprintf(name,"%s","gpr_noop"); return; }
	case GPR_FUNCTION_NOOP3: { sprintf(name,"%s","gpr_noop"); return; }
	case GPR_FUNCTION_NOOP4: { sprintf(name,"%s","gpr_noop"); return; }
	case GPR_FUNCTION_SQUARE_ROOT: {
		sprintf(name,"%s","gpr_square_root");
		return;
	}
	case GPR_FUNCTION_SQUARE_ROOT_COMPLEX: {
		sprintf(name,"%s","gpr_square_root_complex");
		return;
	}
	case GPR_FUNCTION_ABS: { sprintf(name,"%s","gpr_abs"); return; }
	case GPR_FUNCTION_ABS_COMPLEX: { sprintf(name,"%s","gpr_abs_complex"); return; }
	case GPR_FUNCTION_SINE: { sprintf(name,"%s","gpr_sine"); return; }
	case GPR_FUNCTION_ARCSINE: {
		sprintf(name,"%s","gpr_arcsine");
		return;
	}
	case GPR_FUNCTION_COSINE: {
		sprintf(name,"%s","gpr_cosine");
		return;
	}
	case GPR_FUNCTION_ARCCOSINE: {
		sprintf(name,"%s","gpr_arccosine");
		return;
	}
	case GPR_FUNCTION_ADD: { sprintf(name,"%s","gpr_add"); return; }
	case GPR_FUNCTION_ADD_COMPLEX: { sprintf(name,"%s","gpr_add_complex"); return; }
	case GPR_FUNCTION_SUBTRACT: {
		sprintf(name,"%s","gpr_subtract");
		return;
	}
	case GPR_FUNCTION_SUBTRACT_COMPLEX: {
		sprintf(name,"%s","gpr_subtract_complex");
		return;
	}
	case GPR_FUNCTION_MULTIPLY: {
		sprintf(name,"%s","gpr_multiply");
		return;
	}
	case GPR_FUNCTION_MULTIPLY_COMPLEX: {
		sprintf(name,"%s","gpr_multiply_complex");
		return;
	}
	case GPR_FUNCTION_WEIGHT: {
		sprintf(name,"%s","gpr_weight");
		return;
	}
	case GPR_FUNCTION_DIVIDE: {
		sprintf(name,"%s","gpr_divide");
		return;
	}
	case GPR_FUNCTION_DIVIDE_COMPLEX: {
		sprintf(name,"%s","gpr_divide_complex");
		return;
	}
	case GPR_FUNCTION_MODULUS: {
		sprintf(name,"%s","gpr_modulus");
		return;
	}
	case GPR_FUNCTION_FLOOR: {
		sprintf(name,"%s","gpr_floor");
		return;
	}
	case GPR_FUNCTION_GREATER_THAN: {
		sprintf(name,"%s","gpr_greater_than");
		return;
	}
	case GPR_FUNCTION_LESS_THAN: {
		sprintf(name,"%s","gpr_less_than");
		return;
	}
	case GPR_FUNCTION_EQUALS: {
		sprintf(name,"%s","gpr_equals");
		return;
	}
	case GPR_FUNCTION_AND: { sprintf(name,"%s","gpr_and"); return; }
	case GPR_FUNCTION_OR: { sprintf(name,"%s","gpr_or"); return; }
	case GPR_FUNCTION_XOR: { sprintf(name,"%s","gpr_xor"); return; }
	case GPR_FUNCTION_NOT: { sprintf(name,"%s","gpr_not"); return; }
	case GPR_FUNCTION_SET: { sprintf(name,"%s","gpr_store"); return; }
	case GPR_FUNCTION_GET: { sprintf(name,"%s","gpr_fetch"); return; }
	case GPR_FUNCTION_HEBBIAN: {
		sprintf(name,"%s","gpr_hebbian");
		return;
	}
	case GPR_FUNCTION_DEFUN: { sprintf(name,"%s","gpr_noop"); return; }
	case GPR_FUNCTION_MAIN: { sprintf(name,"%s","gpr_noop"); return; }
	case GPR_FUNCTION_ADF: {
		sprintf(name,"gpr_adf%d_",
				((int)f->value)%GPR_MAX_ARGUMENTS);
		return;
	}
	case GPR_FUNCTION_ARG: {
		sprintf(name,
				"temp_arg[call_depth][%d]",
				abs((int)f->value));
		return;
	}
	case GPR_FUNCTION_PROGRAM: {
		sprintf(name,"%s","gpr_noop");
		return;
	}
	default: {
		fprintf(stderr,"%s","\nERROR: Unknown function number %d\n");
	}
	}
}

static void gpr_dot_label(gpr_function * f, FILE * fp, int * ctr)
{
	int i = 0;
	char name[16];

	*ctr = *ctr + 1;

	if (f->function_type == GPR_FUNCTION_NONE) return;

	gpr_get_function_name(f->function_type,f->value,0,name);
	fprintf(fp,"  a%d [label=\"%s\"];\n", *ctr-1, name);

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			gpr_dot_label((gpr_function*)f->argv[i], fp, ctr);
		}
	}
}

static void gpr_dot_links(gpr_function * f, FILE * fp, int * ctr)
{
	int i = 0;
	int parent = *ctr;

	*ctr = *ctr + 1;

	if (f->function_type == GPR_FUNCTION_NONE) return;

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			if (f->argv[i]->function_type!=GPR_FUNCTION_NONE) {
				fprintf(fp,"  a%d -> a%d;\n", parent, *ctr);
			}
			gpr_dot_links((gpr_function*)f->argv[i], fp, ctr);
		}
	}
}

void gpr_dot(gpr_function * f, FILE * fp)
{
	int ctr;

	fprintf(fp,"%s","digraph graphname {\n");
	ctr = 0;
	gpr_dot_label(f, fp, &ctr);
	ctr = 0;
	gpr_dot_links(f, fp, &ctr);
	fprintf(fp,"%s","}\n");
}

/* save function parameters to file */
static void gpr_save_nodes(gpr_function * f, FILE * fp, int * ctr)
{
	int i = 0;

	*ctr = *ctr + 1;

	if (f->function_type == GPR_FUNCTION_NONE) return;

	fprintf(fp,"N %d %d %d %.8f\n",
			*ctr-1, f->function_type,
			f->argc, f->value);

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			gpr_save_nodes((gpr_function*)f->argv[i], fp, ctr);
		}
	}
}

/* save links to file */
static void gpr_save_links(gpr_function * f, FILE * fp, int * ctr)
{
	int i = 0;
	int parent = *ctr;

	*ctr = *ctr + 1;

	if (f->function_type == GPR_FUNCTION_NONE) return;

	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			if (f->argv[i]->function_type!=GPR_FUNCTION_NONE) {
				fprintf(fp,"%d %d %d\n", parent, i, *ctr);
			}
			gpr_save_links((gpr_function*)f->argv[i], fp, ctr);
		}
	}
}

/* save the given program to file */
void gpr_save(gpr_function * f, FILE * fp)
{
	int ctr = 0;

	gpr_save_nodes(f, fp, &ctr);
	ctr = 0;
	gpr_save_links(f, fp, &ctr);
}

/* loads the program state from file */
static void gpr_load_state(gpr_state * state,
						   int sensors, int actuators,
						   FILE * fp)
{
	int i;
	char line[256];

	/* load the sensor sources */
	if (fgets(line, 255, fp) != NULL ) {
		state->no_of_sensor_sources = atoi(line);
		if (state->no_of_sensor_sources>0) {
			state->sensor_source = (int*)malloc(sensors*sizeof(int));
			for (i = 0; i < sensors; i++) {
				if (fgets(line, 255, fp) != NULL ) {
					state->sensor_source[i] = atoi(line);
				}
			}
		}
	}
	else {
		printf("\nUnable to load number of software sources\n");
	}

	/* load the actuator destinations */
	if (fgets(line, 255, fp) != NULL ) {
		state->no_of_actuator_destinations = atoi(line);
		if (state->no_of_actuator_destinations > 0) {
			state->actuator_destination =
				(int*)malloc(actuators*sizeof(int));
			for (i = 0; i < actuators; i++) {
				if (fgets(line, 255, fp) != NULL ) {
					state->actuator_destination[i] = atoi(line);
				}
			}
		}
	}
	else {
		printf("\nUnable to load number of actuator destinations\n");
	}
}

/* saves the program state to file */
static void gpr_save_state(gpr_state * state, FILE * fp)
{
	int i;

	/* save the sensor sources */
	fprintf(fp,"%d\n",state->no_of_sensor_sources);
	if (state->no_of_sensor_sources > 0) {
		for (i = 0; i < state->no_of_sensors; i++) {
			fprintf(fp,"%d\n",state->sensor_source[i]);
		}
	}

	/* save the actuator_destinations */
	fprintf(fp,"%d\n",state->no_of_actuator_destinations);
	if (state->no_of_actuator_destinations > 0) {
		for (i = 0; i < state->no_of_actuators; i++) {
			fprintf(fp,"%d\n",state->actuator_destination[i]);
		}
	}
}

/* save the population to file */
void gpr_save_population(gpr_population *population, FILE * fp)
{
	int i;

	/* save population parameters */
	fprintf(fp,"%d\n",population->size);
	fprintf(fp,"%d\n",(&population->state[0])->no_of_registers);
	fprintf(fp,"%d\n",(&population->state[0])->no_of_sensors);
	fprintf(fp,"%d\n",(&population->state[0])->no_of_actuators);
	fprintf(fp,"%d\n",population->history.index);
	fprintf(fp,"%d\n",population->history.tick);
	fprintf(fp,"%d\n",population->history.interval);
	if ((&population->state[0])->ADF[0]==0) {
		fprintf(fp,"%d\n",0);
	}
	else {
		fprintf(fp,"%d\n",1);
	}
	for (i = 0; i < population->history.index; i++) {
		fprintf(fp,"%.10f\n",population->history.log[i]);
	}
	for (i = 0; i < population->size; i++) {
		/* save an individual */
		gpr_save((gpr_function*)&population->individual[i],fp);
		/* mark the end of the program */
		fprintf(fp,"%s",".\n");
		gpr_save_state((gpr_state*)&population->state[i],fp);
	}
}

/* save the system to file */
void gpr_save_system(gpr_system *system, FILE * fp)
{
	int i;

	/* save population parameters */
	fprintf(fp,"%d\n",system->size);
	fprintf(fp,"%d\n",system->migration_tick);
	for (i = 0; i < system->size; i++) {
		/* save the population */
		gpr_save_population((gpr_population*)&system->island[i],fp);
	}
}

/* load a population from file */
int gpr_load_population(gpr_population * population, FILE * fp)
{
	char line[256];
	int ctr=0,size=0,registers=0,sensors=0,actuators=0,index=0;
	int retval=0;
	int instruction_set[64], no_of_instructions=0;
	unsigned int random_seed = 123;
	int history_index=0,history_interval=0,history_tick=0;
	int integers_only=0;
	int ADFs=0;

	retval = GPR_LOAD_OK;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				/* population size */
				if (ctr==0) {
					size = atoi(line);
				}
				/* number of registers */
				if (ctr==1) {
					registers = atoi(line);
				}
				/* number of sensors */
				if (ctr==2) {
					sensors = atoi(line);
				}
				/* number of actuators */
				if (ctr==3) {
					actuators = atoi(line);
				}
				/* index in the fitness history */
				if (ctr==4) {
					history_index = atoi(line);
				}
				/* tick in the fitness history */
				if (ctr==5) {
					history_tick = atoi(line);
				}
				/* interval in the fitness history */
				if (ctr==6) {
					history_interval = atoi(line);
				}
				/* whether to use ADFs */
				if (ctr==7) {
					ADFs = atoi(line);
					break;
				}
				ctr++;
			}
		}
	}

	if (ctr==7) {
		/* Create an instruction set
		   It actually doesn't matter what this is, since it's
 		   only used by the init function and the randomly
		   generated individuals are then overwritten */
		no_of_instructions =
			gpr_default_instruction_set((int*)instruction_set);

		/* create the population */
		gpr_init_population(population,size,
							registers,sensors,actuators,
							5,-1,1, integers_only,ADFs,
							&random_seed,
							(int*)instruction_set, no_of_instructions);

		/* load the fitness history */
		population->history.index = history_index;
		population->history.tick = history_tick;
		population->history.interval = history_interval;
		for (index = 0; index < history_index; index++) {
			if (fgets(line , 255 , fp) != NULL ) {
				if (strlen(line) > 0) {
					population->history.log[index] = atof(line);
				}
			}
		}

		for (index = 0; index < size; index++) {
			/* clear the existing individual */
			gpr_free(&population->individual[index]);
			/* load the program */
			retval = gpr_load(&population->individual[index],fp);
			if (retval != GPR_LOAD_OK) {
				printf("Individual: %d\n", index);
				break;
			}
			/* load the program state */
			gpr_load_state(&population->state[index],
						   sensors, actuators, fp);

			if (ADFs>0) {
				/* enforce ADF structure */
				gpr_enforce_ADFs(&population->individual[index],
								 &population->state[index]);
			}
		}
	}

	if ((retval==GPR_LOAD_OK) && (index != size)) {
		retval = GPR_LOAD_POPULATION_SIZE;
	}

	return retval;
}

/* save the environment to file */
void gpr_save_environment(gpr_environment *population, FILE * fp)
{
	int i;

	/* save environment parameters */
	fprintf(fp,"%d\n",population->max_population_size);
	fprintf(fp,"%d\n",population->population_size);
	fprintf(fp,"%d\n",(&population->state[0])->no_of_registers);
	fprintf(fp,"%d\n",(&population->state[0])->no_of_sensors);
	fprintf(fp,"%d\n",(&population->state[0])->no_of_actuators);
	if ((&population->state[0])->ADF[0]==0) {
		fprintf(fp,"%d\n",0);
	}
	else {
		fprintf(fp,"%d\n",1);
	}
	for (i = 0; i < population->population_size; i++) {
		/* save an individual */
		gpr_save((gpr_function*)&population->individual[i],fp);
		/* mark the end of the program */
		fprintf(fp,"%s",".\n");
		gpr_save_state((gpr_state*)&population->state[i],fp);
	}
}

/* load an environment from file */
int gpr_load_environment(gpr_environment * population, FILE * fp)
{
	char line[256];
	int ctr=0,max_population_size=0,population_size=0;
	int registers=0,sensors=0,actuators=0,index=0;
	int retval=0;
	int instruction_set[64], no_of_instructions=0;
	unsigned int random_seed = 123;
	int integers_only=0;
	int ADFs=0;

	retval = GPR_LOAD_OK;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				/* maximum population size */
				if (ctr==0) {
					max_population_size = atoi(line);
				}
				/* population size */
				if (ctr==1) {
					population_size = atoi(line);
				}
				/* number of registers */
				if (ctr==2) {
					registers = atoi(line);
				}
				/* number of sensors */
				if (ctr==3) {
					sensors = atoi(line);
				}
				/* number of actuators */
				if (ctr==4) {
					actuators = atoi(line);
				}
				/* whether to use ADFs */
				if (ctr==5) {
					ADFs = atoi(line);
					break;
				}
				ctr++;
			}
		}
	}

	if (ctr==5) {
		/* Create an instruction set
		   It actually doesn't matter what this is, since it's
 		   only used by the init function and the randomly
		   generated individuals are then overwritten */
		no_of_instructions =
			gpr_default_instruction_set((int*)instruction_set);

		/* create the population */
		gpr_init_environment(population,
							 max_population_size,
							 population_size,
							 registers,sensors,actuators,
							 5,-1,1, integers_only,ADFs,
							 &random_seed,
							 (int*)instruction_set,
							 no_of_instructions);

		for (index = 0; index < population_size; index++) {
			/* clear the existing individual */
			gpr_free(&population->individual[index]);
			/* load the program */
			retval = gpr_load(&population->individual[index],fp);
			if (retval != GPR_LOAD_OK) {
				printf("Individual: %d\n", index);
				break;
			}
			/* load the program state */
			gpr_load_state(&population->state[index],
						   sensors, actuators, fp);

			if (ADFs>0) {
				/* enforce ADF structure */
				gpr_enforce_ADFs(&population->individual[index],
								 &population->state[index]);
			}
		}
	}

	if ((retval==GPR_LOAD_OK) && (index != population_size)) {
		retval = GPR_LOAD_POPULATION_SIZE;
	}

	return retval;
}

/* load a system from file */
void gpr_load_system(gpr_system * system,
					 FILE * fp,
					 int * instruction_set, int no_of_instructions)
{
	char line[256];
	int i,ctr=0,islands=0,population_per_island=10;
	int migration_tick=0;
	int registers=4;
	int sensors=1, actuators=1;
	int max_tree_depth=6;
	float min_value=-10, max_value=10;
	unsigned int random_seed = 123;
	int integers_only=0;
	int ADFs=0;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				switch(ctr) {
				case 0: {
					islands = atoi(line);
					break;
				}
				case 1: {
					migration_tick = atoi(line);
					break;
				}
				}
				if (ctr==1) break;
				ctr++;
			}
		}
	}

	if (islands==0) return;

	/* create a system.
	   It doesn't matter what the parameters are here, because
	   they will be overwritten later */
	gpr_init_system(system,
					islands,
					population_per_island,
					registers,
					sensors, actuators,
					max_tree_depth,
					min_value, max_value,
					integers_only,
					ADFs,
					&random_seed,
					instruction_set, no_of_instructions);

	/* set the current tick in the migration cycle */
	system->migration_tick = migration_tick;

	/* load each population */
	for (i = 0; i < islands; i++) {
		/* deallocate the existing population */
		gpr_free_population(&system->island[i]);
		/* load the population */
		gpr_load_population(&system->island[i], fp);
	}
}

/* load a program from file */
int gpr_load(gpr_function *f, FILE * fp)
{
	char line[256], param[128];
	int i, ctr, parent, argindex, child, result, no_of_nodes;
	gpr_function * fn[GPR_MAX_NODES];
	gpr_function * fn_parent, * fn_child;
	int node_number[GPR_MAX_NODES];
	const char separator = ' ';

	result = GPR_LOAD_OK;
	no_of_nodes = 0;

	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line) > 0) {
				if (line[0]=='.') break; /* end of program */
				if (line[0]=='N') {
					/* load node */
					if (no_of_nodes==0) {
						fn[no_of_nodes] = f;
					}
					else {
						fn[no_of_nodes] =
							(gpr_function*)malloc(sizeof(gpr_function));
					}
#ifdef DEBUG
					assert(fn[no_of_nodes] != 0);
#endif

					/* clear arguments */
					for (i = 0; i < GPR_MAX_ARGUMENTS; i++) {
						fn[no_of_nodes]->argv[i] = 0;
					}

					i = 2;

					/* get the node number */
					ctr = 0;
					while ((i < strlen(line)) && (ctr < 127)) {
						if (line[i] == separator) break;
						param[ctr++] = line[i++];
					}
					param[ctr] = 0;
					if (strlen(param)==0) {
						result = GPR_LOAD_NODE_NOT_FOUND;
						break;
					}
					node_number[no_of_nodes] = atoi(param);

					/* get the function type */
					ctr = 0;
					i++;
					while ((i < strlen(line)) && (ctr < 127)) {
						if (line[i] == separator) break;
						param[ctr++] = line[i++];
					}
					param[ctr] = 0;
					if (strlen(param)==0) {
						result = GPR_LOAD_FUNCTION_TYPE_NOT_FOUND;
						break;
					}
					fn[no_of_nodes]->function_type = atoi(param);

					/* get the number of arguments */
					ctr = 0;
					i++;
					while ((i < strlen(line)) && (ctr < 127)) {
						if (line[i] == separator) break;
						param[ctr++] = line[i++];
					}
					param[ctr] = 0;
					if (strlen(param)==0) {
						result = GPR_LOAD_ARGC_NOT_FOUND;
						break;
					}
					fn[no_of_nodes]->argc = atoi(param);

					/* initialise terminals */
					if (fn[no_of_nodes]->function_type ==
						GPR_FUNCTION_VALUE) {
						gpr_init(fn[no_of_nodes]);
					}

					if (fn[no_of_nodes]->function_type ==
						GPR_FUNCTION_ARG) {
						gpr_init(fn[no_of_nodes]);
						fn[no_of_nodes]->function_type =
							GPR_FUNCTION_ARG;
					}

					/* get the function value */
					ctr = 0;
					i++;
					while ((i < strlen(line)) && (ctr < 127)) {
						if (line[i] == separator) break;
						param[ctr++] = line[i++];
					}
					param[ctr] = 0;
					if (strlen(param)==0) {
						result = GPR_LOAD_VALUE_NOT_FOUND;
						break;
					}
					fn[no_of_nodes]->value = atof(param);

					/* increment the node count */
					no_of_nodes++;

					/* check for too many nodes */
					if (no_of_nodes >= GPR_MAX_NODES) {
						/* free allocated memory */
						for (i = 1; i < no_of_nodes; i++) {
							free(fn[i]);
						}
						result = GPR_LOAD_MAX_NODES_REACHED;
						break;
					}
				}
				else {
					/* load links between nodes */
					i = 0;

					/* get the parent node number */
					ctr = 0;
					while ((i < strlen(line)) && (ctr < 127)) {
						if (line[i] == separator) break;
						param[ctr++] = line[i++];
					}
					param[ctr] = 0;
					if (strlen(param)==0) {
						result = GPR_LOAD_PARENT_NOT_FOUND;
						break;
					}
					parent = atoi(param);

					/* get the argument index */
					ctr = 0;
					i++;
					while ((i < strlen(line)) && (ctr < 127)) {
						if (line[i] == separator) break;
						param[ctr++] = line[i++];
					}
					param[ctr] = 0;
					if (strlen(param)==0) {
						result = GPR_LOAD_ARG_NOT_FOUND;
						break;
					}
					argindex = atoi(param);

					/* get the child node number */
					ctr = 0;
					i++;
					while ((i < strlen(line)) && (ctr < 127)) {
						if (line[i] == separator) break;
						param[ctr++] = line[i++];
					}
					param[ctr] = 0;
					if (strlen(param)==0) {
						result = GPR_LOAD_CHILD_NOT_FOUND;
						break;
					}
					child = atoi(param);

					fn_parent = 0;
					fn_child = 0;

					/* get the function pointers for parent
					   and child nodes */
					for (i = 0; i < no_of_nodes; i++) {
						if (node_number[i] == parent) fn_parent = fn[i];
						if (node_number[i] == child) fn_child = fn[i];
					}

					/* connect the child to the parent */
					if ((fn_parent != 0) && (fn_child != 0)) {
						fn_parent->argv[argindex] = fn_child;
					}
					else {
						result = GPR_LOAD_LINK_NOT_FOUND;
					}
					
				}
			}
		}
	}
	return result;
}

/* saves the given program as an S-expression */
void gpr_S_expression(gpr_function * f, FILE * fp)
{
	int i = 0;
	char name[16];

	if (f->function_type == GPR_FUNCTION_NONE) return;

	gpr_get_function_name(f->function_type,f->value,0,name);

	fprintf(fp,"(%s", name);

	if (f->function_type != GPR_FUNCTION_VALUE) {
		fprintf(fp, "%s", " ");
	}

	/* add */
	if ((f->function_type==GPR_FUNCTION_EXP) ||
		(f->function_type==GPR_FUNCTION_ABS) ||
		(f->function_type==GPR_FUNCTION_ABS_COMPLEX) ||
		(f->function_type==GPR_FUNCTION_COSINE) ||
		(f->function_type==GPR_FUNCTION_ARCCOSINE) ||
		(f->function_type==GPR_FUNCTION_SINE) ||
		(f->function_type==GPR_FUNCTION_ARCSINE) ||
		(f->function_type==GPR_FUNCTION_SQUARE_ROOT) ||
		(f->function_type==GPR_FUNCTION_SQUARE_ROOT_COMPLEX)) {
		fprintf(fp,"%s","(+ ");
	}


	for (i = 0; i < f->argc; i++) {
		if (f->argv[i]!=0) {
			gpr_S_expression((gpr_function*)f->argv[i], fp);
		}
	}
	fprintf(fp,"%s", ") ");
	if ((f->function_type==GPR_FUNCTION_EXP) ||
		(f->function_type==GPR_FUNCTION_ABS) ||
		(f->function_type==GPR_FUNCTION_ABS_COMPLEX) ||
		(f->function_type==GPR_FUNCTION_COSINE) ||
		(f->function_type==GPR_FUNCTION_ARCCOSINE) ||
		(f->function_type==GPR_FUNCTION_SINE) ||
		(f->function_type==GPR_FUNCTION_ARCSINE) ||
		(f->function_type==GPR_FUNCTION_SQUARE_ROOT) ||
		(f->function_type==GPR_FUNCTION_SQUARE_ROOT_COMPLEX)) {
		fprintf(fp,"%s",") ");
	}
}

/* returns a C program */
void gpr_c(gpr_function * f, FILE * fp)
{
	int i = 0;
	char name[64];

	get_c_function_name(f,name);

	if (is_terminal(f->function_type)==0) {
		fprintf(fp, "%s%d(", name, f->argc);
		for (i = 0; i < f->argc; i++) {
			if (f->argv[i]!=0) {
				gpr_c((gpr_function*)f->argv[i], fp);
				if (i < f->argc-1) {
					fprintf(fp,"%s",",");
				}
				else {
					if (f->function_type==GPR_FUNCTION_ADF) {
						fprintf(fp,"%s",", call_depth");
					}
				}
			}
		}
		fprintf(fp,"%s", ")");
	}
	else {
		fprintf(fp, "%s", name);
	}
}

static void gpr_arduino_get_inputs(FILE * fp,
								   int no_of_digital_inputs,
								   int no_of_analog_inputs,
								   int digital_high)
{
	fprintf(fp,"%s","void get_inputs() {\n");

	if (no_of_digital_inputs>0) {
		fprintf(fp,"%s","  for (int i = 0; i < ");
		fprintf(fp,"%s","digital_inputs; i++) {\n");
		fprintf(fp,"%s","    if (dInput[i] >= 0) {\n");
		fprintf(fp,"%s","      if (digitalRead(dInput[i]) == HIGH) {\n");
		fprintf(fp,"        sensor[i] = %d;\n",digital_high);
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","      else {\n");
		fprintf(fp,"%s","        sensor[i] = 0;\n");
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","    else {\n");
		fprintf(fp,"%s","      // Input the current time step\n");
		fprintf(fp,"%s","      sensor[i] = tick;\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	if (no_of_analog_inputs>0) {
		fprintf(fp,"%s","  for (int i = 0; i < analog_inputs; i++) {\n");
		fprintf(fp,"%s","    if (aInput[i] >= 0) {\n");
		if (no_of_digital_inputs==0) {
			fprintf(fp,"%s","      sensor[i] = ");
			fprintf(fp,"%s","analogRead(aInput[i]);\n");
		}
		else {
			fprintf(fp,"%s",
					"      sensor[i+digital_inputs] = ");
			fprintf(fp,"%s","analogRead(aInput[i]);\n");
		}
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","    else {\n");
		fprintf(fp,"%s","      // Input the current time step\n");
		if (no_of_digital_inputs==0) {
			fprintf(fp,"%s","      sensor[i] = tick;\n");
		}
		else {
			fprintf(fp,"%s","      sensor[i+digital_inputs] = tick;\n");
		}
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	fprintf(fp,"%s","}\n\n");
}

static void gpr_c_get_inputs(FILE * fp)
{
	fprintf(fp,"%s",
			"static void get_inputs(int argc, char *argv[])\n{\n");
	fprintf(fp,"%s","  int i;\n\n");

	fprintf(fp,"%s","  for (i = 0; i < sensors; i++) {\n");
	fprintf(fp,"%s","    sensor[i] = atof(argv[i]);\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","}\n\n");
}

static void gpr_arduino_set_outputs(FILE * fp,
									int no_of_digital_outputs,
									int no_of_analog_outputs)
{
	fprintf(fp,"%s","void set_outputs() {\n");

	if (no_of_digital_outputs>0) {
		fprintf(fp,"%s",
				"  for (int i = 0; i < digital_outputs; i++) {\n");
		fprintf(fp,"%s","    if (dOutput[i] >= 0) {\n");
		fprintf(fp,"%s","      if (actuator[i] > 0) {\n");
		fprintf(fp,"%s","        digitalWrite(dOutput[i],HIGH);\n");
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","      else {\n");
		fprintf(fp,"%s","        digitalWrite(dOutput[i],LOW);\n");
		fprintf(fp,"%s","      }\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	if (no_of_analog_outputs>0) {
		fprintf(fp,"%s",
				"  for (int i = 0; i < analog_outputs; i++) {\n");
		fprintf(fp,"%s","    if (aOutput[i] >= 0) {\n");
		if (no_of_digital_outputs==0) {
			fprintf(fp,"%s",
					"      analogWrite(aOutput[i],");
			fprintf(fp,"%s","abs((int)actuator[i])%256);\n");
		}
		else {
			fprintf(fp,"%s",
					"      analogWrite(aOutput[i],");
			fprintf(fp,"%s",
					"abs((int)actuator[i+digital_outputs])%256);\n");
		}
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n");
	}

	fprintf(fp,"%s","}\n\n");
}

static void gpr_c_set_outputs(FILE * fp)
{
	fprintf(fp,"%s","static void set_outputs() {\n");
	fprintf(fp,"%s","  int i;\n\n");

	fprintf(fp,"%s","  for (i = 0; i < actuators; i++) {\n");
	fprintf(fp,"%s","    if (i > 0) printf(\" \");\n");
	fprintf(fp,"%s","    printf(\"%.3f\",actuator[i]);\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  printf(\"\\n\");\n");
	fprintf(fp,"%s","}\n\n");
}

/* arduino setup */
static void gpr_arduino_setup(FILE * fp,
							  int no_of_digital_inputs,
							  int no_of_analog_inputs,
							  int no_of_digital_outputs,
							  int no_of_analog_outputs,							  
							  int baud_rate)
{
	fprintf(fp,"%s","void setup() {\n");

	fprintf(fp,"%s","  analogReference(DEFAULT);\n\n");

	fprintf(fp,"%s","  // initialize serial communication:\n");
	fprintf(fp,"  Serial.begin(%d);\n\n", baud_rate);

	fprintf(fp,"%s","  // Clear the sensor values\n");
	fprintf(fp,"%s","  for (int i = 0; i < sensors; i++) {\n");
	fprintf(fp,"%s","    sensor[i] = 0;\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  // Clear the actuator values\n");
	fprintf(fp,"%s","  for (int i = 0; i < actuators; i++) {\n");
	fprintf(fp,"%s","    actuator[i] = 0;\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  // Clear the register values\n");
	fprintf(fp,"%s","  for (int i = 0; i < registers; i++) {\n");
	fprintf(fp,"%s","    registerlist[i] = 0;\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  // Clear the temporary arguments used by ADFs\n");
	fprintf(fp,"  for (int i = 0; i < %d; i++) {\n",
			GPR_MAX_CALL_DEPTH);
	fprintf(fp,"    for (int j = 0; j < %d; j++) {\n",
			GPR_MAX_ARGUMENTS);
	fprintf(fp,"%s","      temp_arg[i][j] = 0;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n\n");

	if (no_of_digital_inputs>0) {
		fprintf(fp,"%s","  // Define digital inputs\n");
		fprintf(fp,"%s",
				"  for (int i = 0; i < digital_inputs; i++) {\n");
		fprintf(fp,"%s","    if (dInput[i]>=0) {\n");
		fprintf(fp,"%s","      pinMode(dInput[i], INPUT);\n");
		fprintf(fp,"%s","    }\n");
		fprintf(fp,"%s","  }\n\n");
	}
	if (no_of_digital_outputs>0) {
		fprintf(fp,"%s","  // Define digital outputs\n");
		fprintf(fp,"%s",
				"  for (int i = 0; i < digital_inputs; i++) {\n");
		fprintf(fp,"%s","    pinMode(dOutput[i], OUTPUT);\n");
		fprintf(fp,"%s","  }\n\n");
	}
	if (no_of_analog_outputs>0) {
		fprintf(fp,"%s","  // Define analog (PWM) outputs\n");
		fprintf(fp,"%s",
				"  for (int i = 0; i < analog_outputs; i++) {\n");
		fprintf(fp,"%s","    pinMode(aOutput[i], OUTPUT);\n");
		fprintf(fp,"%s","  }\n\n");
	}

	fprintf(fp,"%s","}\n\n");
}

static void gpr_c_setup(FILE * fp)
{
	fprintf(fp,"%s","static void setup()\n{\n");
	fprintf(fp,"%s","  int i, j;\n\n");

	fprintf(fp,"%s","  /* Clear the sensor values */\n");
	fprintf(fp,"%s","  for (i = 0; i < sensors; i++) {\n");
	fprintf(fp,"%s","    sensor[i] = 0;\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  /* Clear the actuator values */\n");
	fprintf(fp,"%s","  for (i = 0; i < actuators; i++) {\n");
	fprintf(fp,"%s","    actuator[i] = 0;\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  /* Clear the register values */\n");
	fprintf(fp,"%s","  for (i = 0; i < registers; i++) {\n");
	fprintf(fp,"%s","    registerlist[i] = 0;\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  /* Clear the temporary arguments ");
	fprintf(fp,"%s","used by ADFs */\n");
	fprintf(fp,"  for (i = 0; i < %d; i++) {\n",GPR_MAX_CALL_DEPTH);
	fprintf(fp,"    for (j = 0; j < %d; j++) {\n",GPR_MAX_ARGUMENTS);
	fprintf(fp,"%s","      temp_arg[i][j] = 0;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n\n");
	fprintf(fp,"%s","}\n\n");
}

/* save ADFs to a C program */
static void gpr_c_ADFs(gpr_function * f, FILE * fp, int argc)
{
	int i, j, active;
	char str[32];
	gpr_function * f2;

	for (i = 0; i < f->argc-1; i++) {
		f2 = f->argv[i];

		fprintf(fp,"/* Automatically Defined Function number %d\n",i);
		fprintf(fp,"   with %d arguments */\n",argc);

		sprintf(str,"gpr_adf%d_",i);
		gpr_function_c(fp, str, argc, ", int call_depth");

		active = 0;
		if (f2->argv[0]!=0) {
			if (f2->argv[0]->function_type!=GPR_FUNCTION_NONE) {

				fprintf(fp,"  if (call_depth >= %d) return 0;\n",
						GPR_MAX_CALL_DEPTH-1);

				fprintf(fp,"%s","  /* Store arguments in a temporary array */\n");
				fprintf(fp,"%s","  call_depth++;\n");
				for (j = 0; j < argc; j++) {
					fprintf(fp,"  temp_arg[call_depth][%d] = v%d;\n",
							j,j+1);
				}

				fprintf(fp,"%s","\n  return ");    
				gpr_c(f2, fp);
				fprintf(fp,"%s",";\n");
				active = 1;
			}
		}

		if (active == 0) {
			fprintf(fp,"%s","  /* This function does nothing */\n");
			fprintf(fp,"%s","  return 0;\n");
		}

		fprintf(fp,"%s","}\n\n");
	}
}

/* arduino main loop */
static void gpr_arduino_main(gpr_function * f,FILE * fp, int ADFs)
{
	int argc;

	/* save ADFs */
	if (ADFs>0) {
		for (argc = 2; argc < GPR_MAX_ARGUMENTS; argc++) {
			gpr_c_ADFs(f, fp, argc);
		}
	}

	fprintf(fp,"%s","void loop() {\n");
	fprintf(fp,"%s","  get_inputs();\n  "); 
   
	if (ADFs<=0) {
		/* run the entire tree */
		gpr_c(f, fp);
	}
	else {
		/* ADF mode - run the main program */
		fprintf(fp,"%s","call_depth=0;\n  "); 
		gpr_c(f->argv[f->argc-1], fp);
	}
	fprintf(fp,"%s",";\n");

	fprintf(fp,"%s","  set_outputs();\n\n");
	fprintf(fp,"%s","  // Increment the time step\n");
	fprintf(fp,"%s","  tick++;\n");
	fprintf(fp,"%s","  if (tick>32000) tick=0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* read arguments from stdin */
static void gpr_c_stdin_args(FILE * fp)
{
	fprintf(fp,"%s","static int read_stdin_args(char *argv[])\n{\n");
	fprintf(fp,"%s","  char argstr[1024], value[1024], *retval;\n");
	fprintf(fp,"%s","  int i,ctr=0,field_index=0;\n");
	fprintf(fp,"%s","  FILE * fp = fopen(\"stdin\",\"r\");\n\n");

	fprintf(fp,"%s","  if (!fp) return 0;\n");
	fprintf(fp,"%s","  while (!feof(fp)) {\n");
	fprintf(fp,"%s","    retval = fgets(argstr,1023,fp);\n");
	fprintf(fp,"%s","    if (retval) {\n");
	fprintf(fp,"%s","      for (i = 0; i < strlen(argstr); i++) {\n");
	fprintf(fp,"%s","        if ((argstr[i]==' ') ||\n");
	fprintf(fp,"%s","            (argstr[i]==',') ||\n");
	fprintf(fp,"%s","            (argstr[i]==';') ||\n");
	fprintf(fp,"%s","            (i == strlen(argstr)-1)) {\n");
	fprintf(fp,"%s","          if (i == strlen(argstr)-1) {\n");
	fprintf(fp,"%s","            value[ctr++] = argstr[i];\n");
	fprintf(fp,"%s","          }\n");
	fprintf(fp,"%s","          value[ctr] = 0;\n");
	fprintf(fp,"%s","          ctr = 0;\n");
	fprintf(fp,"%s","          argv[field_index] = ");
	fprintf(fp,"%s","(char*)malloc(strlen(value)+2);\n");
	fprintf(fp,"%s","          sprintf(argv[field_index++],");
	fprintf(fp,"%s","\"%s\",value);\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","        else {\n");
	fprintf(fp,"%s","          value[ctr++] = argstr[i];\n");
	fprintf(fp,"%s","        }\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  fclose(fp);\n");
	fprintf(fp,"%s","  return field_index;\n");
	fprintf(fp,"%s","}\n\n");
}

/* Call a remote classifier using the given sensor values
   and return the actuator values.
   This is fairly crude, and a C xmlrpc client would be possbile
   but the ruby script is only three lines of code
 */
int gpr_xmlrpc_client(char * service_name,
					  int port, char * hostname,
					  int sensors, int actuators,
					  float * sensor_value,
					  float * actuator_value)
{
	FILE * fp;
	char command_str[2048];
	char ruby_script_filename[256];
	char result_filename[256],line[256],valuestr[256];
	int i,ctr=0,field_index=0;

	/* temporary filenames */
	sprintf(ruby_script_filename,
			"%stemp_libgpr_xmlrpc.rb",
			GPR_TEMP_DIRECTORY);
	sprintf(result_filename,"%slibgpr_result.txt",GPR_TEMP_DIRECTORY);

	/* create a ruby script */
	fp = fopen(ruby_script_filename,"w");
	if (!fp) return -1;

	fprintf(fp,"%s","require 'xmlrpc/client'\n\n");
	fprintf(fp,"server = XMLRPC::Client.new2(\"%s:%d\")\n",
			hostname,port);
	fprintf(fp,"puts server.call(\"%s.classify\", ", service_name);
	for (i  =0; i < sensors; i++) {
		if (i>0) fprintf(fp,"%s",",");
		fprintf(fp,"%.5f",sensor_value[i]);
	}
	fprintf(fp,"%s",").inspect\n");   
	fclose(fp);

	/* run the ruby script */
	sprintf(command_str,"ruby %s > %s",
			ruby_script_filename,
			result_filename);
	if (system(command_str)!=0) return -2;

	/* open the results file */
	fp = fopen(result_filename,"r");
	if (!fp) return -3;

	/* read the results into the actuator_value array */
	while (!feof(fp)) {
		if (fgets(line , 255 , fp) != NULL ) {
			if (strlen(line)>0) {
				for (i = 0; i < strlen(line); i++) {
					if (line[i]==',') {
						valuestr[ctr]=0;
						actuator_value[field_index++] = atof(valuestr);
						ctr=0;
					}
					else {
						if (((line[i]>='0') && (line[i]<='9')) ||
							(line[i]=='.')) {
							valuestr[ctr++] = line[i];
						}
					}
				}
				if (ctr>0) {
					valuestr[ctr]=0;
					ctr=0;
					actuator_value[field_index++] = atof(valuestr);
				}
			}
		}
	}

	fclose(fp);

	/* delete the temporary files */
	sprintf(command_str,"rm %s",
			ruby_script_filename);
	i = system(command_str);
	sprintf(command_str,"rm %s",
			result_filename);
	i = system(command_str);

	return 0;
}

/* Saves a ruby script which implements an XMLRPC server
   so that an exported C program can be remotely called.
   Ruby is used here because it requires very few lines of code */
void gpr_xmlrpc_server(char * ruby_script_filename,
					   char * service_name, int port,
					   char * c_program_filename,
					   int sensors, int actuators)
{
	int i;
	FILE * fp;

	fp = fopen(ruby_script_filename,"w");
	if (!fp) return;

	fprintf(fp,"%s","#!/usr/bin/env ruby\n\n");
	fprintf(fp,"%s","require \"xmlrpc/server\"\n\n");

	fprintf(fp,"%s","class Libgpr_server\n");
	fprintf(fp,"  INTERFACE = XMLRPC::interface(\"%s\") {\n",
			service_name);

	fprintf(fp,"%s","    meth 'int sensors()', ");
	fprintf(fp,"%s","'Returns the number of sensors', 'sensors'\n");
	fprintf(fp,"%s","    meth 'int actuators()', ");
	fprintf(fp,"%s","'Returns the number of actuators', 'actuators'\n");

	if (actuators<2) {
		fprintf(fp,"%s","    meth 'float classify(");
	}
	else {
		fprintf(fp,"%s","    meth 'array classify(");
	}
	for (i = 0; i < sensors; i++) {
		if (i > 0) fprintf(fp,"%s",", ");
		fprintf(fp,"%s","float");
	}
	fprintf(fp,"%s",")', 'Runs the classifier', 'classify'\n");
	fprintf(fp,"%s","  }\n\n");

	fprintf(fp,"%s","  def sensors()\n");
	fprintf(fp,"    return %d\n",sensors);
	fprintf(fp,"%s","  end\n\n");

	fprintf(fp,"%s","  def actuators()\n");
	fprintf(fp,"    return %d\n",actuators);
	fprintf(fp,"%s","  end\n\n");

	fprintf(fp,"%s","  def classify(");
	for (i = 0; i < sensors; i++) {
		if (i > 0) fprintf(fp,"%s",", ");
		fprintf(fp,"v%d",i);
	}
	fprintf(fp,"%s",")\n");
	fprintf(fp,"    value = `%s ",c_program_filename);
	for (i = 0; i < sensors; i++) {
		if (i > 0) fprintf(fp,"%s"," ");
		fprintf(fp,"#{v%d}",i);
	}
	fprintf(fp,"%s","`\n");

	if (actuators<2) {
		fprintf(fp,"%s","    return value\n");
	}
	else {
		fprintf(fp,"%s","    return value.split(' ')\n");
	}
    fprintf(fp,"%s","  end\n");
	fprintf(fp,"%s","end\n\n");

	fprintf(fp,"server = XMLRPC::Server.new(%d)\n",port);
	fprintf(fp,"%s","server.add_handler(Libgpr_server::");
	fprintf(fp,"%s","INTERFACE, Libgpr_server.new)\n");
	fprintf(fp,"%s","server.serve\n\n");
	fclose(fp);
}

/* C program main loop */
static void gpr_c_main(gpr_function * f,FILE * fp, int ADFs)
{
	int argc;

	/* save ADFs */
	if (ADFs>0) {
		for (argc = 2; argc < GPR_MAX_ARGUMENTS; argc++) {
			gpr_c_ADFs(f, fp, argc);
		}
	}

	fprintf(fp,"%s","int main(int argc, char *argv[])\n{\n");
	fprintf(fp,"%s","  char * stdin_args[1024];\n");
	fprintf(fp,"%s","  int i,no_of_stdin_args = ");
	fprintf(fp,"%s","read_stdin_args(stdin_args);\n\n");
	fprintf(fp,"%s","  setup();\n\n"); 
	fprintf(fp,"%s","  if (no_of_stdin_args>0) {\n");
	fprintf(fp,"%s","    if (no_of_stdin_args != sensors) {\n");
	fprintf(fp,"%s","      printf(\"Invalid number of ");
	fprintf(fp,"%s","arguments %d/%d\\n\",no_of_stdin_args,sensors);\n");
	fprintf(fp,"%s","      for (i = 0; i < no_of_stdin_args; i++) {\n");
	fprintf(fp,"%s","        free(stdin_args[i]);\n");
	fprintf(fp,"%s","      }\n");
	fprintf(fp,"%s","      return -1;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","    get_inputs(no_of_stdin_args,stdin_args);\n");
	fprintf(fp,"%s","    for (i = 0; i < no_of_stdin_args; i++) {\n");
	fprintf(fp,"%s","      free(stdin_args[i]);\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","  }\n");
	fprintf(fp,"%s","  else {\n");
	fprintf(fp,"%s","    if (argc-1 != sensors) {\n");
	fprintf(fp,"%s","      printf(\"Invalid number of ");
	fprintf(fp,"%s","arguments %d/%d\\n\",argc-1,sensors);\n");
	fprintf(fp,"%s","      return -1;\n");
	fprintf(fp,"%s","    }\n");
	fprintf(fp,"%s","    get_inputs(argc-1,&argv[1]);\n");
	fprintf(fp,"%s","  }\n\n");
   
	if (ADFs<=0) {
		/* run the entire tree */
		fprintf(fp,"%s","  ");
		gpr_c(f, fp);
	}
	else {
		/* ADF mode - run the main program */
		fprintf(fp,"%s","  call_depth=0;\n  "); 
		gpr_c(f->argv[f->argc-1], fp);
	}
	fprintf(fp,"%s",";\n\n");

	fprintf(fp,"%s","  set_outputs();\n\n");
	fprintf(fp,"%s","  /* Increment the time step */\n");
	fprintf(fp,"%s","  tick++;\n");
	fprintf(fp,"%s","  if (tick>32000) tick=0;\n");
	fprintf(fp,"%s","  return 0;\n");
	fprintf(fp,"%s","}\n\n");
}

/* export the given individual as an Arduino program */
void gpr_arduino(gpr_function * f,
				 int baud_rate,
				 int digital_high,
				 int no_of_sensors, int no_of_actuators,
				 int no_of_registers,
				 int * digital_inputs, int no_of_digital_inputs,
				 int * analog_inputs, int no_of_analog_inputs,
				 int * digital_outputs, int no_of_digital_outputs,
				 int * analog_outputs, int no_of_analog_outputs,
				 int ADFs, FILE * fp)
{
	int i;

	/* check that the number of inputs matches the number of sensors */
	if (no_of_digital_inputs+no_of_analog_inputs!=
		no_of_sensors) {
		fprintf(stderr,"%s",
				"Number of digital and analog inputs " \
				"should equal the number of sensors\n");
		return;
	}

	/* check that the number of outputs matches
	   the number of actuators */
	if (no_of_digital_outputs+no_of_analog_outputs!=
		no_of_actuators) {
		fprintf(stderr,"%s",
				"Number of digital and analog outputs " \
				"should equal the number of actuators\n");
		return;
	}

	/* comment header */
	fprintf(fp,"%s","// Genetic Program\n");
	fprintf(fp,"%s","// Evolved using libgpr\n");
	fprintf(fp,"%s","// https://launchpad.net/libgpr\n\n");

	fprintf(fp,"const int sensors = %d;\n",no_of_sensors);
	fprintf(fp,"const int actuators = %d;\n",no_of_actuators);
	fprintf(fp,"const int registers = %d;\n\n",no_of_registers);

	if (ADFs > 0) {
		fprintf(fp,"%s","int call_depth;\n\n");
	}
	
	fprintf(fp,"float sensor[%d];\n",no_of_sensors);
	fprintf(fp,"float actuator[%d];\n",no_of_actuators);
	fprintf(fp,"float registerlist[%d];\n",no_of_registers);
	fprintf(fp,"float temp_arg[%d][%d];\n\n",
			GPR_MAX_CALL_DEPTH,GPR_MAX_ARGUMENTS);
	if (no_of_digital_inputs > 0) {
		fprintf(fp,"const int digital_inputs = %d;\n",
				no_of_digital_inputs);
	}
	if (no_of_analog_inputs > 0) {
		fprintf(fp,"const int analog_inputs = %d;\n",
				no_of_analog_inputs);
	}
	if (no_of_digital_outputs > 0) {
		fprintf(fp,"const int digital_outputs = %d;\n",
				no_of_digital_outputs);
	}
	if (no_of_analog_outputs > 0) {
		fprintf(fp,"const int analog_outputs = %d;\n",
				no_of_analog_outputs);
	}
	fprintf(fp,"%s","int tick = 0;\n");
	fprintf(fp,"%s","\n");

	/* digital inputs */
	if (no_of_digital_inputs>0) {
		fprintf(fp,"%s","const int dInput[] = {");
		for (i = 0; i < no_of_digital_inputs; i++) {
			fprintf(fp,"%d",digital_inputs[i]);
			if (i<no_of_digital_inputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}

	/* analog inputs */
	if (no_of_analog_inputs>0) {
		fprintf(fp,"%s","const int aInput[] = {");
		for (i = 0; i < no_of_analog_inputs; i++) {
			fprintf(fp,"A%d",analog_inputs[i]);
			if (i<no_of_analog_inputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}

	/* digital outputs */
	if (no_of_digital_outputs>0) {
		fprintf(fp,"%s","int dOutput[] = {");
		for (i = 0; i < no_of_digital_outputs; i++) {
			fprintf(fp,"%d",digital_outputs[i]);
			if (i<no_of_digital_outputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}

	/* analog outputs */
	if (no_of_analog_outputs>0) {
		fprintf(fp,"%s","int aOutput[] = {");
		for (i = 0; i < no_of_analog_outputs; i++) {
			fprintf(fp,"%d",analog_outputs[i]);
			if (i<no_of_analog_outputs-1) {
				fprintf(fp,"%s",",");
			}
		}
		fprintf(fp,"%s","};\n");
	}
	fprintf(fp,"%s","\n");

	/* function library */
	for (i = 2; i <= GPR_MAX_ARGUMENTS; i++) {
		if (gpr_contains_function(f, GPR_FUNCTION_CUSTOM, i)==1) {
			gpr_custom_function_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_AVERAGE, i)==1) {
			gpr_average_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ADD, i)==1) {
			gpr_add_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ADD_COMPLEX, i)==1) {
			gpr_add_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SUBTRACT, i)==1) {
			gpr_subtract_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SUBTRACT_COMPLEX, i)==1) {
			gpr_subtract_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MULTIPLY, i)==1) {
			gpr_multiply_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MULTIPLY_COMPLEX, i)==1) {
			gpr_multiply_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_DIVIDE, i)==1) {
			gpr_divide_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_DIVIDE_COMPLEX, i)==1) {
			gpr_divide_complex_c(fp, i);
		}
		if ((gpr_contains_function(f, GPR_FUNCTION_NOOP1, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_NOOP2, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_NOOP3, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_NOOP4, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_PROGRAM, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_MAIN, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_DEFUN, i)==1)) {
			gpr_noop_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ABS, i)==1) {
			gpr_abs_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ABS_COMPLEX, i)==1) {
			gpr_abs_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_POW, i)==1) {
			gpr_pow_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_EXP, i)==1) {
			gpr_exp_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SIGMOID, i)==1) {
			gpr_sigmoid_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MIN, i)==1) {
			gpr_min_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MAX, i)==1) {
			gpr_max_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_NEGATE, i)==1) {
			gpr_negate_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MODULUS, i)==1) {
			gpr_modulus_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_FLOOR, i)==1) {
			gpr_floor_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SQUARE_ROOT, i)==1) {
			gpr_square_root_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SQUARE_ROOT_COMPLEX, i)==1) {
			gpr_square_root_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SINE, i)==1) {
			gpr_sine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ARCSINE, i)==1) {
			gpr_arcsine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_COSINE, i)==1) {
			gpr_cosine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ARCCOSINE, i)==1) {
			gpr_arccosine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_GREATER_THAN, i)==1) {
			gpr_greater_than_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_LESS_THAN, i)==1) {
			gpr_less_than_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_EQUALS, i)==1) {
			gpr_equals_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_AND, i)==1) {
			gpr_and_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_OR, i)==1) {
			gpr_or_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_XOR, i)==1) {
			gpr_xor_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_NOT, i)==1) {
			gpr_not_c(fp, i);
		}
	}

	if (gpr_contains_function(f, GPR_FUNCTION_WEIGHT, 1)==1) {
		gpr_weight_c(fp, 1, f);
	}
	if (gpr_contains_function(f, GPR_FUNCTION_SET, 2)==1) {
		gpr_store_c(fp,2);
	}
	if (gpr_contains_function(f, GPR_FUNCTION_GET, 2)==1) {
		gpr_fetch_c(fp,2);
	}

	/* inputs and outputs */
	gpr_arduino_get_inputs(fp, no_of_digital_inputs,
						   no_of_analog_inputs,
						   digital_high);
	gpr_arduino_set_outputs(fp, no_of_digital_outputs,
							no_of_analog_outputs);

	gpr_arduino_setup(fp,
					  no_of_digital_inputs,
					  no_of_analog_inputs,
					  no_of_digital_outputs,
					  no_of_analog_outputs,							  
					  baud_rate);
	gpr_arduino_main(f, fp, ADFs);
}

/* export the given individual as a C program */
void gpr_c_program(gpr_function * f,
				   int no_of_sensors, int no_of_actuators,
				   int no_of_registers,
				   int ADFs, FILE * fp)
{
	int i;

	/* comment header */
	fprintf(fp,"%s","/* Genetic Program\n");
	fprintf(fp,"%s","   Evolved using libgpr\n");
	fprintf(fp,"%s","   https://launchpad.net/libgpr\n\n");

	fprintf(fp,"%s","   To compile:\n     ");
	fprintf(fp,"%s","gcc -Wall -std=c99 -pedantic ");
	fprintf(fp,"%s","-o agent agent.c -lm\n*/\n\n");

	fprintf(fp,"%s","#include<stdio.h>\n");
	fprintf(fp,"%s","#include<stdlib.h>\n");
	fprintf(fp,"%s","#include<string.h>\n");
	fprintf(fp,"%s","#include<math.h>\n\n");

	fprintf(fp,"const int sensors = %d;\n",no_of_sensors);
	fprintf(fp,"const int actuators = %d;\n",no_of_actuators);
	fprintf(fp,"const int registers = %d;\n\n",no_of_registers);

	if (ADFs > 0) {
		fprintf(fp,"%s","int call_depth;\n\n");
	}
	
	fprintf(fp,"float sensor[%d];\n",no_of_sensors);
	fprintf(fp,"float actuator[%d];\n",no_of_actuators);
	fprintf(fp,"float registerlist[%d];\n",no_of_registers);
	fprintf(fp,"float temp_arg[%d][%d];\n\n",
			GPR_MAX_CALL_DEPTH,GPR_MAX_ARGUMENTS);
	fprintf(fp,"%s","int tick = 0;\n");
	fprintf(fp,"%s","\n");

	/* function library */
	for (i = 2; i <= GPR_MAX_ARGUMENTS; i++) {
		if (gpr_contains_function(f, GPR_FUNCTION_CUSTOM, i)==1) {
			gpr_custom_function_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_AVERAGE, i)==1) {
			gpr_average_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ADD, i)==1) {
			gpr_add_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ADD_COMPLEX, i)==1) {
			gpr_add_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SUBTRACT, i)==1) {
			gpr_subtract_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SUBTRACT_COMPLEX, i)==1) {
			gpr_subtract_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MULTIPLY, i)==1) {
			gpr_multiply_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MULTIPLY_COMPLEX, i)==1) {
			gpr_multiply_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_DIVIDE, i)==1) {
			gpr_divide_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_DIVIDE_COMPLEX, i)==1) {
			gpr_divide_complex_c(fp, i);
		}
		if ((gpr_contains_function(f, GPR_FUNCTION_NOOP1, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_NOOP2, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_NOOP3, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_NOOP4, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_PROGRAM, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_MAIN, i)==1) ||
			(gpr_contains_function(f, GPR_FUNCTION_DEFUN, i)==1)) {
			gpr_noop_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ABS, i)==1) {
			gpr_abs_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ABS_COMPLEX, i)==1) {
			gpr_abs_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_POW, i)==1) {
			gpr_pow_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_EXP, i)==1) {
			gpr_exp_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SIGMOID, i)==1) {
			gpr_sigmoid_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MIN, i)==1) {
			gpr_min_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MAX, i)==1) {
			gpr_max_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_NEGATE, i)==1) {
			gpr_negate_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_MODULUS, i)==1) {
			gpr_modulus_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_FLOOR, i)==1) {
			gpr_floor_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SQUARE_ROOT, i)==1) {
			gpr_square_root_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SQUARE_ROOT_COMPLEX, i)==1) {
			gpr_square_root_complex_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_SINE, i)==1) {
			gpr_sine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ARCSINE, i)==1) {
			gpr_arcsine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_COSINE, i)==1) {
			gpr_cosine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_ARCCOSINE, i)==1) {
			gpr_arccosine_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_GREATER_THAN, i)==1) {
			gpr_greater_than_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_LESS_THAN, i)==1) {
			gpr_less_than_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_EQUALS, i)==1) {
			gpr_equals_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_AND, i)==1) {
			gpr_and_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_OR, i)==1) {
			gpr_or_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_XOR, i)==1) {
			gpr_xor_c(fp, i);
		}
		if (gpr_contains_function(f, GPR_FUNCTION_NOT, i)==1) {
			gpr_not_c(fp, i);
		}
	}

	if (gpr_contains_function(f, GPR_FUNCTION_WEIGHT, 1)==1) {
		gpr_weight_c(fp, 1, f);
	}
	if (gpr_contains_function(f, GPR_FUNCTION_SET, 2)==1) {
		gpr_store_c(fp,2);
	}
	if (gpr_contains_function(f, GPR_FUNCTION_GET, 2)==1) {
		gpr_fetch_c(fp,2);
	}

	/* inputs and outputs */
	gpr_c_get_inputs(fp);
	gpr_c_set_outputs(fp);

	gpr_c_setup(fp);
	gpr_c_stdin_args(fp);
	gpr_c_main(f, fp, ADFs);
}

/* uses gnuplot to plot the fitness history for the given population */
int gpr_plot_history(gpr_population * population,
					 int history_type,
					 char * filename, char * title,
					 int image_width, int image_height)
{
	int index,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float value, min_fitness = 0;
	float max_fitness = 0.01f;

	sprintf(data_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 0; index < population->history.index; index++) {

		switch(history_type) {
		case GPR_HISTORY_FITNESS: {
			value = population->history.log[index];
			if (value<0) value=0;
			break;
		}
		case GPR_HISTORY_AVERAGE: {
			value = population->history.average[index];
			break;
		}
		case GPR_HISTORY_DIVERSITY: {
			value = population->history.diversity[index];
			break;
		}
		}

		fprintf(fp,"%d    %.10f\n",
				index*population->history.interval,
				value);
		/* record the maximum fitnes value */
		if (value > max_fitness) {
			max_fitness = value;
		}
		if ((index==0) ||
			(value < min_fitness)) {
			min_fitness = value;
		}
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [0:%d]\n",
			population->history.index*
			population->history.interval);
	fprintf(fp,"set yrange [%f:%f]\n",
			min_fitness,max_fitness*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Generation\"\n");
	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set ylabel \"Fitness\"\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set ylabel \"Average Fitness\"\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set ylabel \"Population Diversity\"\n");
		break;
	}
	}
	fprintf(fp,"%s","set grid\n");

	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set key right top\n");
		break;
	}
	}

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"plot \"%s\" using 1:2 notitle with lines\n",
			data_filename);
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* uses gnuplot to plot the fitness histogram
   for the given population */
int gpr_plot_fitness(gpr_population * population,
					 char * filename, char * title,
					 int image_width, int image_height)
{
	int index,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float min_fitness = 0;
	float max_fitness = 0.01f;
	int histogram_max=1;
	int histogram[GPR_HISTOGRAM_LEVELS];

	/* create the histogram */
	gpr_fitness_histogram(population, (int*)histogram,
						  GPR_HISTOGRAM_LEVELS,
						  &min_fitness, &max_fitness);

	if (max_fitness <= min_fitness) return 0;

	sprintf(data_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	for (index = 1; index < GPR_HISTOGRAM_LEVELS; index++) {
		if (index==1) {
			min_fitness = population->fitness[index];
			max_fitness = population->fitness[index]+0.001f;
		}
		else {
			if (population->fitness[index] > max_fitness) {
				max_fitness = population->fitness[index];
			}
			if (population->fitness[index] < min_fitness) {
				min_fitness = population->fitness[index];
			}
		}
	}

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 1; index < GPR_HISTOGRAM_LEVELS; index++) {
		fprintf(fp,"%f    %d\n",
				min_fitness + (index*(max_fitness-min_fitness)/
							   GPR_HISTOGRAM_LEVELS), histogram[index]);
		if (histogram[index]>histogram_max) {
			histogram_max = histogram[index];
		}
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [%f:%f]\n",min_fitness,max_fitness);
	fprintf(fp,"set yrange [0:%d]\n",histogram_max*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Fitness\"\n");
	fprintf(fp,"%s","set ylabel \"Instances\"\n");
	fprintf(fp,"%s","set grid\n");
	fprintf(fp,"%s","set key right top\n");

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"plot \"%s\" using 1:2 notitle with lines\n",
			data_filename);
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* uses gnuplot to plot the fitness histogram for the given system */
int gpr_plot_fitness_system(gpr_system * sys,
							char * filename, char * title,
							int image_width, int image_height)
{
	int index,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float min_fitness = 0;
	float max_fitness = 0.01f;
	int histogram_max=1;
	int histogram[GPR_HISTOGRAM_LEVELS];

	/* create the histogram */
	gpr_fitness_histogram_system(sys, (int*)histogram,
								 GPR_HISTOGRAM_LEVELS,
								 &min_fitness, &max_fitness);

	if (max_fitness <= min_fitness) return 0;

	sprintf(data_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 0; index < GPR_HISTOGRAM_LEVELS; index++) {
		fprintf(fp,"%f    %d\n",
				min_fitness + (index*(max_fitness-min_fitness)/
							   GPR_HISTOGRAM_LEVELS),
				histogram[index]);
		if (histogram[index]>histogram_max) {
			histogram_max = histogram[index];
		}
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [%f:%f]\n",min_fitness,max_fitness);
	fprintf(fp,"set yrange [0:%d]\n",histogram_max*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Fitness\"\n");
	fprintf(fp,"%s","set ylabel \"Instances\"\n");
	fprintf(fp,"%s","set grid\n");
	fprintf(fp,"%s","set key right top\n");

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"plot \"%s\" using 1:2 notitle with lines\n",
			data_filename);
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* uses gnuplot to plot the fitness history for the given system */
int gpr_plot_history_system(gpr_system * sys,
							int history_type,
							char * filename, char * title,
							int image_width, int image_height)
{
	int index,i,retval;
	FILE * fp;
	char data_filename[256];
	char plot_filename[256];
	char command_str[256];
	float min_value = 0;
	float max_value = 0.0001f;
	float value;

	sprintf(data_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.dat");
	sprintf(plot_filename,"%s%s",GPR_TEMP_DIRECTORY,"libgpr_data.plot");

	/* save the data */
	fp = fopen(data_filename,"w");
	if (!fp) return -1;
	for (index = 0; index < sys->island[0].history.index; index++) {
		fprintf(fp,"%d", index*sys->island[0].history.interval);
		for (i = 0; i < sys->size; i++) {

			switch(history_type) {
			case GPR_HISTORY_FITNESS: {
				value = sys->island[i].history.log[index];
				if (value<0) value = 0;
				break;
			}
			case GPR_HISTORY_AVERAGE: {
				value = sys->island[i].history.average[index];
				break;
			}
			case GPR_HISTORY_DIVERSITY: {
				value = sys->island[i].history.diversity[index];
				break;
			}
			}

			fprintf(fp,"    %.10f",value);
			/* record the maximum value */
			if (value > max_value) {
				max_value = value;
			}
			/* record the minimum value */
			if (((index==0) && (i==0)) ||
				(value < min_value)) {
				min_value = value;
			}
		}
		fprintf(fp,"%s","\n");
	}
	fclose(fp);

	/* create a plot file */
	fp = fopen(plot_filename,"w");
	if (!fp) return -1;
	fprintf(fp,"%s","reset\n");
	fprintf(fp,"set title \"%s\"\n",title);
	fprintf(fp,"set xrange [0:%d]\n",
			sys->island[0].history.index*
			sys->island[0].history.interval);
	fprintf(fp,"set yrange [%f:%f]\n",min_value,max_value*102/100);
	fprintf(fp,"%s","set lmargin 9\n");
	fprintf(fp,"%s","set rmargin 2\n");
	fprintf(fp,"%s","set xlabel \"Generation\"\n");
	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set ylabel \"Fitness\"\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set ylabel \"Average Fitness\"\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set ylabel \"Population Diversity\"\n");
		break;
	}
	}
	fprintf(fp,"%s","set grid\n");

	switch(history_type) {
	case GPR_HISTORY_FITNESS: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_AVERAGE: {
		fprintf(fp,"%s","set key right bottom\n");
		break;
	}
	case GPR_HISTORY_DIVERSITY: {
		fprintf(fp,"%s","set key right top\n");
		break;
	}
	}

	fprintf(fp,"set terminal png size %d,%d\n",
			image_width, image_height);
	fprintf(fp,"set output \"%s\"\n", filename);
	fprintf(fp,"%s","plot");

	for (i = 0; i < sys->size; i++) {
		fprintf(fp," \"%s\" using 1:%d title \"Island %d\" with lines",
				data_filename, i+2, i+1);
		if (i < sys->size-1) {
			fprintf(fp,"%s",",");
		}
	}
	fclose(fp);

	/* run gnuplot using the created files */
	sprintf(command_str,"gnuplot %s", plot_filename);
	retval = system(command_str); /* I assume this is synchronous */

	/* remove temporary files */
	sprintf(command_str,"rm %s %s", data_filename,plot_filename);
	retval = system(command_str);
	return retval;
}

/* returns a non zero value if the given function type is in the set */
int gpr_function_in_set(int function_type,
						int * instruction_set, int no_of_instructions)
{
	int i;

	for (i = 0; i < no_of_instructions; i++) {
		if (instruction_set[i]==function_type) {
			return 1;
		}
	}
	return 0;
}

/* Writes the given image buffer to a file.
   The image is expected to have three bytes per pixel */
int write_png_file(char * filename,
				   int width, int height,
				   unsigned char * buffer)
{
    png_t png;
    FILE * fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fprintf(stderr,
				"Could not open file %s for writing\n",
				filename);
        return 1;
    }
    fclose(fp);

    png_init(0,0);
    png_open_file_write(&png, filename);
    png_set_data(&png, width, height, 8, PNG_TRUECOLOR, buffer);
    png_close_file(&png);

    return 0;
}

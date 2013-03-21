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

#ifndef GPR_GLOBALS_H
#define GPR_GLOBALS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#define GPR_VERSION           "1.03"

/* define this to enable asserts */
#undef DEBUG

#define GPR_WEB               "https://github.com/fuzzgun/libgpr"

#define GPR_MAX_ARGUMENTS     4
#define GPR_DEFAULT_ARGUMENTS 2
#define GPR_TERMINALS         2

/* used to indicate a missing value */
#define GPR_MISSING_VALUE     -9999

/* number of fitness histogram levels */
#define GPR_HISTOGRAM_LEVELS  100

/* the maximum call depth for ADFs */
#define GPR_MAX_CALL_DEPTH    4

/* maximum age of a program in generations */
#define GPR_MAX_AGE 1000

/* the default branching probability for random trees */
#define GPR_DEFAULT_BRANCHING_PROB 0.5

/* values returned by conditional nodes (and, or, not) */
#define GPR_TRUE  1
#define GPR_FALSE 0

/* types of validation on the program tree */
#define GPR_VALIDATE_OK                                   0
#define GPR_VALIDATE_TREE_NOT_TERMINATED                 -1
#define GPR_VALIDATE_TREE_TOO_DEEP                       -2
#define GPR_VALIDATE_TERMINATOR_AT_ROOT                  -3
#define GPR_VALIDATE_FUNCTION_LESS_THAN_MIN              -4
#define GPR_VALIDATE_FUNCTION_TYPE_TOO_LARGE             -5
#define GPR_VALIDATE_FUNCTION_TYPE_NOT_IN_SET            -6
#define GPR_VALIDATE_CONNECTION_LESS_THAN_ZERO           -7
#define GPR_VALIDATE_CONNECTION_TOO_LARGE                -8
#define GPR_VALIDATE_ACTUATOR_CONNECTION_LESS_THAN_ZERO  -9
#define GPR_VALIDATE_ACTUATOR_CONNECTION_OUT_OF_RANGE    -10
#define GPR_VALIDATE_STATE_VALUE_OUT_OF_RANGE            -11
#define GPR_VALIDATE_ADF_PATTERN                         -12
#define GPR_VALIDATE_ADF_CONTAINS_VALUE                  -13
#define GPR_VALIDATE_ADF_CONTAINS_ARG                    -14
#define GPR_VALIDATE_ADF_ARG_OUT_OF_RANGE                -15
#define GPR_VALIDATE_ADF_WITHIN_ADF                      -16
#define GPR_VALIDATE_ADF_CONNECTIONS_OUT_OF_RANGE        -17
#define GPR_VALIDATE_ADF_NO_OF_ARGS                      -18

/* the maximum percentage by which terminal values can be mutated */
#define GPR_MUTATE_VALUE_PERCENT 2

/* the maximum/minimum terminal value */
#define GPR_MAX_CONSTANT         4096

/* maximum nodes when loading or saving a function */
#define GPR_MAX_NODES 5096

/* return values when loading a function */
#define GPR_LOAD_OK                        0
#define GPR_LOAD_MAX_NODES_REACHED        -1
#define GPR_LOAD_ARG_NOT_FOUND            -2
#define GPR_LOAD_CHILD_NOT_FOUND          -3
#define GPR_LOAD_PARENT_NOT_FOUND         -4
#define GPR_LOAD_VALUE_NOT_FOUND          -5
#define GPR_LOAD_FUNCTION_TYPE_NOT_FOUND  -6
#define GPR_LOAD_LINK_NOT_FOUND           -7
#define GPR_LOAD_NODE_NOT_FOUND           -8
#define GPR_LOAD_ARGC_NOT_FOUND           -9
#define GPR_LOAD_POPULATION_SIZE          -10

/* maximum number of trials during unit tests */
#define GPR_MAX_TESTS 100

/* temporary directory */
#define GPR_TEMP_DIRECTORY "/tmp/"

/* maximum recorded history in time steps */
#define GPR_MAX_HISTORY 10000

/* The radius of block coppies in self-modifying
   Cartesian genetic programming */
#define GPR_BLOCK_WIDTH 2

/* probability of function type mutations being ADFs */
#define GPR_ADF_PROB 0.1

/* learning rate for hebbian connections */
#define GPR_HEBBIAN_LEARNING_RATE 0.001

/* types of history plot */
#define GPR_HISTORY_FITNESS   0
#define GPR_HISTORY_AVERAGE   1
#define GPR_HISTORY_DIVERSITY 2

/* maximum number of modules in CGP */
#define GPRC_MAX_ADF_MODULES        10

/* maximum number of sensors for a CGP module */
#define GPRC_MAX_ADF_MODULE_SENSORS 10

/* the maximum number of genes when an ADF is created */
#define GPRC_MAX_ADF_GENES        30

/* the minimum number of genes for an ADF */
#define GPRC_MIN_ADF_GENES        2

#endif

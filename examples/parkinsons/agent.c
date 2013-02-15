/* Cartesian Genetic Program
   Evolved using libgpr
   https://launchpad.net/libgpr

   To compile:
     gcc -Wall -std=c99 -pedantic -o agent agent.c -lm
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

const int sensors = 19;
const int actuators = 2;
const int rows = 9;
const int columns = 10;
const int connections_per_gene = 5;
float * genome[5];
float * state[5];
int tick = 0;

/* Genome for ADF_module 0 */
float gene0[] = {22,85.600,17,8,9,6,12,-1,-1,-1,-1,-1,-1,-1,11,-35.080,17,0,18,14,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,-7.015,9,6,8,15,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,28,-2.820,6,22,25,22,24,-1,-1,-1,-1,-1,-1,-1,9,46.100,0,27,5,19,27,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,44,0.000,5,1,19,36,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,6,-81.420,19,21,33,34,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,-5.860,75,56,11,68,39,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,38,83};

/* State array for ADF_module 0 */
float state0[111];

/* Genome for ADF_module 1 */
float gene1[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,24,-26.040,7,8,8,8,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,10,-5.880,0,5,5,7,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,23,2.380,66,48,39,18,17,12,76.800,12,34,13,18,50,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,29};

/* State array for ADF_module 1 */
float state1[101];

/* Genome for ADF_module 2 */
float gene2[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,12,-2.260,5,6,7,2,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,29,70.476,8,18,3,9,14,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,-5.860,66,2,2,59,30,3,72.671,1,0,8,11,47,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,74};

/* State array for ADF_module 2 */
float state2[101];

/* Genome for ADF_module 3 */
float gene3[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,28,-74.300,2,4,6,4,18,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,70.840,45,20,22,27,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,28,-82.751,3,48,7,55,1,27,-74.740,6,34,4,5,68,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,84};

/* State array for ADF_module 3 */
float state3[101];

/* Genome for ADF_module 4 */
float gene4[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,26,-32.654,19,12,3,5,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,13,-36.639,5,19,2,0,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,26,96.130,32,8,17,6,8,11,6.499,21,3,2,9,19,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,64};

/* State array for ADF_module 4 */
float state4[101];

void get_inputs(int argc, char* argv[])
{
  int i;

  for (i = 0; i < argc; i++) {
    state[0][i] = atof(argv[i]);
  }
}

void run(int ADF_module)
{
  int row,col,n=0,i=0,j,k,g,ctr,src,dest,block_from,block_to,no_of_args;
  int sens,act;
  int call_ADF_module,itt;
  float * gp;

  if (ADF_module == 0) {
    sens = sensors;
    act = actuators;
  }
  else {
    sens = 10;
    act = 1;
  }

  for (col = 0; col < columns; col++) {
    for (row = 0; row < rows; row++, i++,n+=7) {
      if (genome[ADF_module][n] < 0) continue;
      gp = &genome[ADF_module][n];
      switch((int)gp[0]) {
      case 44: {
        if (ADF_module == 0) {
          call_ADF_module = 1 + (abs((int)gp[1])%4);
          no_of_args = 1 + (((int)gp[2])%10);
          if (no_of_args >= 5) {
            no_of_args = 5-1;
          }
          for (j = 0; j < no_of_args; j++) {
            state[call_ADF_module][j] = state[ADF_module][(int)gp[3+j]];
          }
          for (itt = 0; itt < 2; itt++) {
            run(call_ADF_module);
          }
          state[ADF_module][sens+i] =
            state[call_ADF_module][10+(rows*columns)];
        }
        break;
      }
      case 1: {
        state[ADF_module][sens+i] = (int)gp[1];
        break;
      }
      case 30: {
        no_of_args = 1 + (abs((int)gp[1])%4);
        state[ADF_module][sens+i] = 0;
        for (j = 0; j < no_of_args; j++) {
          state[ADF_module][sens+i] += state[ADF_module][(int)gp[2+j]];
        }
        state[ADF_module][sens+i] = 1.0f / (1.0f + exp(state[ADF_module][sens+i]));
        break;
      }
      case 2: {
        no_of_args = 1 + (abs((int)gp[1])%4);
        state[ADF_module][sens+i] = 0;
        for (j = 0; j < no_of_args; j++) {
          state[ADF_module][sens+i] += state[ADF_module][(int)gp[2+j]];
        }
        break;
      }
      case 3: {
        no_of_args = 1 + (abs((int)gp[1])%4);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 0; j < no_of_args; j++) {
          state[ADF_module][sens+i] -= state[ADF_module][(int)gp[2+j]];
        }
        break;
      }
      case 4: {
        state[ADF_module][sens+i] = -(int)state[ADF_module][(int)gp[2]];
        break;
      }
      case 5: {
        no_of_args = 1 + (abs((int)gp[1])%4);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 0; j < no_of_args; j++) {
          state[ADF_module][sens+i] *= state[ADF_module][(int)gp[2+j]];
        }
        break;
      }
      case 6: {
        if (state[ADF_module][(int)gp[3]] == 0) {
          state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        }
        else {
          state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]] / state[ADF_module][(int)gp[3]];
        }
        break;
      }
      case 7: {
        if ((int)state[ADF_module][(int)gp[3]] == 0) {
          state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        }
        else {
          state[ADF_module][sens+i] = (int)state[ADF_module][(int)gp[2]] % (int)state[ADF_module][(int)gp[3]];
        }
        break;
      }
      case 8: {
        state[ADF_module][sens+i] = floor(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 9: {
        no_of_args = 1 + (abs((int)gp[1])%4);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 1; j < no_of_args; j++) {
          state[ADF_module][sens+i] += state[ADF_module][(int)gp[2+j]];
        }
        state[ADF_module][sens+i] /= no_of_args;
        break;
      }
      case 10: {
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        break;
      }
      case 11: {
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        break;
      }
      case 12: {
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        break;
      }
      case 13: {
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        break;
      }
      case 14: {
        if (state[ADF_module][(int)gp[2]] > state[ADF_module][(int)gp[3]]) {
          state[ADF_module][sens+i] = (int)gp[1];
        }        else {          state[ADF_module][sens+i] = 0;
        }        break;
      }
      case 15: {
        if (state[ADF_module][(int)gp[2]] < state[ADF_module][(int)gp[3]]) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 16: {
        if (state[ADF_module][(int)gp[2]] == state[ADF_module][(int)gp[3]]) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 17: {
        if ((state[ADF_module][(int)gp[2]]>0) &&
            (state[ADF_module][(int)gp[3]]>0)) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 18: {
        if ((state[ADF_module][(int)gp[2]]>0) ||
            (state[ADF_module][(int)gp[3]]>0)) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 19: {
        if ((state[ADF_module][(int)gp[2]]>0) !=
            (state[ADF_module][(int)gp[3]]>0)) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 41: {
          state[ADF_module][sens+i] = 0.10 * state[ADF_module][(int)gp[2]] * state[ADF_module][(int)gp[3]];
        break;
      }
      case 22: {
        state[ADF_module][sens+i] = (float)exp(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 23: {
        state[ADF_module][sens+i] = (float)sqrt(fabs(state[ADF_module][(int)gp[2]]));
        break;
      }
      case 24: {
        state[ADF_module][sens+i] = (float)fabs(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 25: {
        state[ADF_module][sens+i] = (float)(sin(state[ADF_module][(int)gp[2]])*256);
        break;
      }
      case 26: {
        state[ADF_module][sens+i] = (float)asin(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 27: {
        state[ADF_module][sens+i] = (float)(cos(state[ADF_module][(int)gp[2]])*256);
        break;
      }
      case 28: {
        state[ADF_module][sens+i] = (float)acos(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 29: {
        state[ADF_module][sens+i] = (float)pow(state[ADF_module][(int)gp[2]],state[ADF_module][(int)gp[3]]);
        break;
      }
      case 31: {
        no_of_args = 1 + (abs((int)gp[1])%4);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 1; j < no_of_args; j++) {
          if (state[ADF_module][(int)gp[2+j]] < state[ADF_module][sens+i]) {
            state[ADF_module][sens+i] = state[ADF_module][(int)gp[2+j]];
          }
        }
        break;
      }
      case 32: {
        no_of_args = 1 + (abs((int)gp[1])%4);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 1; j < no_of_args; j++) {
          if (state[ADF_module][(int)gp[2+j]] > state[ADF_module][sens+i]) {
            state[ADF_module][sens+i] = state[ADF_module][(int)gp[2+j]];
          }
        }
        break;
      }
      case 33: {
        if ((gp[2] > sens) && (gp[3] > sens)) {
          src = ((int)gp[2]-sens)*7;
          dest = ((int)gp[3]-sens)*7;
          genome[ADF_module][dest] = genome[ADF_module][src];
        }
        break;
      }
      case 34: {
        if ((gp[2] > sens) && (gp[3] > sens)) {
          src = ((int)gp[2]-sens)*7;
          dest = ((int)gp[3]-sens)*7;
          genome[ADF_module][dest+1] = genome[ADF_module][src+1];
        }
        break;
      }
      case 35: {
        state[ADF_module][(int)gp[3]] = state[ADF_module][(int)gp[2]];
        break;
      }
      case 36: {
        block_from = (int)gp[2];
        block_to = (int)gp[3];
        if (block_from<block_to) {
          block_from = (int)gp[3];
          block_to = (int)gp[2];
        }
        k = block_to - 2;
        for (j = block_from - 2;
             j <= block_from + 2; j++,k++) {
          if ((j>sens) &&
              (k>sens) &&
              (j<i) && (k<i)) {
            for (g = 0; g < 7; g++) {
              genome[ADF_module][(j-sens)*7 + g] =
                genome[ADF_module][(k-sens)*7 + g];
            }
          }
        }
        break;
      }
      case 37: {
        if (gp[2] > sens) {
          src = ((int)gp[2] - sens) * 7;
          gp[3] = genome[ADF_module][src+2];
        };
        break;
      }
      case 38: {
        if (gp[3] > sens) {
          src = ((int)gp[3] - sens) * 7;
          gp[2] = genome[ADF_module][src+2];
        };
        break;
      }
      case 39: {
        if (gp[3] > sens) {
          src = ((int)gp[3] - sens) * 7;
          gp[2] = genome[ADF_module][src+3];
        };
        break;
      }
      case 40: {
        if (gp[2] > sens) {
          src = ((int)gp[2] - sens) * 7;
          gp[3] = genome[ADF_module][src+3];
        };
        break;
      }
      }
      /* prevent values from going out of range */
      if ((isnan(state[ADF_module][sens+i])) || (isinf(state[ADF_module][sens+i]))) {
        state[ADF_module][sens+i] = 0;
      }
      if (state[ADF_module][sens+i] > 4096) {
        state[ADF_module][sens+i] = 4096;
      }
      if (state[ADF_module][sens+i] < -4096) {
        state[ADF_module][sens+i] = -4096;
      }
    }
  }

  /* set the actuator values */
  ctr = sens + i;
  for (i = 0; i < act; i++, ctr++, n++) {
    state[ADF_module][ctr] = state[ADF_module][(int)genome[ADF_module][n]];
  }
}

void set_outputs()
{
  int i;

  for (i = 0; i < actuators; i++) {
    if (i > 0) printf(" ");
    printf("%.3f",state[0][sensors+(rows*columns)+i]);
  }
  printf("\n");
}

static void setup()
{
  int i,j;

  genome[0] = gene0;
  state[0] = state0;
  genome[1] = gene1;
  state[1] = state1;
  genome[2] = gene2;
  state[2] = state2;
  genome[3] = gene3;
  state[3] = state3;

  /* Clear the state array */
  for (i = 0; i < sensors + (rows*columns) + actuators; i++) {
    for (j = 0; j < 4; j++) {
      state[j][i] = 0;
    }
  }

}

static int read_stdin_args(char *argv[])
{
  char argstr[1024], value[1024], *retval;
  int i,ctr=0,field_index=0;
  FILE * fp = fopen("stdin","r");

  if (!fp) return 0;
  while (!feof(fp)) {
    retval = fgets(argstr,1023,fp);
    if (retval) {
      for (i = 0; i < strlen(argstr); i++) {
        if ((argstr[i]==' ') ||
            (argstr[i]==',') ||
            (argstr[i]==';') ||
            (i == strlen(argstr)-1)) {
          if (i == strlen(argstr)-1) {
            value[ctr++] = argstr[i];
          }
          value[ctr] = 0;
          ctr = 0;
          argv[field_index] = (char*)malloc(strlen(value)+2);
          sprintf(argv[field_index++],"%s",value);
        }
        else {
          value[ctr++] = argstr[i];
        }
      }
    }
  }
  fclose(fp);
  return field_index;
}

int main(int argc, char* argv[])
{
  char * stdin_args[1024];
  int i,no_of_stdin_args = read_stdin_args(stdin_args);

  setup();

  if (no_of_stdin_args>0) {
    if (no_of_stdin_args != sensors) {
      printf("Invalid number of arguments %d/%d\n",no_of_stdin_args,sensors);
      for (i = 0; i < no_of_stdin_args; i++) {
        free(stdin_args[i]);
      }
      return -1;
    }
    get_inputs(no_of_stdin_args,stdin_args);
    for (i = 0; i < no_of_stdin_args; i++) {
      free(stdin_args[i]);
    }
  }
  else {
    if (argc-1 != sensors) {
      printf("Invalid number of arguments %d/%d\n",argc-1,sensors);
      return -1;
    }
    get_inputs(argc-1,&argv[1]);
  }

  for (i = 0; i < 2; i++) {
    run(0);
  }
  set_outputs();

  /* Increment the time step */
  tick++;
  if (tick>32000) tick=0;

  return 0;
}


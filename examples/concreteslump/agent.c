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

const int sensors = 7;
const int actuators = 3;
const int rows = 9;
const int columns = 16;
const int connections_per_gene = 11;
float * genome[1];
float * state[1];
int tick = 0;

/* Genome for ADF_module 0 */
float gene0[] = {11,-11.140,5,2,0,3,3,4,2,4,1,3,2,8,49.678,0,1,6,4,0,1,1,6,3,2,5,9,30.692,1,1,6,3,1,4,2,6,0,4,3,25,-48.141,2,3,4,3,1,0,2,1,6,1,4,29,3.500,0,1,1,0,5,4,5,0,2,0,0,1,-20.398,3,4,0,6,3,5,2,5,0,5,4,19,37.858,5,2,2,4,4,3,4,1,3,1,0,23,91.768,1,6,6,3,2,5,2,2,5,6,3,29,-36.235,3,5,1,4,0,5,4,0,0,1,5,7,89.680,7,6,9,0,3,11,5,0,2,0,5,8,8.280,3,1,11,13,4,7,5,10,9,8,8,7,63.051,9,14,6,3,14,10,13,10,4,7,12,25,-27.140,12,6,13,4,3,2,9,7,11,6,5,29,18.595,11,2,11,4,12,6,3,15,12,13,0,16,-89.976,12,15,0,15,14,11,11,5,14,8,0,8,-43.851,10,11,0,14,5,6,13,12,14,8,10,27,13.100,6,0,5,11,0,12,7,11,6,10,0,25,31.239,14,5,13,8,13,2,14,15,11,1,14,22,-46.002,20,7,16,3,19,21,12,15,23,12,23,7,29.465,21,14,15,5,1,1,15,0,12,12,19,29,81.964,15,1,23,5,15,2,0,19,0,13,7,28,-17.340,9,0,21,18,20,10,14,21,21,18,11,17,-66.080,2,9,20,8,16,17,13,0,9,9,18,10,-59.435,1,14,19,12,11,15,23,21,9,8,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,22,93.043,20,13,21,19,15,12,12,18,6,6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,23,7.620,22,18,26,12,4,31,22,0,7,8,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,9,32.971,3,7,15,1,14,1,0,5,7,19,16,22,-74.138,28,0,4,8,18,30,28,1,12,1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,18,-31.624,20,1,23,21,30,21,24,30,33,11,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,21,28.131,0,24,1,29,2,33,5,6,22,2,2,19,75.199,40,21,7,6,6,0,2,8,3,6,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,43,-58.580,29,24,36,20,17,12,42,25,39,27,33,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,17,35.257,49,28,24,24,5,8,13,26,4,32,30,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,8,-26.636,4,24,12,41,47,15,0,19,16,29,3,23,-4.920,12,10,14,43,7,2,40,9,5,19,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,16,37.960,60,1,11,35,13,45,29,9,60,4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,37,65.680,2,27,6,8,7,45,20,29,5,59,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,8,90.759,0,67,6,19,34,48,29,69,26,32,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,30,-36.862,26,62,56,43,20,3,62,60,58,26,47,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,-33.915,6,14,70,76,73,18,5,11,44,63,70,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,31,-85.305,16,86,37,4,29,63,28,31,15,67,72,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,20,-15.902,90,30,0,34,9,48,24,24,22,24,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,22,20.885,36,60,33,49,47,86,69,19,83,60,71,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-30.774,26,49,75,43,41,29,44,107,34,36,92,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,19,-19.891,105,60,4,50,128,31,113,87,54,56,53,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,15.533,56,139,46,78,100,59,95,23,59,66,21,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,139,143,123};

/* State array for ADF_module 0 */
float state0[154];

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
    for (row = 0; row < rows; row++, i++,n+=13) {
      if (genome[ADF_module][n] < 0) continue;
      gp = &genome[ADF_module][n];
      switch((int)gp[0]) {
      case 46: {
        break;
      }
      case 23: {
        j = abs((int)state[ADF_module][(int)gp[2]] + (int)state[ADF_module][(int)gp[3]])
            %(rows*columns);
        state[ADF_module][sens+i] = state[ADF_module][sens+j];
        break;
      }
      case 22: {
        j = abs((int)state[ADF_module][(int)gp[3]])
            %(rows*columns);
        state[ADF_module][sens+i] = gp[1]*state[ADF_module][(int)gp[2]];
        state[ADF_module][sens+j] = state[ADF_module][sens+i];

        if (state[ADF_module][sens+j] > 4096) {
          state[ADF_module][sens+j] = 4096;
        }
        if (state[ADF_module][sens+j] < -4096) {
          state[ADF_module][sens+j] = -4096;
        }
        break;
      }
      case 1: {
        state[ADF_module][sens+i] = (int)gp[1];
        break;
      }
      case 32: {
        no_of_args = 1 + (abs((int)gp[1])%10);
        state[ADF_module][sens+i] = 0;
        for (j = 0; j < no_of_args; j++) {
          state[ADF_module][sens+i] += state[ADF_module][(int)gp[2+j]];
        }
        state[ADF_module][sens+i] = 1.0f / (1.0f + exp(state[ADF_module][sens+i]));
        break;
      }
      case 2: {
        no_of_args = 1 + (abs((int)gp[1])%10);
        state[ADF_module][sens+i] = 0;
        for (j = 0; j < no_of_args; j++) {
          state[ADF_module][sens+i] += state[ADF_module][(int)gp[2+j]];
        }
        break;
      }
      case 3: {
        no_of_args = 1 + (abs((int)gp[1])%10);
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
        no_of_args = 1 + (abs((int)gp[1])%10);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 0; j < no_of_args; j++) {
          state[ADF_module][sens+i] *= state[ADF_module][(int)gp[2+j]];
        }
        break;
      }
      case 6: {
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]] * gp[1];
        break;
      }
      case 7: {
        if (state[ADF_module][(int)gp[3]] == 0) {
          state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        }
        else {
          state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]] / state[ADF_module][(int)gp[3]];
        }
        break;
      }
      case 8: {
        if ((int)state[ADF_module][(int)gp[3]] == 0) {
          state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        }
        else {
          state[ADF_module][sens+i] = (int)state[ADF_module][(int)gp[2]] % (int)state[ADF_module][(int)gp[3]];
        }
        break;
      }
      case 9: {
        state[ADF_module][sens+i] = floor(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 10: {
        no_of_args = 1 + (abs((int)gp[1])%10);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 1; j < no_of_args; j++) {
          state[ADF_module][sens+i] += state[ADF_module][(int)gp[2+j]];
        }
        state[ADF_module][sens+i] /= no_of_args;
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
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        break;
      }
      case 15: {
        if (state[ADF_module][(int)gp[2]] > state[ADF_module][(int)gp[3]]) {
          state[ADF_module][sens+i] = (int)gp[1];
        }        else {          state[ADF_module][sens+i] = 0;
        }        break;
      }
      case 16: {
        if (state[ADF_module][(int)gp[2]] < state[ADF_module][(int)gp[3]]) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 17: {
        if (state[ADF_module][(int)gp[2]] == state[ADF_module][(int)gp[3]]) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 18: {
        if ((state[ADF_module][(int)gp[2]]>0) &&
            (state[ADF_module][(int)gp[3]]>0)) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 19: {
        if ((state[ADF_module][(int)gp[2]]>0) ||
            (state[ADF_module][(int)gp[3]]>0)) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 20: {
        if ((state[ADF_module][(int)gp[2]]>0) !=
            (state[ADF_module][(int)gp[3]]>0)) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 21: {
        if (((int)state[ADF_module][(int)gp[2]]) !=
            ((int)state[ADF_module][(int)gp[3]])) {
          state[ADF_module][sens+i] = (int)gp[1];
        }
        else {
          state[ADF_module][sens+i] = 0;
        }
        break;
      }
      case 43: {
          state[ADF_module][sens+i] = 0.10 * state[ADF_module][(int)gp[2]] * state[ADF_module][(int)gp[3]];
        break;
      }
      case 24: {
        state[ADF_module][sens+i] = (float)exp(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 25: {
        state[ADF_module][sens+i] = (float)sqrt(fabs(state[ADF_module][(int)gp[2]]));
        break;
      }
      case 26: {
        state[ADF_module][sens+i] = (float)fabs(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 27: {
        state[ADF_module][sens+i] = (float)(sin(state[ADF_module][(int)gp[2]])*256);
        break;
      }
      case 28: {
        state[ADF_module][sens+i] = (float)asin(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 29: {
        state[ADF_module][sens+i] = (float)(cos(state[ADF_module][(int)gp[2]])*256);
        break;
      }
      case 30: {
        state[ADF_module][sens+i] = (float)acos(state[ADF_module][(int)gp[2]]);
        break;
      }
      case 31: {
        state[ADF_module][sens+i] = (float)pow(state[ADF_module][(int)gp[2]],state[ADF_module][(int)gp[3]]);
        break;
      }
      case 33: {
        no_of_args = 1 + (abs((int)gp[1])%10);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 1; j < no_of_args; j++) {
          if (state[ADF_module][(int)gp[2+j]] < state[ADF_module][sens+i]) {
            state[ADF_module][sens+i] = state[ADF_module][(int)gp[2+j]];
          }
        }
        break;
      }
      case 34: {
        no_of_args = 1 + (abs((int)gp[1])%10);
        state[ADF_module][sens+i] = state[ADF_module][(int)gp[2]];
        for (j = 1; j < no_of_args; j++) {
          if (state[ADF_module][(int)gp[2+j]] > state[ADF_module][sens+i]) {
            state[ADF_module][sens+i] = state[ADF_module][(int)gp[2+j]];
          }
        }
        break;
      }
      case 35: {
        if (((int)gp[2] > sens) && ((int)gp[3] > sens)) {
          src = ((int)gp[2]-sens)*13;
          dest = ((int)gp[3]-sens)*13;
          genome[ADF_module][dest] = genome[ADF_module][src];
        }
        break;
      }
      case 36: {
        if (((int)gp[2] > sens) && ((int)gp[3] > sens)) {
          src = ((int)gp[2]-sens)*13;
          dest = ((int)gp[3]-sens)*13;
          genome[ADF_module][dest+1] = genome[ADF_module][src+1];
        }
        break;
      }
      case 37: {
        state[ADF_module][(int)gp[3]] = state[ADF_module][(int)gp[2]];
        break;
      }
      case 38: {
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
            for (g = 0; g < 13; g++) {
              genome[ADF_module][(j-sens)*13 + g] =
                genome[ADF_module][(k-sens)*13 + g];
            }
          }
        }
        break;
      }
      case 39: {
        if (gp[2] > sens) {
          src = ((int)gp[2] - sens) * 13;
          gp[3] = genome[ADF_module][src+2];
        };
        break;
      }
      case 40: {
        if (gp[3] > sens) {
          src = ((int)gp[3] - sens) * 13;
          gp[2] = genome[ADF_module][src+2];
        };
        break;
      }
      case 41: {
        if (gp[3] > sens) {
          src = ((int)gp[3] - sens) * 13;
          gp[2] = genome[ADF_module][src+3];
        };
        break;
      }
      case 42: {
        if (gp[2] > sens) {
          src = ((int)gp[2] - sens) * 13;
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


  /* Clear the state array */
  for (i = 0; i < sensors + (rows*columns) + actuators; i++) {
    for (j = 0; j < 0; j++) {
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


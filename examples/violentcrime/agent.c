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

const int sensors = 2;
const int actuators = 10;
const int rows = 9;
const int columns = 10;
const int connections_per_gene = 2;
int tick = 0;

/* The evolved genome.  Unused functions and parameters have been set to -1 */
float genome[] = {17,41.140,0,1,10,28.360,0,0,-1,-1,-1,-1,1,63.760,0,1,12,39.471,1,0,2,-19.880,1,0,3,2.260,0,1,11,-57.440,0,1,14,-4.260,0,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,7,96.160,5,5,2,6.029,0,7,-1,-1,-1,-1,7,-94.086,7,0,-1,-1,-1,-1,9,-69.720,8,6,-1,-1,-1,-1,-1,-1,-1,-1,10,-66.060,14,2,12,-16.640,11,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,82.360,3,15,4,16.580,5,13,4,-89.540,22,2,17,1.540,2,8,-1,-1,-1,-1,1,-99.540,14,18,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,5,24.300,5,32,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,9,49.460,8,27,-1,-1,-1,-1,-1,-1,-1,-1,12,90.300,35,19,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,14,-61.796,28,30,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,14,89.764,10,58,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,61,17,9,58,86,43,82,23,29,22};

/* State array */
float state[102];

void get_inputs(int argc, char* argv[])
{
  int i;

  for (i = 0; i < argc; i++) {
    state[i] = atof(argv[i]);
  }
}

void run()
{
  int row,col,n=0,i=0,j,k,g,ctr,src,dest,block_from,block_to;
  float * gp;

  for (col = 0; col < columns; col++) {
    for (row = 0; row < rows; row++, i++,n+=4) {
      if (genome[n] < 0) continue;
      gp = &genome[n];
      switch((int)gp[0]) {
      case 1: {
        state[sensors+i] = (int)gp[1];
        break;
      }
      case 2: {
        state[sensors+i] = state[(int)gp[2]] + state[(int)gp[3]];
        break;
      }
      case 3: {
        state[sensors+i] = state[(int)gp[2]] - state[(int)gp[3]];
        break;
      }
      case 4: {
        state[sensors+i] = -(int)state[(int)gp[2]];
        break;
      }
      case 5: {
        state[sensors+i] = (int)state[(int)gp[2]] * (int)state[(int)gp[3]];
        break;
      }
      case 6: {
        if (state[(int)gp[3]] == 0) {
          state[sensors+i] = state[(int)gp[2]];
        }
        else {
          state[sensors+i] = state[(int)gp[2]] / state[(int)gp[3]];
        }
        break;
      }
      case 7: {
        if ((int)state[(int)gp[3]] == 0) {
          state[sensors+i] = state[(int)gp[2]];
        }
        else {
          state[sensors+i] = (int)state[(int)gp[2]] % (int)state[(int)gp[3]];
        }
        break;
      }
      case 8: {
        state[sensors+i] = floor(state[(int)gp[2]] + state[(int)gp[3]]);
        break;
      }
      case 9: {
        state[sensors+i] = (int)((state[(int)gp[2]] + state[(int)gp[3]])/2);
        break;
      }
      case 10: {
        state[sensors+i] = state[(int)gp[2]];
        break;
      }
      case 11: {
        state[sensors+i] = state[(int)gp[2]];
        break;
      }
      case 12: {
        state[sensors+i] = state[(int)gp[3]];
        break;
      }
      case 13: {
        state[sensors+i] = state[(int)gp[3]];
        break;
      }
      case 14: {
        if (state[(int)gp[2]] > state[(int)gp[3]]) {
          state[sensors+i] = (int)gp[1];
        }        else {          state[sensors+i] = 0;
        }        break;
      }
      case 15: {
        if (state[(int)gp[2]] < state[(int)gp[3]]) {
          state[sensors+i] = (int)gp[1];
        }
        else {
          state[sensors+i] = 0;
        }
        break;
      }
      case 16: {
        if (state[(int)gp[2]] == state[(int)gp[3]]) {
          state[sensors+i] = (int)gp[1];
        }
        else {
          state[sensors+i] = 0;
        }
        break;
      }
      case 17: {
        if ((state[(int)gp[2]]>0) &&
            (state[(int)gp[3]]>0)) {
          state[sensors+i] = (int)gp[1];
        }
        else {
          state[sensors+i] = 0;
        }
        break;
      }
      case 18: {
        if ((state[(int)gp[2]]>0) ||
            (state[(int)gp[3]]>0)) {
          state[sensors+i] = (int)gp[1];
        }
        else {
          state[sensors+i] = 0;
        }
        break;
      }
      case 19: {
        if ((state[(int)gp[2]]>0) !=
            (state[(int)gp[3]]>0)) {
          state[sensors+i] = (int)gp[1];
        }
        else {
          state[sensors+i] = 0;
        }
        break;
      }
      case 40: {
          state[sensors+i] = 0.10 * state[(int)gp[2]] * state[(int)gp[3]];
        break;
      }
      case 22: {
        state[sensors+i] = (float)exp(state[(int)gp[2]]);
        break;
      }
      case 23: {
        state[sensors+i] = (float)sqrt(fabs(state[(int)gp[2]]));
        break;
      }
      case 24: {
        state[sensors+i] = (float)fabs(state[(int)gp[2]]);
        break;
      }
      case 25: {
        state[sensors+i] = (float)(sin(state[(int)gp[2]])*256);
        break;
      }
      case 26: {
        state[sensors+i] = (float)asin(state[(int)gp[2]]);
        break;
      }
      case 27: {
        state[sensors+i] = (float)(cos(state[(int)gp[2]])*256);
        break;
      }
      case 28: {
        state[sensors+i] = (float)acos(state[(int)gp[2]]);
        break;
      }
      case 29: {
        state[sensors+i] = (float)pow(state[(int)gp[2]],state[(int)gp[3]]);
        break;
      }
      case 30: {
        if (state[(int)gp[2]] < state[(int)gp[3]]) {
           state[sensors+i] = state[(int)gp[2]];
        }
        else {
          state[sensors+i] = state[(int)gp[3]];
        }
        break;
      }
      case 31: {
        if (state[(int)gp[2]] > state[(int)gp[3]]) {
          state[sensors+i] = state[(int)gp[2]];
        }
        else {
          state[sensors+i] = state[(int)gp[3]];
        }
        break;
      }
      case 32: {
        if ((gp[2] > sensors) && (gp[3] > sensors)) {
          src = ((int)gp[2]-sensors)*4;
          dest = ((int)gp[3]-sensors)*4;
          genome[dest] = genome[src];
        }
        break;
      }
      case 33: {
        if ((gp[2] > sensors) && (gp[3] > sensors)) {
          src = ((int)gp[2]-sensors)*4;
          dest = ((int)gp[3]-sensors)*4;
          genome[dest+1] = genome[src+1];
        }
        break;
      }
      case 34: {
        state[(int)gp[3]] = state[(int)gp[2]];
        break;
      }
      case 35: {
        block_from = (int)gp[2];
        block_to = (int)gp[3];
        if (block_from<block_to) {
          block_from = (int)gp[3];
          block_to = (int)gp[2];
        }
        k = block_to - 2;
        for (j = block_from - 2;
             j <= block_from + 2; j++,k++) {
          if ((j>sensors) &&
              (k>sensors) &&
              (j<i) && (k<i)) {
            for (g = 0; g < 4; g++) {
              genome[(j-sensors)*4 + g] =
                genome[(k-sensors)*4 + g];
            }
          }
        }
        break;
      }
      case 36: {
        if (gp[2] > sensors) {
          src = ((int)gp[2] - sensors) * 4;
          gp[3] = genome[src+2];
        };
        break;
      }
      case 37: {
        if (gp[3] > sensors) {
          src = ((int)gp[3] - sensors) * 4;
          gp[2] = genome[src+2];
        };
        break;
      }
      case 38: {
        if (gp[3] > sensors) {
          src = ((int)gp[3] - sensors) * 4;
          gp[2] = genome[src+3];
        };
        break;
      }
      case 39: {
        if (gp[2] > sensors) {
          src = ((int)gp[2] - sensors) * 4;
          gp[3] = genome[src+3];
        };
        break;
      }
      }
      /* prevent values from going out of range */
      if ((isnan(state[sensors+i])) || (isinf(state[sensors+i]))) {
        state[sensors+i] = 0;
      }
      if (state[sensors+i] > 4096) {
        state[sensors+i] = 4096;
      }
      if (state[sensors+i] < -4096) {
        state[sensors+i] = -4096;
      }
    }
  }

  /* set the actuator values */
  ctr = sensors + i;
  for (i = 0; i < actuators; i++, ctr++, n++) {
    state[ctr] = state[(int)genome[n]];
  }
}

void set_outputs()
{
  int i;

  for (i = 0; i < actuators; i++) {
    if (i > 0) printf(" ");
    printf("%.3f",state[sensors+(rows*columns)+i]);
  }
  printf("\n");
}

static void setup()
{
  int i;

  /* Clear the state array */
  for (i = 0; i < sensors + (rows*columns) + actuators; i++) {
    state[i] = 0;
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
    run();
  }
  set_outputs();

  /* Increment the time step */
  tick++;
  if (tick>32000) tick=0;

  return 0;
}


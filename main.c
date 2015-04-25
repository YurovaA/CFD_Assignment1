#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){
  int xSize = 3;
  int ySize = 4;
  double **U = matrix(0, xSize, 0, ySize);
  double **V = matrix(0, xSize, 0, ySize);
  double **P = matrix(0, xSize, 0, ySize);

  double U_default = 1.0;
  double V_default = 2.0;
  double P_default = 3.0;
  printf("Running init_uvp with test values: %f, %f, %f\n", U_default, V_default, P_default);
  init_uvp(1, 2, 3, xSize, ySize, U, V, P);

  printf("\n  U Contains:\n---------------\n");
  for (int i = 0; i < ySize; i++)
  {
    printf("[");
    for (int j = 0; j < xSize; j++)
    {
        printf(" %f ", U[i][j]);
    }
    printf("]\n");
  }
  printf("\n  V Contains:\n---------------\n");
  for (int i = 0; i < ySize; i++)
  {
    printf("[");
    for (int j = 0; j < xSize; j++)
    {
        printf(" %f ", V[i][j]);
    }
    printf("]\n");
  }
  printf("\n  P Contains:\n---------------\n");
  for (int i = 0; i < ySize; i++)
  {
    printf("[");
    for (int j = 0; j < xSize; j++)
    {
        printf(" %f ", P[i][j]);
    }
    printf("]\n");
  }

  printf("\nDone!\n");
  return -1;
}

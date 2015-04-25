#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "sor.h"
#include "boundary_val.h"
#include <stdio.h>

int parseCommandLine(int argn, char** args, char** inputFile)
{
    if ( argn < 2 ) // no arguments given, print usage statement
    {
        printf( "usage: %s filename", args[0] );
        return 0;
    }
    else
    {
        *inputFile = args[1]; // TODO: don't assume argv[1] always contains a filename since we may have other options in the future...
    }
    
    return 1;
}

void output_uvp(double **U, double **V, double **P, int imax, int jmax)
{
  printf("\n  U Contains:\n---------------\n");
  for (int i = 0; i < imax; i++)
  {
    printf("[");
    for (int j = 0; j < jmax; j++)
    {
        printf(" %f ", U[i][j]);
    }
    printf("]\n");
  }
  printf("\n  V Contains:\n---------------\n");
  for (int i = 0; i < imax; i++)
  {
    printf("[");
    for (int j = 0; j < jmax; j++)
    {
        printf(" %f ", V[i][j]);
    }
    printf("]\n");
  }
  printf("\n  P Contains:\n---------------\n");
  for (int i = 0; i < imax; i++)
  {
    printf("[");
    for (int j = 0; j < jmax; j++)
    {
        printf(" %f ", P[i][j]);
    }
    printf("]\n");
  }
}

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
  char *inputFile; // file name to read configuration parameters from

  // initial parameters given by inputFile
  double Re;
  double UI;
  double VI;
  double PI;
  double GX;
  double GY;
  double t_end;
  double xlength;
  double ylength;
  double dt;
  double dx;
  double dy;
  int imax;
  int jmax;
  double alpha;
  double omg;
  double tau;
  int itermax;
  double eps;
  double dt_value;

  // Data about the flow
  double **U; // speed in x-dir
  double **V; // speed in y-dir
  double **P; // pressure
  double **F; // Fn
  double **G; // Gn
  double ** RS; // right-hand side of pressure equation

  // current time
  double t = 0.0;

  if (!parseCommandLine(argn, args, &inputFile))
  {
    return -1; // exit if no arguments were given
  }
  else
  {
    read_parameters(inputFile, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength,
                    &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
                    &itermax, &eps, &dt_value);
  }

  // allocate space for U,V,P, F,G, RS
  U  = matrix(0, imax+1, 0, jmax+1);
  V  = matrix(0, imax+1, 0, jmax+1);
  P  = matrix(0, imax+1, 0, jmax+1);
  F  = matrix(0, imax+1, 0, jmax+1);
  G  = matrix(0, imax+1, 0, jmax+1);
  RS = matrix(0, imax+1, 0, jmax+1);

  // inititialize U,V,P
  init_uvp(UI, VI, PI, imax, jmax, U, V, P);
  

  // -------------------- //
  // main simulation loop //
  // -------------------- //

  while (t < t_end)
  {
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V); // get new timestep, dt
    boundaryvalues(imax,jmax,U,V);                // set boundary values in U & V for this timestep
    calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);  // compute Fn & Gn
    calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);      // compute right-hand side of pressure equation

    int it = 0;
    double res = 1.0; /// should this initial value depend on eps?
    while (it < itermax && res > eps)
    {
      sor(omg, dx, dy, imax, jmax, P, RS, &res);  // perform SOR
      it += 1;
    }

    calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);   // compute new velocities

    // ------------------------------------------------- //
    /* output U,V,P for current timestep here, if needed */
    // ------------------------------------------------- //

    t += dt;
  }

  
  // ----------------------- //
  /* output final U,V,P here */
  // ----------------------- //

  // free memory
  free_matrix(U, 0, imax+1, 0, jmax+1);
  free_matrix(V, 0, imax+1, 0, jmax+1);
  free_matrix(P, 0, imax+1, 0, jmax+1);
  free_matrix(F, 0, imax+1, 0, jmax+1);
  free_matrix(G, 0, imax+1, 0, jmax+1);
  free_matrix(RS, 0, imax+1, 0, jmax+1);


  printf("\nDone!\n");
  return -1;
}

#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "sor.h"
#include "boundary_val.h"
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>

int parseCommandLine(int argn, char** args, char** inputFile, char** outputLocation, int* outputAllFrames, int* outputBoundary)
{
    char* usage = "-f inputFile [-o outputDir] [-a]\n     -a will cause a .vtk file to be produced for every timestep\n     If no output directory is specified, the current directory will be used by default";
    char* outputDir;
    int opt;
    int setOutput = 0;
    int setInput = 0;

    if (argn > 2)
    {
        while ((opt = getopt(argn, args, "f:o:ab")) != -1)
        {
            switch (opt)
            {
                case 'f': // input filename
                    *inputFile = optarg;
                    setInput = 1;
                    break;
                case 'o': // output filename
                    outputDir = optarg;
                    setOutput = 1;
                    break;
                case 'a': // output VTK file for every timestep
                    *outputAllFrames = 1;
                    break;
                case 'b': // include boundary layer in output
                    *outputBoundary = 1;
                    break;
                case '?':
                    switch (optopt)
                    {
                        case 'f':
                        case 'o':
                            fprintf(stderr, "Option -%c requires a filename.\n", optopt);
                            break;
                        default:
                            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                            break;
                    }
                    return 0;
                default:  // unexpected option
                    printf("ERROR: Unexpected option \'%x\'\n", opt);
                    return 0;
                    break;
            }
        }
    }
    else if (argn == 2 && args[1][0] != '-') // assume only input file was given
    {
        *inputFile = args[1];
        setInput = 1;
    }
    else    // insufficient options given
    {
        printf( "usage: %s %s", args[0], usage );
        return 0;
    }

    if (!setInput)
    {
        fprintf(stderr, "You must provide an input file!\n");
        return 0;
    }

    // configure output location
    if (!setOutput)
    {
        *outputLocation = *inputFile;
    }
    else
    {
        asprintf(outputLocation, "%s/%s", outputDir, *inputFile);
    }

    printf("Input File: %s\nOutput Results To: %s.*.vtk\nOutput Result Every Step: %s\nInclude Boundaries in output: %s\n\n",
                        *inputFile, *outputLocation, *outputAllFrames ? "YES" : "NO",
                        *outputBoundary ? "YES" : "NO");
    return 1;
}

void output_uvp(char *outputLocation, int n, double xlength, double ylength, int imax, int jmax, double dx, double dy, double **U, double **V, double **P, int includeBoundary)
{
    if (includeBoundary)
    {
        write_vtkFile(outputLocation, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
    }
    else
    {
        write_vtkFile(outputLocation, n, xlength, ylength, imax-1, jmax-1, dx, dy, U+1, V+1, P+1);
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
  char *outputLocation; // where to store the results
  int outputAllSteps  = 0; // if nonzero, we will generate VTKs for ALL timesteps, not jsut the last one
  int outputBoundary = 0; // if non, we will include the boundary layer in our output

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

  if (!parseCommandLine(argn, args, &inputFile, &outputLocation, &outputAllSteps, &outputBoundary))
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
  int n = 0;
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
    if (outputAllSteps)
    {
        output_uvp(outputLocation, n, xlength, ylength, imax, jmax, dx, dy, U, V, P, outputBoundary);
    }

    t += dt;
    n++;
  }

  
  // ----------------------- //
  /* output final U,V,P here */
  // ----------------------- //
  output_uvp(outputLocation, n, xlength, ylength, imax, jmax, dx, dy, U, V, P, outputBoundary);

  // free memory
  free_matrix(U, 0, imax+1, 0, jmax+1);
  free_matrix(V, 0, imax+1, 0, jmax+1);
  free_matrix(P, 0, imax+1, 0, jmax+1);
  free_matrix(F, 0, imax+1, 0, jmax+1);
  free_matrix(G, 0, imax+1, 0, jmax+1);
  free_matrix(RS, 0, imax+1, 0, jmax+1);


  printf("\nDone!\n  Total Iterations: %d\n", n);
  return -1;
}

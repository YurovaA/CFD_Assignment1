#include "boundary_val.h"

void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
){
    int i, j;
    int imax_plus_one = imax + 1;
    int jmax_plus_one = jmax + 1;

    for (j = 0; j <= jmax_plus_one; j++){
        /* Left Wall */
        U[0][j] = 0;        /// Possibly redundant - always remains zero. Can be removed?
        V[0][j] = -V[1][j];

        /* Right Wall */
        U[imax][j] = 0;
        V[imax_plus_one][j] = -V[imax][j];
    }

    for (i = 0; i <= imax + 1; i++){
        /* Bottom Wall */
        U[i][0] = -U[i+1][0];
        V[i][0] = 0;        /// Possibly redundant - always remains zero. Can be removed?

        /* Top (Moving) Wall */
        U[i][jmax_plus_one] = 2 - U[i][jmax];   /// Corresponding to U_top = 1
        V[i][jmax] = 0;
    }
}

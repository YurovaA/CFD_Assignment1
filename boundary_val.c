#include <boundary_val.h>

void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
){
    int i, j;
    int imax_plus_two = imax + 2;
    int imax_plus_one = imax + 1;
    int jmax_plus_two = jmax + 2;
    int jmax_plus_one = jmax + 1;

    /* Left Wall */
    for (j = 0; j <= jmax + 1; j++){
        *U[j*imax_plus_two] = 0;       /// Possibly redundant - always remains zero. Can be removed?
        *V[j*imax_plus_two] = -*V[j*imax_plus_two + 1];
    }

    /* Right Wall */
    for (j = 0; j <= jmax + 1; j++){
        *U[j*imax_plus_two + imax] = 0;  /// Possibly redundant - always remains zero. Can be removed?
        *V[j*imax_plus_two + imax_plus_one] = -*V[j*imax_plus_two + imax];
    }

    /* Bottom Wall */
    for (i = 0; i <= imax + 1; i++){
        U[i] = -U[imax_plus_one+i];
        V[i] = 0;                       /// Possibly redundant - always remains zero. Can be removed?
    }

    /* Top (Moving) Wall */
    for (i = 0; i <= imax + 1; i++){
        U[jmax_plus_one*imax_plus_two + i] = 2.0 - U[jmax*imax_plus_two + i];
        V[jmax*imax_plus_two + i] = 0;  /// Possibly redundant - always remains zero. Can be removed?
    }
}

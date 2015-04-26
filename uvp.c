/**
 * Determines the value of U and G according to the formula
 */
void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
)
{
    return;
}


/**
 * This operation computes the right hand side of the pressure poisson equation.
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
    // start loops at 1 since the 0th entries correspond to boundary values outside the domain
    // NOTE: We don't need to set the 0th entries of any row/column because SOR also skips them
    for (int i = 1; i < imax; i++)
    {
        for (int j = 1; j < jmax; j++)
        {
            RS[i][j] = (1/dt) * ( ((F[i][j] - F[i-1][j]) / dx) 
                                + ((G[i][j] - G[i][j-1]) / dy) );
        }
    }
}


/**
 * Determines the maximal time step size.
 */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
    return;
}


/**
 * Calculates the new velocity values
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
    return;
}

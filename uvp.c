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
    return;
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

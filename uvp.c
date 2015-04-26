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
	double** d2udx2;
	double** d2udy2;
	double** du2dx;
	double** duvdy;
	d2udx2 = matrix(0, jmax+1, 0, imax+1);
	d2udy2 = matrix(0, jmax+1, 0, imax+1);
	du2dx = matrix(0, jmax+1, 0, imax+1);
	duvdy = matrix(0, jmax+1, 0, imax+1);
	int i;
	int j;
	// calculating d^2u/dx^2 in all points
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			d2udx2[i][j] = (U[i+1][j] - 2*U[i][j] + U[i-1][j])/(dx*dx);
		};
	};
	// calculating d^2u/dy^2 in all points
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			d2udy2[i][j] = (U[i][j+1] - 2*U[i][j] + U[i][j-1])/(dy*dy);
		};
	};
	// calculating du^2/dx in all points
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			du2dx[i][j] = (pow(((U[i][j] + U[i+1][j])/2),2)-pow(((U[i-1][j] + U[i][j])/2),2))/dx;
			du2dx[i][j] += alpha/dx*((U[i][j]+U[i+1][j])*(U[i][j] - U[i+1][j]) - (U[i-1][j]+U[i][j])*(U[i-1][j] - U[i][j]))/4;
		};
	};
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			duvdy[i][j] = 0.25*((V[i][j] + V[i+1][j])*(U[i][j] + V[i][j+1])-(V[i][j-1] + V[i+1][j-1])*(U[i][j-1] + U[i][j]))/dy;
			duvdy[i][j] += alpha/dy*((V[i][j]+V[i+1][j])*(U[i][j] - U[i][j+1]) - (V[i][j-1]+V[i+1][j-1])*(U[i][j-1] - U[i][j]))/4;
		};
	};
	F[i][j] = U[i][j] + dt*(1/Re*(d2udx2[i][j] + d2udy2[i][j])-du2dx[i][j] - duvdy[i][j] + GX);

	double** d2vdx2;
	double** d2vdy2;
	double** dv2dy;
	double** duvdx;
	d2vdx2 = matrix(0, jmax+1, 0, imax+1);
	d2vdy2 = matrix(0, jmax+1, 0, imax+1);
	dv2dy = matrix(0, jmax+1, 0, imax+1);
	duvdx = matrix(0, jmax+1, 0, imax+1);
	// calculating d^2v/dx^2 in all points
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			d2vdx2[i][j] = (V[i+1][j] - 2*V[i][j] + V[i-1][j])/(dx*dx);
		};
	};
	// calculating d^2v/dy^2 in all points
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			d2vdy2[i][j] = (V[i][j+1] - 2*V[i][j] + V[i][j-1])/(dy*dy);
		};
	};
	// calculating dv^2/dy in all points
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			dv2dy[i][j] = (pow(((V[i][j] + V[i][j+1])/2),2)-pow(((V[i][j-1] + V[i][j])/2),2))/dy;
			dv2dy[i][j] += alpha/dy*((V[i][j]+V[i][j+1])*(V[i][j] - V[i][j+1]) - (V[i][j-1]+V[i][j])*(V[i][j-1] - V[i][j]))/4;
		};
	};
	// calculating duv/dx in all points
	for (i = 1; i<imax-1; i++){
		for (j = 1; j<jmax; j++){
			duvdx[i][j] = 0.25*((U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j])-(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j]))/dx;
			duvdx[i][j] += alpha/dy*((U[i][j]+U[i][j+1])*(V[i][j] - V[i+1][j]) - (U[i-1][j]+U[i-1][j+1])*(V[i-1][j] - V[i][j]))/4;
		};
	};
	G[i][j] = V[i][j] + dt*(1/Re*(d2vdx2[i][j] + d2vdy2[i][j])-dv2dy[i][j] - duvdx[i][j] + GY);
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

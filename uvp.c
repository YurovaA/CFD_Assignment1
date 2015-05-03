/**
 * Determines the value of U and G according to the formula
 */

#include "helper.h"
#include "math.h"

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

	double inv_dxdx = 1/(dx*dx);
	double inv_dydy = 1/(dy*dy);
	double inv_dx = 1/dx;
	double inv_dy = 1/dy;
	double alpha_by_dx = alpha/dx;
	double alpha_by_dy = alpha/dy;
	double inv_Re = 1/Re;

	// calculating d^2u/dx^2 in all points
	for (i = 1; i<=imax-1; i++){
		for (j = 1; j<=jmax; j++){
			d2udx2[i][j] = (U[i+1][j] - 2*U[i][j] + U[i-1][j])*inv_dxdx;
		};
	};
	// calculating d^2u/dy^2 in all points
	for (i = 1; i<=imax-1; i++){
		for (j = 1; j<=jmax; j++){
			d2udy2[i][j] = (U[i][j+1] - 2*U[i][j] + U[i][j-1])*inv_dydy;
		};
	};
	// calculating du^2/dx in all points
	for (i = 1; i<=imax-1; i++){
		for (j = 1; j<=jmax; j++){
			du2dx[i][j] = (pow(((U[i][j] + U[i+1][j])/2),2)-pow(((U[i-1][j] + U[i][j])/2),2))*inv_dx;
			du2dx[i][j] += (0.25 * alpha_by_dx)*(fabs(U[i][j]+U[i+1][j])*(U[i][j] - U[i+1][j]) - fabs(U[i-1][j]+U[i][j])*(U[i-1][j] - U[i][j]));
		};
	};
	for (i = 1; i<=imax-1; i++){
		for (j = 1; j<=jmax; j++){
			duvdy[i][j] = 0.25*((V[i][j] + V[i+1][j])*(U[i][j] + U[i][j+1])-(V[i][j-1] + V[i+1][j-1])*(U[i][j-1] + U[i][j]))*inv_dy;
			duvdy[i][j] += (0.25*alpha_by_dy)*(fabs(V[i][j]+V[i+1][j])*(U[i][j] - U[i][j+1]) - fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1] - U[i][j]));
		};
	};
	for (i = 1; i<=imax-1; i++){
		for (j = 1; j<=jmax; j++){
			F[i][j] = U[i][j] + dt*(inv_Re*(d2udx2[i][j] + d2udy2[i][j])-du2dx[i][j] - duvdy[i][j] + GX);
		};
	};

	// update boundaries
	for (j = 1; j <= jmax; j++)
	{
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}


	double** d2vdx2;
	double** d2vdy2;
	double** dv2dy;
	double** duvdx;
	d2vdx2 = matrix(0, jmax+1, 0, imax+1);
	d2vdy2 = matrix(0, jmax+1, 0, imax+1);
	dv2dy = matrix(0, jmax+1, 0, imax+1);
	duvdx = matrix(0, jmax+1, 0, imax+1);
	// calculating d^2v/dx^2 in all points
	for (i = 1; i<=imax; i++){
		for (j = 1; j<=jmax-1; j++){
			d2vdx2[i][j] = (V[i+1][j] - 2*V[i][j] + V[i-1][j])*inv_dxdx;
		};
	};
	// calculating d^2v/dy^2 in all points
	for (i = 1; i<=imax; i++){
		for (j = 1; j<=jmax-1; j++){
			d2vdy2[i][j] = (V[i][j+1] - 2*V[i][j] + V[i][j-1])*inv_dydy;
		};
	};
	// calculating dv^2/dy in all points
	for (i = 1; i<=imax; i++){
		for (j = 1; j<=jmax-1; j++){
			dv2dy[i][j] = (pow(((V[i][j] + V[i][j+1])/2),2)-pow(((V[i][j-1] + V[i][j])/2),2))*inv_dy;
			dv2dy[i][j] += (0.25*alpha_by_dy)*(fabs(V[i][j]+V[i][j+1])*(V[i][j] - V[i][j+1]) - fabs(V[i][j-1]+V[i][j])*(V[i][j-1] - V[i][j]));
		};
	};
	// calculating duv/dx in all points
	for (i = 1; i<=imax; i++){
		for (j = 1; j<=jmax-1; j++){
			duvdx[i][j] = 0.25*((U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j])-(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j]))*inv_dx;
			duvdx[i][j] += (0.25*alpha_by_dy)*(fabs(U[i][j]+U[i][j+1])*(V[i][j] - V[i+1][j]) - fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j] - V[i][j]));
		};
	};
	for (i = 1; i<=imax; i++){
		for (j = 1; j<=jmax-1; j++){
			G[i][j] = V[i][j] + dt*(inv_Re*(d2vdx2[i][j] + d2vdy2[i][j])-dv2dy[i][j] - duvdx[i][j] + GY);
		};
	};

	// update boundaries
	for (i = 1; i <= imax; i++)
	{
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}

	// cleanup
	free_matrix(d2udx2, 0, jmax+1, 0, imax+1);
	free_matrix(d2udy2, 0, jmax+1, 0, imax+1);
	free_matrix(du2dx,  0, jmax+1, 0, imax+1);
	free_matrix(duvdy,  0, jmax+1, 0, imax+1);
	free_matrix(d2vdx2, 0, jmax+1, 0, imax+1);
	free_matrix(d2vdy2, 0, jmax+1, 0, imax+1);
	free_matrix(dv2dy,  0, jmax+1, 0, imax+1);
	free_matrix(duvdx,  0, jmax+1, 0, imax+1);
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
    double inv_dt = 1/dt;
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            RS[i][j] = inv_dt * ( ((F[i][j] - F[i-1][j]) / dx)
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
	double umax = 0;
	double vmax = 0;

	for (int i = 1; i <= imax; i++)
	{
		for (int j = 1; j <= jmax; j++)
		{
			if (fabs(U[i][j]) > umax)
				umax = fabs(U[i][j]);
			if (fabs(V[i][j]) > vmax)
				vmax = fabs(V[i][j]);
		}
	}

	*dt = tau * fmin( (Re / 2) * (1 / (((1 / (dx*dx)) + (1 / (dy*dy))))),
               fmin( ((double)dx/umax), ((double)dy/vmax) ));
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
	int i;
	int j;
	double dt_by_dx = dt/dx;
	double dt_by_dy = dt/dy;
	for (i = 1; i<=imax-1; i++){
		for (j = 1; j<=jmax; j++){
			U[i][j] = F[i][j] - dt_by_dx*(P[i+1][j] - P[i][j]);
	    }
	}
	for (i = 1; i<=imax; i++){
		for (j = 1; j<=jmax-1; j++){
			V[i][j] = G[i][j] - dt_by_dy*(P[i][j+1] - P[i][j]);
		}
	}
}

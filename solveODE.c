#include "bulkODE.h"
#include "solveODE.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#define ARGS t,y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],w
const int dim=8;

int func (double t, const double y[], double f[], void *params)
{
	double w=*(double *)params;
	f[0] = y[1];
	f[1] = f1(ARGS);
	f[2] = y[3];
	f[3] = f3(ARGS);
	f[4] = y[5];
	f[5] = f5(ARGS);
	f[6] = y[6];
	f[7] = f7(ARGS);

	return GSL_SUCCESS;
}/*
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	double w=*(double *)params;

	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, dim, dim);
	gsl_matrix * m = &dfdy_mat.matrix; 

	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 0, 2, 0.0);
	gsl_matrix_set (m, 0, 3, 0.0);
	gsl_matrix_set (m, 0, 4, 0.0);
	gsl_matrix_set (m, 0, 5, 0.0);

	gsl_matrix_set (m, 1, 0, j1_0(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 1, 1, j1_1(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 1, 2, j1_2(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 1, 3, j1_3(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 1, 4, j1_4(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 1, 5, j1_5(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );

	gsl_matrix_set (m, 2, 0, 0.0);
	gsl_matrix_set (m, 2, 1, 0.0);
	gsl_matrix_set (m, 2, 2, 0.0);
	gsl_matrix_set (m, 2, 3, 1.0);
	gsl_matrix_set (m, 2, 4, 0.0);
	gsl_matrix_set (m, 2, 5, 0.0);

	gsl_matrix_set (m, 3, 0, j3_0(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 3, 1, j3_1(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 3, 2, j3_2(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 3, 3, j3_3(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 3, 4, j3_4(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 3, 5, j3_5(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );

	gsl_matrix_set (m, 4, 0, 0.0);
	gsl_matrix_set (m, 4, 1, 0.0);
	gsl_matrix_set (m, 4, 2, 0.0);
	gsl_matrix_set (m, 4, 3, 0.0);
	gsl_matrix_set (m, 4, 4, 0.0);
	gsl_matrix_set (m, 4, 5, 1.0);

	gsl_matrix_set (m, 5, 0, j5_0(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 5, 1, j5_1(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 5, 2, j5_2(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 5, 3, j5_3(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 5, 4, j5_4(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );
	gsl_matrix_set (m, 5, 5, j5_5(t,y[0],y[1],y[2],y[3],y[4],y[5],w) );


	dfdt[0] =0.0; 
	dfdt[1] = dfdz1(t,y[0],y[1],y[2],y[3],y[4],y[5],w);
	dfdt[2] =0.0;
	dfdt[3] = dfdz3(t,y[0],y[1],y[2],y[3],y[4],y[5],w);
	dfdt[4] = 0.0;
	dfdt[5] = dfdz5(t,y[0],y[1],y[2],y[3],y[4],y[5],w);

	return GSL_SUCCESS;
}*/

int solveODE(double * start, double eps, int n, double *out,double *params)
{

	int i,j;
	//for(i=0; i<30; i++)
	//	printf ("skit: %f\n", out[i]);
	//                       yprim,jac ,dim, param
	gsl_odeiv2_system sys = {func, NULL, dim, params};
	//gsl_odeiv2_system sys = {func, jac, dim, params};
	//                             sys,driver type           ,initial h, abs err, rel err
	gsl_odeiv2_driver * d = 
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, -1e-8, 1e-8, 0.0);
	//gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, -1e-8, 1e-8, 0.0);
	//gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp, -1e-8, 1e-8, 0.0);
	double z0 = 1.0-eps, z1 = eps;
	double z=z0;
	double y[dim];
	for (i=0; i<dim; i++)
	{
		y[i]=start[i];
		out[i]=start[i];
	}	
	for (i = 1; i < n; i++)
	{
		double zi = z0+ i * (z1-z0) / ((double) (n-1));
		int status = gsl_odeiv2_driver_apply (d, &z, zi, y);

		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			  break;
		}
		for(j=0;j<dim;j++)
		{
			out[i*dim+j]=y[j];
		}
			//printf ("%.5e %.5e %.5e\n", z, y[0], y[1]);
	}

	gsl_odeiv2_driver_free (d);
	return 0;
}

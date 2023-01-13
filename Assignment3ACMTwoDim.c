#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define grid 129
#define gridheight 513

int main (void)
{
	double u[grid][gridheight+1], un[grid][gridheight+1], uc[grid][gridheight+1], un2[grid][gridheight+1];
	double v[grid+1][gridheight], vn[grid+1][gridheight], vc[grid+1][gridheight], vn2[grid+1][gridheight];
	double p[grid+1][gridheight+1], pn[grid+1][gridheight+1], pc[grid+1][gridheight+1], pn2[grid+1][gridheight+1];
	double m[grid+1][gridheight+1];
	int i, j, step;
	double dx, dy, dt, tau, delta, error, Re;
	int maxstep = 1000000;
	double residual[maxstep], umomresidual[maxstep], vmomresidual[maxstep];
	step =1;
	dx = 1.0/(grid-1);
	dy = dx;
	dt = 0.0005;
	delta = 1.0;
	error = 1.0;
	Re = 3200.0;
	clock_t start, end;
	double cpu_time_used;
	double CFL1, CFL2;
	CFL1 = dt/dx;
	CFL2 = (1./Re)*dt/(dx*dx);
	printf("CFL numbers are %lf and %lf\n", CFL1, CFL2);
	double WU;
	// Initializing u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(gridheight); j++)
			{
				u[i][j] = 0.0;
				u[i][gridheight] = 1.0;
				u[i][gridheight-1] = 1.0;
			}
		}
		
	// Initializing v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(gridheight-1); j++)
			{
				v[i][j] = 0.0;
			}
		}
		
	// Initializing p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(gridheight); j++)
			{
				p[i][j] = 1.0;
			}
		}
	
	start = clock();
	while (error > 0.000001)
	{
		// Solve u-momentum
		for (i=1; i<=(grid-2); i++)
		{
			for (j=1; j<=(gridheight-1); j++)
			{
				un[i][j] = u[i][j] - dt*(  (u[i+1][j]*u[i+1][j]-u[i-1][j]*u[i-1][j])/2.0/dx 
							+0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]) )/dy  )
								- dt/dx*(p[i+1][j]-p[i][j]) 
									+ dt*1.0/Re*( (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx/dx +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy/dy );
			}
		}
		
		for (j=1; j<=(gridheight-1); j++)
		{
			un[0][j] = 0.0;
			un[grid-1][j] = 0.0;
		}
		
		for (i=0; i<=(grid-1); i++)
		{
			un[i][0] = -un[i][1];
			un[i][gridheight] = 2.0 - un[i][gridheight-1];
		}
		
		
		// Solves v-momentum
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(gridheight-2); j++)
			{
				vn[i][j] = v[i][j] - dt* ( 0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j]) )/dx 
							+(v[i][j+1]*v[i][j+1]-v[i][j-1]*v[i][j-1])/2.0/dy ) 
								- dt/dy*(p[i][j+1]-p[i][j]) 
									+ dt*1.0/Re*( (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy/dy );
			}
		}
		
		for (j=1; j<=(gridheight-2); j++)
		{
			vn[0][j] = -vn[1][j];
			vn[grid][j] = -vn[grid-1][j];
		}		

		for (i=0; i<=(grid); i++)
		{
			vn[i][0] = 0.0;
			vn[i][gridheight-1] = 0.0;
		}		
	
		// Solves continuity equation
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(gridheight-1); j++)
			{
				pn[i][j] = p[i][j]-dt*delta*(  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] ) /dy  );
			}
		}
		
		for (i=1; i<=(grid-1); i++)
		{
			pn[i][0] = pn[i][1];
			pn[i][gridheight] = pn[i][gridheight-1];
		}
		
		for (j=0; j<=(gridheight); j++)
		{
			pn[0][j] = pn[1][j];
			pn[grid][j] = pn[grid-1][j];
		}		
		
		// Displaying error
		error = 0.0;
		
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(gridheight-1); j++)
			{
				m[i][j] = (  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] )/dy  );
				error = error + fabs(m[i][j]);
			}
		}
		// residual[step] = log10(error);
		
		if (step%10000 ==1)
		{
	    printf("Error is %5.10lf for the step %d\n", error, step);
		}
		
		
		// Iterating u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(gridheight); j++)
			{
				u[i][j] = un[i][j];
			}
		}
		
		// Iterating v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(gridheight-1); j++)
			{
				v[i][j] = vn[i][j];
			}
		}
		
		// Iterating p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(gridheight); j++)
			{
				p[i][j] = pn[i][j];
			}
		}

		step = step + 1;
		WU = WU + 1;
	}
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	
	for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(gridheight-1); j++)
			{	
			    uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
                vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
                pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
			}
		}
	
	
	printf("CPU Time = %lfs\n", cpu_time_used);
	printf("Workunits = %lf\n", WU);
	// OUTPUT DATA
	FILE *fout2, *fout3, *fout4, *fout5, *fout6, *fout7;
	fout2 = fopen("UVP.plt","w+t");
	fout3 = fopen("Central_U.plt","w+t");
	fout7 = fopen("Central_V.plt","w+t");
	fout4 = fopen("Residual.plt","w+t");
	fout5 = fopen("ResidualUMom.plt","w+t");
	fout6 = fopen("ResidualVMom.plt","w+t");
	
	if ( fout2 == NULL )
	{
    printf("\nERROR when opening file\n");
    fclose( fout2 );
	}

  else
	{
	fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
	fprintf( fout2, "ZONE  F=POINT\n");
	fprintf( fout2, "I=%d, J=%d\n", grid, gridheight );

	for ( j = 0 ; j < (gridheight) ; j++ )
	{
    for ( i = 0 ; i < (grid) ; i++ )
    {
		double xpos, ypos;
		xpos = i*dx;
		ypos = j*dy;

		fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, uc[i][j], vc[i][j], pc[i][j] );
    }
	}
	}

	fclose( fout2 );
	
	// CENTRAL --U
	  fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
	  fprintf(fout3, "ZONE F=POINT\n");
	  fprintf(fout3, "I=%d\n", gridheight );

	  for ( j = 0 ; j < gridheight ; j++ )
	  {
		double ypos;
		ypos = (double) j*dy;

		fprintf( fout3, "%5.8lf\t%5.8lf\n", (uc[(grid-1)/2][j]), ypos );
	  }
	 
	// CENTRAL --V
	  fprintf(fout7, "VARIABLES=\"V\",\"X\"\n");
	  fprintf(fout7, "ZONE F=POINT\n");
	  fprintf(fout7, "I=%d\n", grid );

	  for ( i = 0 ; i < grid ; i++ )
	  {
		double xpos;
		xpos = (double) i*dx;

		fprintf( fout7, "%5.8lf\t%5.8lf\n", (vc[i][(gridheight-1)/2]), xpos );
	  }	 
	// Residual
	  // fprintf(fout4, "VARIABLES=\"Workunits\",\"Error\"\n");
	  // fprintf(fout4, "ZONE F=POINT\n");
	  // fprintf(fout4, "I=%d\n", maxstep );

	  // for ( j = 0 ; j < maxstep ; j++ )
	  // {
		// fprintf( fout4, "%d\t%5.18lf\n", j, residual[j] );
	  // }
	  
	// // U Momentum Residual
	  // fprintf(fout5, "VARIABLES=\"Step\",\"U-Res\"\n");
	  // fprintf(fout5, "ZONE F=POINT\n");
	  // fprintf(fout5, "I=%d\n", maxstep );

	  // for ( j = 0 ; j < maxstep ; j++ )
	  // {
		// fprintf( fout5, "%d\t%5.18lf\n", j, umomresidual[j] );
	  // }
	  
	// // V Momentum Residual
	  // fprintf(fout6, "VARIABLES=\"Step\",\"V-Res\"\n");
	  // fprintf(fout6, "ZONE F=POINT\n");
	  // fprintf(fout6, "I=%d\n", maxstep );

	  // for ( j = 0 ; j < maxstep ; j++ )
	  // {
		// fprintf( fout6, "%d\t%5.18lf\n", j, vmomresidual[j] );
	  // }  

}
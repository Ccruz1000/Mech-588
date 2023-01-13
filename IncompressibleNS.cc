// Include Statements 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>



int main()
{
	// Define variables and initial conditions 
int gridsize = 10; // Uniform cartesian mesh 

double dx, dy, dt, tau, delta, error, Re;
int i, j, step;
step = 1; 
dx = 1.0 / (gridsize - 1);
dy = 1.0 / (gridsize - 1);
dt = 1e-3;
delta = 4.5;
error = 1.0;
Re = 100.0;
double errtol = 1e-6;
double u[gridsize][gridsize + 1], un[gridsize][gridsize + 1], uc[gridsize][gridsize]; // u is previous time step, unew is timestep being solved for, uc is u at the center point of the grid
double v[gridsize + 1][gridsize], vn[gridsize + 1][gridsize], vc[gridsize][gridsize];
double p[gridsize + 1][gridsize + 1], pn[gridsize + 1][gridsize + 1], pc[gridsize][gridsize]; 
double m[gridsize + 1][gridsize + 1];
// Initialize u
	for(i = 0; i < gridsize; i++)
	{
		for(j = 0; j <= gridsize; j++)
		{
			u[i][j] = 0.0;
			u[i][gridsize] = 1.0;
			u[i][gridsize - 1] = 1.0;
		}
	}

// Initialize v
	for(i = 0; i <= gridsize; i++)
	{
		for(j = 0; j < gridsize; j++)
		{
			v[i][j] = 0;
		}
	}

// Initialize p
	for(i = 0; i <= gridsize; i++)
	{
		for(j = 0; j <= gridsize; j++)
		{
			p[i][j] = 0;
		}
	}

// Time loop 
while(error > errtol)
{
	// Solve u-momentum
	// Internal points using central difference 
	for(i = 1; i < gridsize - 1; i++)
	{
		for(j = 1; j < gridsize; j++)
		{
			un[i][j] = u[i][j] - dt * ((((u[i + 1][j] * u[i + 1][j]) - (u[i - 1][j] * u[i - 1][j])) / (2.0 * dx)) +
					   0.25 * ((((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j])) - (u[i][j] + u[i][j - 1]) * (v[i + 1][j - 1] + v[i][j - 1])) / dy)) - 
					   (dt / dx) * (p[i + 1][j] - p[i][j]) +
					   dt * (1.0/Re) * ((u[i + 1][j] - 2 * u[i][j] + u[i - 1][j] / (dx * dx)) + (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (dy * dy)); 
		}
	}

	// Boundary conditions
	for(j = 1; j < gridsize; j ++)
	{
		un[0][j] = 0.0;
		un[gridsize - 1][j] = 0.0;
	}

	for(i = 0; i < gridsize; i ++)
	{
		un[i][0] = -un[i][1];
		un[i][gridsize] = 2.0 - un[i][gridsize - 1];
	}

	// Solve v-momentum
	// Internal points using central difference
	for(i = 1; i < gridsize; i++)
	{
		for(j = 1; j < gridsize - 1; j++)
		{
			vn[i][j] = v[i][j] - dt * (0.25 * ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j])) - ( )
		}
	}
step += 1;
printf("Step = %d\n", step);
if(step > 4)
{
	break;
}
}
}
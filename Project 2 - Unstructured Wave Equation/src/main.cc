// Included headers
#include "Mesh_Read.h"
#include "field_functions.h"

/********************************************************************
// Unstructured Mesh Programming Assignment Main File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

int main()
{
// Calculate runtime
time_t start, end;
time(&start);

// Input data
// std::string meshname = "Face-Cell/Unit_square.mesh";
std::string meshname = "Face-Cell/mech511-square-medium.mesh"; // Enter mesh file name
double dt_max = 8.04e-3; // Enter maximum stable timestep here
double max_time = 2; // Enter maximum time here
double t0 = 0.0; // Enter initial time here 
double CFL = 1.0; // Enter CFL number here. Timestep will be CFL*dt_max
int num_fps = 20; // How many solutions to save
int num_sol = max_time / (dt_max * num_fps); // How often do we want to save a solution file? Every num_sol timesteps, a solution will be saved
double dt = CFL * dt_max; // Timestep dt
int num_timestep = max_time / dt; // Number of timesteps before solution finishes
/*
Below, select the time scheme
1 - Euler
2 - Adams Bashforth
3 - Midpoint
4 - RK2 
*/
int timescheme = 4; 

// Read mesh
Mesh mesh = read_mesh(meshname);

// Initialize velocity vector, and initial conditions for temperature
std::vector<double> Temp = std::vector<double>(mesh.iNCell); // Vector to store temperature at each cell
std::vector<double> flux_curr = std::vector<double>(mesh.iNCell); // Store flux integral for each cell centroid for the current timestep
std::vector<double> flux_prev = std::vector<double>(mesh.iNCell); // Store flux integral for each cell centroid for the previous timestep
std::array<std::vector<double>, 2> vel = {std::vector<double>(mesh.iNEdge), std::vector<double>(mesh.iNEdge)}; // Store velocity at each edge midpoint
initial_condition(mesh, Temp, vel);

// Calculate upwind cells and dot products (This doesnt change since velocity is constant)
calc_upwind(mesh, vel);

// Initialize vectors to store cell gradients, and edge fluxes
std::array<std::vector<double>, 2> Cell_Grad = {std::vector<double>(mesh.iNCell), std::vector<double>(mesh.iNCell)};
std::vector<double> Edge_flux = std::vector<double>(mesh.iNEdge); // Store flux at each edge

// Initialize more variables
int itime = 0; // Number of timesteps passed
double time_curr = 0.0; // Current time
double InitNorm = calc_L2Norm(Temp); // Save initial L2Norm
double finalNorm; // Store final norm
double currNorm; // Store norm of current cell
double normratio; // Store residual norm ratio

// Store initial and final temperature to compare later
std::vector<double> Temp_init = std::vector<double>(mesh.iNCell);
std::vector<double> Temp_final = std::vector<double>(mesh.iNCell);
for(int i = 0; i < mesh.iNCell; i++)
{
	Temp_init[i] = Temp[i];
}

// Begin timestepping
while(time_curr < max_time)
{
	// Initialize all values of temperature at the previous step to be the values of the current timestep
	for(int i = 0; i < mesh.iNCell; i++)
	{
		flux_prev[i] = flux_curr[i];
	}

	// Calculate gradient at each cell
	calc_grad(mesh, Temp, Cell_Grad);
	// Calculate flux at each edge
	calc_flux(mesh, Temp, Cell_Grad, Edge_flux);
	// Calculate flux integral for each cell
	calc_flux_integral(mesh, Edge_flux, flux_curr);
	// Step time forward
	switch(timescheme)
	{
	case 1: 
		euler(mesh, Temp, flux_curr, dt);
		break;
	case 2: 
		adams_bashforth(mesh, Temp, flux_curr, flux_prev, dt, itime);
		break;
	case 3:
		midpoint(mesh, Temp, flux_curr, dt);
		break;
	case 4:
		RK2(mesh, Temp, flux_curr, dt);
		break;
	default:
		printf("Invalid time scheme, please check.\n");
		printf("Valid options:\n1. Euler\n2. Adams-Bashforth\n3. Midpoint\n4. RK2\n");
		exit(0);
	}
	// Save VTK file every num_sol timesteps
	if(itime % num_sol == 0)
	{
		printf("VTK File Saved\n");
		std::string fileName = "Results/Timestep-" + std::to_string(itime) + ".vtk";
		save_VTK(fileName, mesh, Temp);
	}

	if(itime == num_timestep)
	{
		printf("*************************************\nSaving Final Data\n*************************************\n");
		finalNorm = calc_L2Norm(Temp);
		// Save temperature vector after 1 period to compare with initial solution
		for(int i = 0; i < mesh.iNCell; i++)
		{
			Temp_final[i] = Temp[i];
		}
	}

	// Increment timestep and time
	itime += 1;
	time_curr += dt;
	currNorm = calc_L2Norm(Temp);
	normratio = currNorm / InitNorm;
	if(normratio > 5)
	{
		printf("************************************************\n");
		printf("             SOLUTION EXPLODING                 \n");
		printf("************************************************\n");
	}
	printf("Timestep: %i Time: %f dt: %4.2e  Timescheme: %i\nNorm: %f Residual Norm Ratio: %f\n", itime, time_curr, dt,  timescheme, currNorm, normratio);
}

switch(timescheme)
{
	case 1: 
		printf("Solution finished after %i timesteps using Euler time marching scheme\n", itime);
		break;
	case 2: 
		printf("Solution finished after %i timesteps using Adams-Bashforth time marching scheme\n", itime);
		break;
	case 3:
		printf("Solution finished after %i timesteps using Midpoint time marching scheme\n", itime);
		break;
	case 4:
		printf("Solution finished after %i timesteps using RK2 time marching scheme\n", itime);
		break;
}

// Calculate error
std::vector<double> Temp_error = std::vector<double>(mesh.iNCell);
for(int i = 0; i < mesh.iNCell; i++)
{
	Temp_error[i] = Temp_final[i] - Temp_init[i];
}

double error = calc_L2Norm(Temp_error);
printf("Final L2Norm Error is %14.12e\n", error);

time(&end);
double time_taken = double(end - start); 
printf("Code succesfully run in %.3f seconds\n", time_taken);
return 0;
}

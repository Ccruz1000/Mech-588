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
std::string meshname = "Face-Cell/mech511-square-verycoarse.mesh"; // Enter mesh file name
double dt_max = 2.5e-3; // Enter maximum stable timestep here
double max_time = 2; // Enter maximum time here
double t0 = 0.0; // Enter initial time here 
double CFL = 1.0; // Enter CFL number here. Timestep will be CFL*dt_max
int num_sol = 30; // How often do we want to save a solution file? Every num_sol timesteps, a solution will be saved
double dt = CFL * dt_max; // Timestep dt

// Read mesh
Mesh mesh = read_mesh(meshname);

// Vector to store temperature at each cell
std::vector<double> Temp = std::vector<double>(mesh.iNCell); 

// Initialize velocity vector, and initial conditions for temperature
std::vector<double> flux_curr = std::vector<double>(mesh.iNCell); // Store temperature at each cell centroid for the current timestep
std::vector<double> flux_prev = std::vector<double>(mesh.iNCell); // Store temperature at each cell centroid for the previous timestep
std::array<std::vector<double>, 2> vel = {std::vector<double>(mesh.iNEdge), std::vector<double>(mesh.iNEdge)}; // Store velocity at each edge midpoint
initial_condition(mesh, Temp, vel);

// Calculate upwind cells and dot products (This doesnt change since velocity is constant)
calc_upwind(mesh, vel);

// Initialize vectors to store cell gradients, and edge fluxes
std::array<std::vector<double>, 2> Cell_Grad = {std::vector<double>(mesh.iNCell), std::vector<double>(mesh.iNCell)};
std::vector<double> Edge_flux = std::vector<double>(mesh.iNEdge); // Store flux at each edge


int itime = 0; // Number of timesteps passed
double time_curr = 0.0; // Current time

double InitNorm = calc_L2Norm(flux_curr); // Save initial L2Norm

while(time_curr < max_time)
{
	// Initialize all values of temperature at the previous step to be the values of the current timestep
	for(int i = 0; i < mesh.iNCell; i++)
	{
		flux_prev[i] = flux_curr[i];
	}

	// Calculate gradient at each cell
	calc_grad(mesh, flux_curr, Cell_Grad);
	// Calculate flux at each edge
	calc_flux(mesh, flux_curr, Cell_Grad, Edge_flux);

	// Timestep using Adams-Bashforth second order
	itime += 1;
	time_curr += dt;
	printf("Time at timestep %i is %f\n", itime, time_curr);
}
// 

// // Read mesh 
// std::string verycoarse = "Face-Cell/mech511-square-verycoarse.mesh";
// std::string coarse = "Face-Cell/mech511-square-coarse.mesh";
// std::string medium = "Face-Cell/mech511-square-medium.mesh";
// std::string fine = "Face-Cell/mech511-square-fine.mesh";
// std::string veryfine = "Face-Cell/mech511-square-veryfine.mesh";
// std::string analytical = "Face-Cell/analytical.mesh";

// Mesh mesh = read_mesh(analytical);

// std::vector<double> temp_cent = std::vector<double>(mesh.iNCell); // Store temperatures at cell centroids
// std::array<std::vector<double>, 2> vel = {std::vector<double>(mesh.iNEdge), std::vector<double>(mesh.iNEdge)}; // Store velocity u and v at edge midpoints
// std::vector<double> Edge_flux = std::vector<double>(mesh.iNEdge); // Store flux at each edge

// initial_condition(mesh, temp_cent, vel);
// std::string file = "Results/Test.vtk";
// save_VTK(file, mesh, temp_cent);

// std::array<std::vector<double>, 2> Cell_Grad = {std::vector<double>(mesh.iNCell), std::vector<double>(mesh.iNCell)};
// calc_grad(mesh, temp_cent, Cell_Grad);
// calc_upwind(mesh, vel);
// calc_flux(mesh, temp_cent, Cell_Grad, Edge_flux);

time(&end);
double time_taken = double(end - start); 
printf("Code succesfully run in %.3f seconds\n", time_taken);
return 0;
}
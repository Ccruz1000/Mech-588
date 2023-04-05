// Included header files
#include "Roe.h"

/********************************************************************
// Compressible Flow Programming Assignment Main File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

int main(int argc, char* argv[])
{
	// Declare run parameters
	double CFL, xbegin, xend, dt, dx, gamma, t_final;
	int nx;

	// We will store the initial condition values in an array. The first 3 entries will be for the left side
	// and will represent rho, u, and P in that order. The next 3 entries will be for the right side, and will 
	// represent rho, u, and P in that order

	double init_cond[6] = {0};

	// Read input file, or provide default values if none given
	assign_input(argc, argv, CFL, nx, xbegin, xend, dt, gamma, init_cond, t_final);

	// printf("\nCFL: %f xbegin: %f xend: %f dt: %f gamma: %f nx: %d\n", CFL, xbegin, xend, dt, gamma, nx);
	// printf("\nInitial conditions:\nLeft: rho - %f u - %f P - %f\nRight: rho - %f u - %f P - %f\n", init_cond[0], init_cond[1], init_cond[2], init_cond[3], init_cond[4], init_cond[5]);

	// Calculate the spatial step
	dx = (xend - xbegin) / (nx - 1);

	// Create array to store primitives. Prim[0][i] is rho, Prim[1][i] is velocity, Prim[2][i] is pressure, and Prim[3][i] is energy
	std::array<std::vector<double>, 4> Prim = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};

	// Assign initial condition values
	initial_conditions(Prim, nx, init_cond, gamma);

	// Save initial condition results
	std::string fileName = "Results/Init_cond.txt";
	save_results(fileName, Prim, dx, nx);

	// Create and populate array to store conserved variables. Calculate speed of sound in each cell at the same time
	std::array<std::vector<double>, 3> U = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};  // Solution vector
	std::array<std::vector<double>, 3> F = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)}; // Flux vector
	// Populate solution and sound speed vector using primitives
	int num_save = 100;
	for(int i = 0; i < nx; i++)
	{
		U[0][i] = Prim[0][i];  // rho
		U[1][i] = Prim[0][i] * Prim[1][i];  // rho * u
		U[2][i] = Prim[3][i]; // Et
	}

	// Calculate initial flux
	Roe_flux(U, Prim, F, dx, gamma, nx);

	// Begin timestepping
	double t_curr = 0;
	int itime = 0;
	while(t_curr < t_final)
	{
		Roe_flux(U, Prim, F, dx, gamma, nx);
		forward_euler(U, F, dt, dx, nx);

		itime += 1;
		t_curr += dt;
		std::string inter_results;
		inter_results = "Results/inter_" + std::to_string(itime) + ".txt";
		if (itime % num_save == 0)
		{
			save_results(inter_results, Prim, dx, nx);
		} 
		printf("Timestep: %i Time: %f\n", itime, t_curr);
	}

	std::string final_results = "Results/Final_cond.txt";
	save_results(final_results, Prim, dx, nx);


	return 0;
}
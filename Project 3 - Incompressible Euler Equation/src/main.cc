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
	double CFL, xbegin, xend, dt, dx, gamma;
	int nx;

	// We will store the initial condition values in an array. The first 3 entries will be for the left side
	// and will represent rho, u, and P in that order. The next 3 entries will be for the right side, and will 
	// represent rho, u, and P in that order

	double init_cond[6] = {0};

	// Read input file, or provide default values if none given
	assign_input(argc, argv, CFL, nx, xbegin, xend, dt, gamma, init_cond);

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
	std::array<std::vector<double>, 3> U = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	std::vector<double> c = std::vector<double>(nx);
	for(int i = 0; i < nx; i++)
	{
		U[0][i] = Prim[0][i];
		U[1][i] = Prim[0][i] * Prim[1][i];
		U[2][i] = Prim[3][i] * Prim[0][i];
		c[i] = std::sqrt(gamma * Prim[2][i] / Prim[0][i]);
	}
	Roe_flux(U, Prim, dx, gamma, c, nx);


	return 0;
}
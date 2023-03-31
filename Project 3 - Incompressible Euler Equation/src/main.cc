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

	printf("\nCFL: %f xbegin: %f xend: %f dt: %f gamma: %f nx: %d\n", CFL, xbegin, xend, dt, gamma, nx);
	printf("\nInitial conditions:\nLeft: rho - %f u - %f P - %f\nRight: rho - %f u - %f P - %f\n", init_cond[0], init_cond[1], init_cond[2], init_cond[3], init_cond[4], init_cond[5]);

	// Calculate the spatial step
	dx = (xend - xbegin) / (nx - 1);

	// Assign vectors to store parameters
	std::vector<double> rho = std::vector<double>(nx);  // Store density
	std::vector<double> u = std::vector<double>(nx);  // Store velocity
	std::vector<double> P = std::vector<double>(nx);  // Store pressure
	
	// Assign initial condition values
	initial_conditions(rho, u, P, nx, init_cond);

	// Save initial condition results
	std::string fileName = "Results/Init_cond.txt";
	save_results(fileName, rho, u, P, dx, nx);



	return 0;
}
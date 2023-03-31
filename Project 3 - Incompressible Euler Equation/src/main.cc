// Included header files
#include "Roe.h"

/********************************************************************
// Compressible Flow Programming Assignment Main File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

int main(int argc, char* argv[])
{
	// Run parameters
	double CFL, xbegin, xend, dt, gamma;
	int nx; // Store number of points in mesh

	// Check if input file was supplied
	if(argc < 2)
	{
		printf("No input file was provided. Please check.\nDefault parameters will be assigned to all input parameters.\n");
		
		std::string input_file = "Null";
		read_input(input_file, CFL, nx, xbegin, xend, dt, gamma);
	}

	// Let user know that input was successfully provided
	else 
	{
		// Convert the inputs file from an array of char to a string
		std::string input_file = argv[1]; // = argv[1].c_str();

		// Inform user that the input file was successfully provided, and let them know which one was provided.
		std::cout << "Input file provided successfully" << std::endl;
		std::cout << "File given: " << input_file << std::endl;
		read_input(input_file, CFL, nx, xbegin, xend, dt, gamma);
	}

	printf("CFL: %f xbegin: %f xend: %f dt: %f gamma: %f nx: %d\n", CFL, xbegin, xend, dt, gamma, nx);

	return 0;
}
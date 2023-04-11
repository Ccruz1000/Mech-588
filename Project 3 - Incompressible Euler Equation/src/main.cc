// Included header files
#include "functions.h"
#include "Steger.h"
#include "Roe.h"

/********************************************************************
// Compressible Flow Programming Assignment Main File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

int main(int argc, char* argv[])
{	
	// Calculate runtime
	clock_t tstart = clock();

	// Run Roe
	run_sim(argc, argv, 0);
	
	// Run Steger-Warming
	run_sim(argc, argv, 1);

	// Run Lax-Friedrich
	run_sim(argc, argv, 2);

	// Let user know runtime
	printf("\nCode succesfully ran in %6.3e seconds\n", (double)(clock() - tstart)/CLOCKS_PER_SEC);
	return 0;
}
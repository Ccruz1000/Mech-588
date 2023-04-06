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
{	// Run Roe
	// Need to refactor Roe before I can try running it
	// run_sim(argc, argv, 0);
	
	// Run Steger-Warming
	run_sim(argc, argv, 1);


	return 0;
}
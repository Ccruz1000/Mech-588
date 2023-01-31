#include <iostream>
#include <string>
#include <fstream>

#ifndef MESHRD_H

#define MESHRD_H

// Taken from base code provided in Mech 587
class Operations{
public:
	void alloc1D(double **V, const unsigned long &N, bool isInit = true, double def = 0.0) const;
	void deAlloc1D(double **V, const unsigned long &N) const;
	
};

class Solution
{
private:
	Operations O;
	

public:
	// Members
	double *X, *Y, *u, *v, *T; // Pointers to arrays for coordinates, and solution vectors
	unsigned long nx, ny; // Value to store number of X and Y
	unsigned long N; // Store total number of points
	// Constructors
	Solution();
	Solution(unsigned long Nx, unsigned long Ny);
	//Solution(unsigned long Nx, unsigned long Ny, double *X, double *Y, double *u, double *v, double *T);
	//void storeVTKSolution;
	

};


#endif
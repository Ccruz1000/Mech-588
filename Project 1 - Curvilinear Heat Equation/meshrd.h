#include <iostream>
#include <string>
#include <fstream>
#include <vector>

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

public:
	// Member Variables
	std::vector<double> X , Y, u, v, T;  // Vectors to store data
	unsigned long nx, ny; // Value to store number of X and Y points
	unsigned long N; // Store total number of points
	// Constructors
	Solution(); // Null constructor, initializes vectors as empty
	Solution(unsigned long iNx, unsigned long iNy); // Void constructor, initializes with all 0's
	Solution(unsigned long iNx, unsigned long iNy, std::vector<double> iX, std::vector<double> iY, 
		     std::vector<double> iu, std::vector<double> iv, std::vector<double> iT); // Constructor, initializes with pre defined vectors
	// Member Functions
	void printmesh();
	//void storeVTKSolution(const Solution &s, const char *fileName);
};


#endif
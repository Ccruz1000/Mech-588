#include "meshrd.h"


// Taken from Mech 587 base code
void Operations::alloc1D(double **V, const unsigned long &N, bool isInit, double def) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}

	*V = new double[N];

	if(isInit)
		for(unsigned long i = 0; i < N; i++)
			(*V)[i] = 0.0;
}

void Operations::deAlloc1D(double **V, const unsigned long &N) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}
}


// Base case - Create empty solution
Solution::Solution() : X(nullptr), Y(nullptr), u(nullptr), v(nullptr), T(nullptr), nx(0), ny(0), N(0) 
{}


// Constructor to initialize solution struct with 0 values
Solution::Solution(unsigned long iNx, unsigned long iNy) : nx(iNx), ny(iNy)
{
	N = nx * ny; // Total number of grid points

	// Allocate memory for solution variables
	O.alloc1D(&X, N);
	O.alloc1D(&Y, N);
	O.alloc1D(&u, N);
	O.alloc1D(&v, N);
	O.alloc1D(&T, N);

	// Initialize variables values to 0
	for(int i = 0; i < N; i++)
	{
		X[i] = 0.0;
		Y[i] = 0.0;
		u[i] = 0.0;
		v[i] = 0.0;
		T[i] = 0.0;
	}
}




int main()
{
	// Different mesh files provided
	std::string coursemesh = "mesh-8x32.vel";
	std::string medmesh = "mesh-16x64.vel";
	std::string finemesh = "mesh-32x128.vel";


	// Read individual lines without whitespace
	double Nx, Ny; // Number of X and Y points in mesh 

	std::ifstream mesh(coursemesh); // Open mesh file for reading  

	int cntr = 0; // counter to check which line we are on

	int ifile, jfile; // i and j to be read from the file
	double x, y, u, v; // coordinates and velocities read from the file 

	if(mesh.is_open()) // Check that file is opened
	{
		if(cntr == 0) // If reading first line, read for number of points in X and Y
		{
			mesh >> Nx >> Ny;
			cntr += 1;
			std::cout << Nx << " Nx\n" << Ny << " Ny\n";
		}
		mesh >> ifile >> jfile >> x >> y >> u >> v;
		std::cout << "i: " << ifile << " j: " << jfile << " x: " << x << " y: " << y << " u: " << u << " v: " << v << std::endl;

	}

	// Test to see if void constructor works properly
	unsigned long nx, ny;
	Solution test(10, 10);
	nx = test.nx;
	ny = test.ny;
	std::cout << std::endl << std::endl << nx << " Nx " << ny << " Ny\n";

	// TODO Figure out how to read arrays from Solution variable that we currently have
	// TODO Allow entry of custom arrays into solution class
	// TODO Learn how to store solution class into VTK file
	// TODO Literally the rest of the entire project

	// double X[10] = {test.X};
	// for(int k = 0; k < 10; k++)
	// {
	// 	std::cout << X[k] << " X\n";
	// }
	return 0;
}


#include "meshrd.h"


// Taken from Mech 587 base code 
// Allocate memory for dynamically sized arrays
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

// Deallocate memory for dynamically sized arrays
void Operations::deAlloc1D(double **V, const unsigned long &N) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}
}


// Base case - Create empty solution
Solution::Solution() : X(nullptr), Y(nullptr), u(nullptr), v(nullptr), T(nullptr), nx(0), ny(0), N(0) 
{

}

// Constructor to initialize solution struct with 0 values
Solution::Solution(unsigned long iNx, unsigned long iNy) : X(nullptr), Y(nullptr), u(nullptr), v(nullptr), T(nullptr), nx(iNx), ny(iNy)
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

// Constructor to initialize solution struct using vectors that already exist
Solution::Solution(unsigned long iNx, unsigned long iNy, double iX, double iY, double iu, double iv, double iT) :
                   X(*iX), Y(*iY), u(*iu), v(*iv), T(*iT), nx(iNx), ny(iNy)
{
	N = nx * ny;
}


// Test functions
void test_void_constructor()
{
	// Test to see if void constructor works properly
	unsigned long nx, ny, N;
	double *X;
	Solution test(4, 3);
	nx = test.nx;
	ny = test.ny;
	N = test.N;
	X = test.X;
	X[5] = 12.42;
	int i = 1;
	std::cout << std::endl << std::endl << nx << " Nx " << ny << " Ny " << N << " N\n";
	for(int k = 0; k < N; k++)
	{
		std::cout << "i: " << i << " X: " << X[k] << std::endl;
		i += 1;
	}
}

void test_mesh_read(std::string meshname)
{
	std::ifstream mesh(meshname); // Open mesh file for reading  
	// Read individual lines without whitespace
	long unsigned int Nx, Ny, numpoints; // Number of X and Y points in mesh 
	int ifile, jfile; // i and j to be read from the file
	//double x, y, u, v;
	double x_store, y_store, u_store, v_store; // Store values to input into array later
	//std::ifstream mesh(falsemesh); // Test that catch for improper mesh works
	if(mesh.is_open()) // Check that file is opened
	{
		// Read Nx and Ny values from first line 
		mesh >> Nx >> Ny;
		numpoints = Nx * Ny; 
		std::cout << numpoints << std::endl;
		std::cout << Nx << " Nx\n" << Ny << " Ny\n";
		// Create arrays to store x, y, u and v from mesh
		double x[numpoints], y[numpoints], u[numpoints], v[numpoints];
		int cntr = 0; // Counter to index position in the array
		// Read through other lines
		while(mesh >> ifile >> jfile >> x_store >> y_store >> u_store >> v_store)
		{
			// Store read values into mesh
			x[cntr] = x_store;
			y[cntr] = y_store;
			u[cntr] = u_store;
			v[cntr] = v_store;
			cntr += 1; 
		}
		// Print out values read from mesh to check
		for(int i = 0; i < numpoints; i++)
		{	
			if(i % Ny == 0)
			{
				std::cout << std::endl;
			}
			std::cout << "X: " << x[i] << " Y: " << y[i] << " u: " << u[i] << " v: " << v[i] << std::endl;
		}
	}

	// Let user know if mesh file isn't found 
	else
	{
		std::cout << "Mesh file not found, please check file name\n";
	}
}


int main()
{
	// Different mesh files provided
	std::string coursemesh = "mesh-8x32.vel";
	std::string medmesh = "mesh-16x64.vel";
	std::string finemesh = "mesh-32x128.vel";
	std::string falsemesh = "mesh-13x2.vel"; // False mesh name to make sure catch works

	//test_void_constructor();
	//test_mesh_read(medmesh);
	
	// unsigned long nx, ny;
	// nx = 4;
	// ny = 4;
	// double *x[4] = {1.0, 2.0, 3.0, 4.0};
	// double *y[4] = {1.0, 2.0, 3.0, 4.0};
	// double *u[4] = {1.0, 2.0, 3.0, 4.0};
	// double *v[4] = {1.0, 2.0, 3.0, 4.0};
	// double *T[4] = {1.0, 2.0, 3.0, 4.0};
	// Solution test(nx, ny, x, y, u, v, T);





	// TODO Allow entry of custom arrays into solution class
	// TODO Learn how to store solution class into VTK file
	// TODO Literally the rest of the entire project

	return 0;
}


	//Solution test(nx, ny, x, y, u, v, T);

#include "meshrd.h"

// Solution class constructors and member classes

// Base case - Create empty solution
Solution::Solution() : nx(0), ny(0), N(0), X(0), Y(0), u(0), v(0), T(0) 
{}

// Constructor to initialize solution struct with 0 values
Solution::Solution(unsigned long iNx, unsigned long iNy) : nx(iNx), ny(iNy), X(iNx * iNy), Y(iNx * iNy), u(iNx * iNy), v(iNx * iNy), T(iNx * iNy)
{
	unsigned long N = nx * ny;
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

// Constructor to initialize solution struct with pre defined arrays
Solution::Solution(unsigned long iNx, unsigned long iNy, std::vector<double> iX, std::vector<double> iY, std::vector<double> iu, std::vector<double> iv, std::vector<double> iT) 
                  : nx(iNx), ny(iNy), X(iX), Y(iY), u(iu), v(iv), T(iT)
{}

// Member functions
// Print out the values read from the mesh
void Solution::printmesh()
{
	for(int i = 0; i < (nx * ny); i++)
	{
		if(i % ny == 0)
		{
			std::cout << std::endl;
		}
		std::cout << "X: " << X[i] << " Y: " << Y[i] << " u: " << u[i] << " v: " << v[i] << " T: " << T[i] << std::endl;
	}
}


/******************************************************************************************
 ******************************************************************************************/

// Useful functions

// Read mesh file and store it as a solution
Solution readmesh(std::string meshname)
{
	std::ifstream mesh(meshname); // Open mesh file
	// Read first line to obtain Nx and Ny values
	long unsigned int nx, ny;
	mesh >> nx >> ny;
	Solution newmesh(nx, ny);
	
	if(mesh.is_open())
	{
		std::cout << "Mesh opened\n";
		//mesh >> nx >> ny;
		// Initialize void solution using number of points from mesh
		//Solution newmesh(nx, ny);
		// Initialize vectors to store data from solution
		double X, Y, u, v, T; // Temporary variable to store value in solution vector
		// Read through rest of lines to populate vectors
		int cntr = 0; // Counter to index position in array
		int i, j; // Store i and j from file
		// Loop through each line
		//while(mesh >> i >> j >> X >> Y >> u >> v)
		while(true)
		{
			mesh >> i >> j >> X >> Y >> u >> v;
			newmesh.X[cntr] = X;
			newmesh.Y[cntr] = Y;
			newmesh.u[cntr] = u;
			newmesh.v[cntr] = v;
			// X[cntr] = xtemp;
			// Y[cntr] = ytemp;
			// u[cntr] = utemp;
			// v[cntr] = vtemp;
			// T[cntr] = Ttemp;
			cntr += 1;
			// Leave loop once all lines of mesh have been read
			if(mesh.eof())
			{
				std::cout << "Last line reached, exiting\n";
				mesh.close();
				std::cout << "Mesh closed\n";
				break;
			}
		}
	}
	// Let user know if file is not read properly
	else
	{
		std::cout << "Mesh file not found\n";
	}
	std::cout << "Mesh gets returned in next line\n";
	return newmesh;
	std::cout << "Mesh returned in last line\n";
}


// void storeVTKSolution(const Solution &s, const char *fileName)
// {
// }

/******************************************************************************************
 ******************************************************************************************/

// Test functions

// Test null constructor to initialize empty solution
void test_null_constructor()
{
	unsigned long nx;
	std::vector<double> X;
	Solution test1;
	X = test1.X;
	nx = test1.nx;
	std::cout << "nx: " << nx << std::endl;
	std::cout << "X: " << X.size() << std::endl;
}

// Test void constructor to initialize solution of 0's
void test_void_constructor()
{
	// Test to see if void constructor works properly
	unsigned long nx, ny, N;
	std::vector<double> X, Y, u, v, T;
	Solution test(6, 1);
	nx = test.nx;
	ny = test.ny;
	N = test.N;
	X = test.X;
	Y = test.Y;
	u = test.u;
	v = test.v;
	T = test.T;
	X[0] = 12.42;
	Y[1] = 41.52;
	u[2] = 2.04;
	v[3] = 1.95;
	T[4] = -1.24;
	int i = 0;
	std::cout << std::endl << nx << " Nx " << ny << " Ny " << N << " N\n";
	for(int k = 0; k < N; k++)
	{
		std::cout << "i: " << i << " X: " << X[k] << " Y: " << Y[k] << " u: " << u[k] << " v: " << v[k] << " T: " << T[k] << std::endl;
		i += 1;
	}
}

// Test constructor to initialize solution with vectors that exist
void test_constructor()
{
	unsigned long nx, ny;
	nx = 4;
	ny = 1;
	std::vector<double> X(nx * ny);
	std::vector<double> Y(nx * ny);
	std::vector<double> u(nx * ny);
	std::vector<double> v(nx * ny);
	std::vector<double> T(nx * ny);
	X = {1, 1, 1, 1};
	Y = {2, 2, 2, 2};
	u = {3, 3, 3, 3};
	v = {4, 4, 4, 4};
	T = {5, 5, 5, 5};
	Solution test(nx, ny, X, Y, u, v, T);
	std::vector<double> Xtest, Ytest, utest, vtest, Ttest;
	unsigned long nxtest, nytest;
	nxtest = test.nx;
	nytest = test.ny;
	Xtest = test.X;
	Ytest = test.Y;
	utest = test.u;
	vtest = test.v;
	Ttest = test.T;

	std::cout << "Nx: " << nxtest << " Ny: " << nytest << std::endl;
	for(int i = 0; i < (nx * ny); i++)
	{
		std::cout << "i: " << i + 1 << " X: " << Xtest[i] << " Y: " << Ytest[i] << " u: " << utest[i] << " v: " << vtest[i] << " T: " << Ttest[i] << std::endl;;  
	}
}

// Test that reading mesh occurs properly
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

	Solution mesh = readmesh(coursemesh);
	mesh.printmesh();

	// const char vtkname[50] = "file.vtk";
	std::string vtkname = "file.vtk";
	std::ofstream vtkfile(vtkname);
	unsigned long i;
	vtkfile << "# vtk DataFile Version 2.0\n";
	vtkfile << "TITLE = \"Quad data\"\n";
	vtkfile << "ASCII\n";
	vtkfile << "DATASET STRUCTURED_GRID\n";
	vtkfile << "Dimensions " << mesh.nx << " " << mesh.ny << " " << 1 << std::endl;
	vtkfile << "POINTS " << mesh.nx * mesh.ny << " FLOAT\n";
	for(i = 0; i < mesh.nx * mesh.ny; i++)
	{
		vtkfile << mesh.X[i] << " " << mesh.Y[i] << " " << 0.0 << std::endl;
	}
	vtkfile.close();
	// unsigned long i, j;
	// fprintf(vtkfile,"# vtk DataFile Version 2.0\n");
	// fprintf(vtkfile,"TITLE = \"Quad data\"\n");
	// fprintf(vtkfile,"ASCII\n");
	// fprintf(vtkfile,"DATASET STRUCTURED_GRID\n");
	// fprintf(vtkfile,"DIMENSIONS %lu %lu %d\n", mesh.nx, mesh.ny, 1);
	// fprintf(vtkfile,"POINTS %lu FLOAT\n",mesh.nx*mesh.ny);
	// for(i = 0; i < mesh.nx * mesh.ny; i++)
	// {
	// 	fprintf(vtkfile, "%f %f %f\n", mesh.X[i], mesh.Y[i], 0.0);
	// }
	// vtkfile.close();


		




	// TODO Learn how to store solution class into VTK file
	// TODO Literally the rest of the entire project

	return 0;
}


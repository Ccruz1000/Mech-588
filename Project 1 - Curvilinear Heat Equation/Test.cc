#include "Fcurv.h"
#include "Test.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Test File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/


// Test void constructor to initialize solution with vectors that exist 
void test_void_constructor()
{
	unsigned long nx = 2;
	unsigned long ny = 2;
	Solution test(nx, ny);
	std::cout << test.N << std::endl;
	for (int i = 0; i < test.nx; i++)
		for(int j = 0; j < test.ny; j++)
		{
			printf("i: %u j: %u X(i, j): %f Y(i,j): %f u(i, j):  %f v(i, j): %f T(i, j): %f\n", i, j, test.X(i, j), test.Y(i, j), test.u(i, j), test.v(i, j), test.T(i, j));
		}
}

// Test constructor to initialize solution with vectors that exist
void test_constructor()
{
	unsigned long nx = 2;
	unsigned long ny = 2;
	double cntr = 1.0;
	Vector iX(nx, ny);
	Vector iY(nx, ny);
	Vector iu(nx, ny);
	Vector iv(nx, ny);
	Vector iT(nx, ny);

	for (unsigned long i = 0; i < nx; i++)
		for(unsigned long j = 0; j < ny; j++)
		{
			iX(i, j) = cntr;
			iY(i, j) = cntr;
			iu(i, j) = cntr;
			iv(i, j) = cntr;
			iT(i, j) = cntr;
			cntr += 1.0;
		}
	Solution test(nx, ny, iX, iY, iu, iv, iT);

	for (unsigned int i = 0; i < nx; i++)
		for(unsigned int j = 0; j < ny; j++)
		{
			printf("i: %u, j: %u, X(i, j): %f, Y(i, j) %f, u(i, j): %f, v(i, j): %f, T(i, j): %f\n", i, j, test.X(i, j), test.Y(i, j), test.u(i, j), test.v(i, j), test.T(i, j));
		}
}

void test_mesh_read(std::string meshname)
{
	std::ifstream mesh(meshname); // Open mesh file for reading  
	// Read individual lines without whitespace
	long unsigned int iNx, iNy, numpoints; // Number of X and Y points in mesh 
	int ifile, jfile; // i and j to be read from the file
	double x_store, y_store, u_store, v_store; // Store values to input into array later
	if(mesh.is_open()) // Check that file is opened
	{
		// Read Nx and Ny values from first line 
		mesh >> iNx >> iNy;
		numpoints = iNx * iNy; 
		std::cout << numpoints << std::endl;
		std::cout << iNx << " Nx\n" << iNy << " Ny\n";
		// Create arrays to store x, y, u and v from mesh
		double x[numpoints], y[numpoints], u[numpoints], v[numpoints], T[numpoints];
		int cntr = 0; // Counter to index position in the array
		double xcheck, ycheck, ucheck, vcheck; // Check 4 random pre selected points
		// Read through other lines
		while(mesh >> ifile >> jfile >> x_store >> y_store >> u_store >> v_store)
		{
			// Store read values into mesh
			x[cntr] = x_store;
			y[cntr] = y_store;
			u[cntr] = u_store;
			v[cntr] = v_store;
			if(ifile == 27 and jfile == 6)
			{
				xcheck = x_store;
			}
			if(ifile == 31 and jfile == 5)
			{
				ycheck = y_store;
			}
			if(ifile == 22 and jfile == 7)
			{
				ucheck = u_store;
			}
			if(ifile == 29 and jfile == 4)
			{
				vcheck = v_store;
			}
			cntr += 1; 
		}


	// Create empty solution to assign values from 1D double array above
	Solution test(iNx, iNy);
	int index = 0;
	for(int i = 0; i < iNx; i++)
		for(int j = 0; j < iNy; j++)
		{
			test.X(i, j) = x[i * iNy + j];
			test.Y(i, j) = y[i * iNy + j];
			test.u(i, j) = u[i * iNy + j];
			test.v(i, j) = v[i * iNy + j];
		}

	// Print out solution to check that mesh is read right
	for(int i = 0; i < iNx; i++)
		for(int j = 0; j < iNy; j++)
		{
			printf("X: %f, Y: %f, u: %f, v: %f\n", test.X(i, j), test.Y(i, j), test.u(i, j), test.v(i, j));
			if(j == iNy - 1)
			{
				printf("\n");
			}
		}

	std::cout << "Comparison tests\n\n";

	// Some comparisons to ensure validation even better
	bool xmesh6_27 = (test.X(27, 6) == xcheck);
	bool ymesh5_31 = (test.Y(31, 5) == ycheck);
	bool umesh7_22 = (test.u(22, 7) == ucheck);
	bool vmesh4_29 = (test.v(29, 4) == vcheck);
	std::cout << "Check 6_27: " << xmesh6_27 << std::endl;
	std::cout << "Check 5_31: " << ymesh5_31 << std::endl;
	std::cout << "Check 7_22: " << umesh7_22 << std::endl;
	std::cout << "Check 4_29: " << vmesh4_29 << std::endl;
	}


	// Let user know if mesh file isn't found 
	else
	{
		std::cout << "Mesh file not found, please check file name\n";
	}
	mesh.close();
}


int main()
{
	// Different mesh files provided
	std::string coursemesh = "mesh-8x32.vel";
	std::string medmesh = "mesh-16x64.vel";
	std::string finemesh = "mesh-32x128.vel";
	std::string falsemesh = "mesh-13x2.vel"; // False mesh name to make sure catch works
	std::string fileName = "file.vtk"; // filename to store vtk file

	// Solution test = readmesh(coursemesh);
	// storeVTKSolution(test, fileName);
	// test_void_constructor();
	/// test_constructor();
	// test_mesh_read(coursemesh);

	// TODO Literally the rest of the entire project

	return 0;
}
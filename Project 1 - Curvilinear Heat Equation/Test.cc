#include "Fcurv.h"
#include "Test.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Test File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

/********************************************************************
// Test functions - Used to validate creation of solution class
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

// Test VTK output
void test_VTK_output(Solution &s)
{
	std::string fileName = "Test.vtk";
	const char * fileChar = fileName.c_str();
	FILE *vtkFile;
	vtkFile = fopen(fileChar, "w");
	unsigned long i, j;
	fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
	fprintf(vtkFile,"TITLE = \"Quad data\"\n");
	fprintf(vtkFile,"ASCII\n");
	fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
	fprintf(vtkFile,"DIMENSIONS %lu %lu %d\n", s.ny, s.nx, 1);
	fprintf(vtkFile,"POINTS %lu FLOAT\n",s.N);
	for(int i = 0; i < s.nx; i++)
		for(int j = 0; j < s.ny; j++)
		{
			fprintf(vtkFile, "%f %f %f\n", s.X(i, j), s.Y(i, j), 0.0);
		}
		// Change below to CELL_DATA
	fprintf(vtkFile,"CELL_DATA %lu\n", (s.nx - 1) * (s.ny - 1));
	fprintf(vtkFile,"SCALARS temperature FLOAT 1\n");
	fprintf(vtkFile,"LOOKUP_TABLE default\n");
	for(int i = 0; i < s.nx - 1; i++)
		for(int j = 0; j < s.ny - 1; j++)
		{
			double PI = 4.0*atan(1.0);
			s.T(i+1,j+1) = sin(PI*s.X(i+1,j+1))*sin(PI*s.Y(i+1,j+1));
			s.T(i,j+1) = sin(PI*s.X(i,j+1))*sin(PI*s.Y(i,j+1));
			s.T(i,j) = sin(PI*s.X(i,j))*sin(PI*s.Y(i,j));
			s.T(i+1,j) = sin(PI*s.X(i+1,j))*sin(PI*s.Y(i+1,j));

			double Tavg = 0.25*(s.T(i+1,j+1) + s.T(i,j) + s.T(i,j+1)+ s.T(i+1,j));
			fprintf(vtkFile, "%lf\n", Tavg);
		}
	// fprintf(vtkFile,"VECTORS velocity FLOAT\n");
	// for(int i = 0; i < s.nx; i++)
	// 	for(int j = 0; j < s.ny; j++)
	// 	{
	// 		fprintf(vtkFile, "%lf %lf %lf\n", s.u(i, j), s.v(i, j), 0.0);
	// 	}
	fclose(vtkFile);
}

/********************************************************************
// Test functions - Used to validate mesh metric calculations
********************************************************************/
// i+1/2, j faces
void itest_rectil_mesh_metric()
{
	double PI = 3.14159265;
	// Interior points for i + 1/2 face
	Solution testcourse = readmesh("rectil_testmesh-10x10.vel");
	Vector dx_deta(testcourse.nx, testcourse.ny);
	Vector dy_deta(testcourse.nx, testcourse.ny);
	Vector dx_dxi(testcourse.nx, testcourse.ny);
	Vector dy_dxi(testcourse.nx, testcourse.ny);

	// Testing along uniform, rectilinear grid. Expect a deta_dx of 0, deta_dy of 1, dxi_dx of 1 and dxi_dy of 0
			
	for(int i = 1; i < testcourse.nx - 1; i++)
		for(int j = 1; j < testcourse.ny - 1; j++)
		{
			dx_deta(i, j) = testcourse.X(i, j) - testcourse.X(i, j - 1);
			dy_deta(i, j) = testcourse.Y(i, j) - testcourse.Y(i, j - 1);
			dx_dxi(i, j) = (testcourse.X(i + 1, j) + testcourse.X(i + 1, j - 1) - testcourse.X(i - 1, j) - testcourse.X(i - 1, j - 1)) / 4;
			dy_dxi(i, j) = (testcourse.Y(i + 1, j) + testcourse.Y(i + 1, j - 1) - testcourse.Y(i - 1, j) - testcourse.Y(i - 1, j - 1)) / 4;
			printf("dx_deta: %f dy_deta: %f dx_dxi %f dy_dxi %f\n", dx_deta(i,j), dy_deta(i,j), dx_dxi(i,j), dy_dxi(i,j));
		}
}


void itest_simple_mesh_metric()
{
	double PI = 3.14159265;
	// Interior points for i + 1/2 face
	Solution testcourse = readmesh("simple_testmesh-10x10.vel");
	Vector dx_deta(testcourse.nx, testcourse.ny);
	Vector dy_deta(testcourse.nx, testcourse.ny);
	Vector d_xiix(testcourse.nx, testcourse.ny);
	Vector dy_dxi(testcourse.nx, testcourse.ny);

	// Testing along uniform, rectilinear grid. Expect a deta_dx of 1/2, deta_dy of 2, dxi_dx of 1/2 and dxi_dy of -2
			
	for(int i = 1; i < testcourse.nx - 1; i++)
		for(int j = 1; j < testcourse.ny - 1; j++)
		{
			dx_deta(i, j) = testcourse.X(i, j) - testcourse.X(i, j - 1);
			dy_deta(i, j) = testcourse.Y(i, j) - testcourse.Y(i, j - 1);
			d_xiix(i, j) = (testcourse.X(i + 1, j) + testcourse.X(i + 1, j - 1) - testcourse.X(i - 1, j) - testcourse.X(i - 1, j - 1)) / 4;
			dy_dxi(i, j) = (testcourse.Y(i + 1, j) + testcourse.Y(i + 1, j - 1) - testcourse.Y(i - 1, j) - testcourse.Y(i - 1, j - 1)) / 4;
			printf("dx_deta: %f dy_deta: %f d_xiix %f dy_dxi %f\n", dx_deta(i,j), dy_deta(i,j), d_xiix(i,j), dy_dxi(i,j));
		}
}

void itest_quad_mesh_metric(std::string(meshname))
{
	std::cout << meshname << std::endl;
	double PI = 3.14159265;
	// Interior points for i + 1/2 face
	Solution testcourse = readmesh(meshname);
	Vector dx_deta(testcourse.nx, testcourse.ny);
	Vector dy_deta(testcourse.nx, testcourse.ny);
	Vector dx_dxi(testcourse.nx, testcourse.ny);
	Vector dy_dxi(testcourse.nx, testcourse.ny);
	Vector dx_dxi_exact(testcourse.nx, testcourse.ny);
	Vector dy_dxi_exact(testcourse.nx, testcourse.ny);
	Vector dx_deta_exact(testcourse.nx, testcourse.ny);
	Vector dy_deta_exact(testcourse.nx, testcourse.ny);

	double xi, eta;
	double dx = 1.0/(testcourse.nx-1);
	// Testing along uniform, rectilinear grid. Expect a deta_dx of -eta, deta_dy of xi, dxi_dx of xi and dxi_dy of eta
			
	for(int i = 1; i < testcourse.nx - 1; i++)
		for(int j = 1; j < testcourse.ny - 1; j++)
		{
			xi = 0.5 + (i+0.5)*dx;
			eta = 0.5 + j*dx;
			// Calculate exact derivatives
			dx_dxi_exact(i, j) = xi*dx; dy_dxi_exact(i,j) = eta*dx;
			dx_deta_exact(i,j) = -eta*dx; dy_deta_exact(i,j) = xi*dx;
			// Calculate derivatives with second order central differencing
			dx_deta(i, j) = testcourse.X(i, j) - testcourse.X(i, j - 1);
			dy_deta(i, j) = testcourse.Y(i, j) - testcourse.Y(i, j - 1);
			dx_dxi(i, j) = (testcourse.X(i + 1, j) + testcourse.X(i + 1, j - 1) - testcourse.X(i - 1, j) - testcourse.X(i - 1, j - 1)) / 4;
			dy_dxi(i, j) = (testcourse.Y(i + 1, j) + testcourse.Y(i + 1, j - 1) - testcourse.Y(i - 1, j) - testcourse.Y(i - 1, j - 1)) / 4;
			if(i == 8 && j == 3)
			{
			printf("X: %f Y: %f\n Xi:%f\n", testcourse.X(i, j), testcourse.Y(i, j), dx_dxi(0,0));
			printf("Approx: d_etaix: %f dy_deta: %f dx_dxi %f dy_dxi %f\n", dx_deta(i,j), dy_deta(i,j), dx_dxi(i,j), dy_dxi(i,j));	
			printf("Exact: d_etaix: %f dy_deta: %f dx_dxi %f dy_dxi %f\n", dx_deta_exact(i,j), dy_deta_exact(i,j), dx_dxi_exact(i,j), dy_dxi_exact(i,j));
			}
		}
	double dx_deta_error = (dx_deta_exact - dx_deta).L2Norm(); 
	double dx_dxi_error = (dx_dxi_exact - dx_dxi).L2Norm();
	double dy_deta_error = (dy_deta_exact - dy_deta).L2Norm();
	double dy_dxi_error = (dy_dxi_exact - dy_dxi).L2Norm();

	printf("L2Norm of dx/deta: %14.12e\nL2Norm of dx/dxi: %14.12e\nL2Norm of dy/deta: %14.12e\nL2Norm of dy/dxi: %14.12e\n", dx_deta_error, dx_dxi_error, dy_deta_error, dy_dxi_error);
}

// i, j+1/2 faces
void jtest_rectil_mesh_metric()
{
	double PI = 3.14159265;
	// Interior points for i + 1/2 face
	Solution testcourse = readmesh("rectil_testmesh-10x10.vel");
	Vector d_etaix_course(testcourse.nx, testcourse.ny);
	Vector d_etaiy(testcourse.nx, testcourse.ny);
	Vector d_xiix(testcourse.nx, testcourse.ny);
	Vector d_xiiy(testcourse.nx, testcourse.ny);

	// Testing along uniform, rectilinear grid. Expect a deta_dx of 0, deta_dy of 1, dxi_dx of 1 and dxi_dy of 0
			
	for(int i = 1; i < testcourse.nx - 1; i++)
		for(int j = 1; j < testcourse.ny - 1; j++)
		{
			d_etaix_course(i, j) = testcourse.X(i, j) - testcourse.X(i, j - 1);
			d_etaiy(i, j) = testcourse.Y(i, j) - testcourse.Y(i, j - 1);
			d_xiix(i, j) = (testcourse.X(i + 1, j) + testcourse.X(i + 1, j - 1) - testcourse.X(i - 1, j) - testcourse.X(i - 1, j - 1)) / 4;
			d_xiiy(i, j) = (testcourse.Y(i + 1, j) + testcourse.Y(i + 1, j - 1) - testcourse.Y(i - 1, j) - testcourse.Y(i - 1, j - 1)) / 4;
			printf("d_etaix: %f d_etaiy: %f d_xiix %f d_xiiy %f\n", d_etaix_course(i,j), d_etaiy(i,j), d_xiix(i,j), d_xiiy(i,j));
		}
}

int main()
{
	// Different mesh files provided
	std::string coursemesh = "mesh-8x32.vel";
	std::string medmesh = "mesh-16x64.vel";
	std::string finemesh = "mesh-32x128.vel";
	std::string falsemesh = "mesh-13x2.vel"; // False mesh name to make sure catch works
	std::string fileName = "file.vtk"; // filename to store vtk file

	Solution test = readmesh(coursemesh);
	// test_void_constructor();
	/// test_constructor();
	// test_mesh_read(coursemesh);
	// test_VTK_output(test);

	// Test Mesh Files
	std::string coursequadmesh = "Quad_testmesh-10x10.vel";
	std::string medquadmesh = "Quad_testmesh-20x20.vel";
	std::string finequadmesh = "Quad_testmesh-40x40.vel";
	// Testing mesh metric generator
	itest_rectil_mesh_metric();
	// itest_simple_mesh_metric();

	// itest_quad_mesh_metric(coursequadmesh);
	// itest_quad_mesh_metric(medquadmesh);
	// itest_quad_mesh_metric(finequadmesh);


	return 0;
} 
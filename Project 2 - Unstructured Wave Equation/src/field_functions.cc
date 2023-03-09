// Included header
#include "field_functions.h"
constexpr double pi() { return std::atan(1)*4; } // Define pi

/********************************************************************
// Unstructured Mesh Programming Assignment Field Functions File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

void initial_condition(const Mesh &mesh, std::vector<double> &temp_cent, std::array<std::vector<double>, 2> &vel)
{
	// This function will initialize both the velocities u, and v at the midpoints. 
	// It will also initialize the temperature T at the cell centroid

	// Loop through edge midpoints to calculate velocity at these points
	for(int i = 0; i < mesh.iNEdge; i++)
	{
		vel[0][i] = mesh.Edge_centroid[1][i] * pi(); 
		vel[1][i] = -mesh.Edge_centroid[0][i] * pi();
	}

	// for(int i = 0; i < mesh.iNEdge; i++)
	// {
	// 	printf("Edge %i has velocity U: %f V: %f\n", i, vel[0][i], vel[1][i]);
	// }

	// Loop through cell centroids to calculate velocity at these points
	for(int i = 0; i < mesh.iNCell; i++)
	{
		temp_cent[i] = exp(-5 * (mesh.Cell_centroid[0][i] * mesh.Cell_centroid[0][i] + (mesh.Cell_centroid[1][i] - 1) * (mesh.Cell_centroid[1][i] - 1)));
	}

	// for(int i = 0; i < mesh.iNCell; i++)
	// {
	// 	printf("Cell %i has initial temp %14.12e\n", i, temp_cent[i]);
	// }
}

void save_VTK(std::string fileName, const Mesh mesh, const std::vector<double> temp_cent)
{
	const char * fileChar = fileName.c_str();
	FILE *vtkFile;
	vtkFile = fopen(fileChar, "w");

	// Print data to VTK file
	fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
	fprintf(vtkFile, "Mech 588 VTK File\n");
	fprintf(vtkFile, "ASCII\n");
	fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(vtkFile, "POINTS %i DOUBLE\n", mesh.iNVert);

	// Print vertex coordinates
	for(int i = 0; i < mesh.iNVert; i++)
	{
		fprintf(vtkFile, "%f %f %f \n", mesh.Vert[0][i], mesh.Vert[1][i], 0.0);
	}

	// Print cell connectivity
	fprintf(vtkFile, "CELLS %i %i\n", mesh.iNCell, 4*mesh.iNCell);
	for(int i = 0; i < mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%i %i %i %i\n", 3, mesh.Cell_vert[0][i], mesh.Cell_vert[1][i], mesh.Cell_vert[2][i]);
	}

	// Print cell type for each cell (triangle = 5)
	fprintf(vtkFile, "CELL_TYPES %i\n", mesh.iNCell);
	for(int i = 0; i < mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%i\n", 5);
	}

	// Create lookup table for temperature
	fprintf(vtkFile, "CELL_DATA %i\n", mesh.iNCell);
	fprintf(vtkFile, "SCALARS temperature DOUBLE 1\n");
	fprintf(vtkFile, "LOOKUP_TABLE temperature\n");
	for(int i = 0; i < mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%14.12e\n", temp_cent[i]);
	}

	fclose(vtkFile);
}

// void storeVTKSolution(Solution &s, std::string fileName)
// {
// 	const char * fileChar = fileName.c_str();
// 	FILE *vtkFile;
// 	vtkFile = fopen(fileChar, "w");
// 	unsigned long i, j;
// 	fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
// 	fprintf(vtkFile,"TITLE = \"Quad data\"\n");
// 	fprintf(vtkFile,"ASCII\n");
// 	fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
// 	fprintf(vtkFile,"DIMENSIONS %lu %lu %d\n", s.ny, s.nx, 1);
// 	fprintf(vtkFile,"POINTS %lu FLOAT\n",s.N);
// 	for(int i = 0; i < s.nx; i++)
// 		for(int j = 0; j < s.ny; j++)
// 		{
// 			fprintf(vtkFile, "%f %f %f\n", s.X(i, j), s.Y(i, j), 0.0);
// 		}
// 		// Change below to CELL_DATA
// 	fprintf(vtkFile,"CELL_DATA %lu\n", (s.nx - 1) * (s.ny - 1));
// 	fprintf(vtkFile,"SCALARS temperature FLOAT 1\n");
// 	fprintf(vtkFile,"LOOKUP_TABLE default\n");
// 	for(int i = 0; i < s.nx - 1; i++)
// 		for(int j = 0; j < s.ny - 1; j++)
// 		{
// 			double PI = 4.0*atan(1.0);
// 			s.T(i+1,j+1) = sin(PI*s.X(i+1,j+1))*sin(PI*s.Y(i+1,j+1));
// 			s.T(i,j+1) = sin(PI*s.X(i,j+1))*sin(PI*s.Y(i,j+1));
// 			s.T(i,j) = sin(PI*s.X(i,j))*sin(PI*s.Y(i,j));
// 			s.T(i+1,j) = sin(PI*s.X(i+1,j))*sin(PI*s.Y(i+1,j));

// 			double Tavg = 0.25*(s.T(i+1,j+1) + s.T(i,j) + s.T(i,j+1)+ s.T(i+1,j));
// 			fprintf(vtkFile, "%lf\n", Tavg);
// 		}
// 	// fprintf(vtkFile,"VECTORS velocity FLOAT\n");
// 	// for(int i = 0; i < s.nx; i++)
// 	// 	for(int j = 0; j < s.ny; j++)
// 	// 	{
// 	// 		fprintf(vtkFile, "%lf %lf %lf\n", s.u(i, j), s.v(i, j), 0.0);
// 	// 	}
// 	fclose(vtkFile);
// }

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

	fprintf(vtkFile, "VECTORS Velocity DOUBLE\n");
	for(int i = 0; i< mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%14.12e %14.12e %14.12e\n", mesh.Cell_centroid[1][i]*pi(), -mesh.Cell_centroid[0][i]*pi(), 0.0);
	}

	fclose(vtkFile);
}

void calc_grad(const Mesh mesh, std::vector<double> temp_cent, std::array<std::vector<double>, 2> Cell_Grad)
{
	/*
	This function will calculate the gradient in each cell. We will use the method derived in the analytical section of this project. 
	The overall order of the function is as follows:
	-> Loop through each cell
	-> For each cell, count how many edges of the cell are boundary edges. 
	-> Once we know how many boundary edges of the cell are boundary edges, we can determine which cells to use for the flux calculation. 
	For each cell that has a boundary edge, we will use another edge from a neighboring cell. For example if cell i has neighboring cells a b and c,
	we will use the a b and c to calculate the gradient in i. However, if cell i is on the boundary, and has neighbors a and b, where b is neighbors with 
	d, then we will use a b and d to calculate the gradient in i. Finally, if cell i is on the boundary, and only has neighbor a (2 edges are boundary edges),
	where a has neighbors e and f, we will use a e and f to calculate the gradient in i.
	*/
	// Loop through cells
	for(int cell = 0; cell < mesh.iNCell; cell++)
	{
		// Check how many cell edges sit on the boundary.
		int bound_edge = 0;
		for(int edge = 0; edge < 3; edge++)
		{
			if (mesh.Edge[0][mesh.Cell_edge[edge][cell]] == -1 || mesh.Edge[1][mesh.Cell_edge[edge][cell]] == -1)
			{
				bound_edge += 1;
			} 	
		}

		// If there are 0 boundary edges, we will use the centroids of each neighbor
		if(bound_edge == 0)
		{
			int neighbor[3];
			for(int i = 0; i < 3; i++)
			{
				neighbor[i] = mesh.Cell_neighbor[i][cell];
			}
			printf("Cell %i has neighbors %i %i %i\n", cell, neighbor[0], neighbor[1], neighbor[2]);
		}

		// If there is 1 boundary edge, we will use the centroids of neighbors a and b, as well as a neighbor of a
		if(bound_edge == 1)
		{
			int neighbor[2];
			for(int i = 0; i < 3; i++)
			{
				if(mesh.Cell_neighbor[i][cell] != -1)
				{
					neighbor[i] = mesh.Cell_neighbor[i][cell];
				}
			}
			printf("Cell %i has neighbors %i %i\n", cell, neighbor[0], neighbor[1]);
		}

		// If there are 2 boundary edges we will use the centroid of neighbor a, as well as two neighbors of a
		if(bound_edge == 2)
		{
			int neighbor = 9;
			for(int i = 0; i < 3; i++)
			{
				//printf("Cell %i Neighbor %i\n", cell, i);
				if(mesh.Cell_neighbor[i][cell] != -1)
				{
					neighbor = mesh.Cell_neighbor[i][cell];
				}

				printf("Cell %i: %i %i %i, %i\n", cell, mesh.Cell_neighbor[0][cell], mesh.Cell_neighbor[1][cell], mesh.Cell_neighbor[2][cell], bound_edge);
			}
			printf("Cell %i has neighbors %i\n", cell, neighbor);
		}		
		printf("Cell %i has %i boundary edges\n", cell, bound_edge);
	}
}


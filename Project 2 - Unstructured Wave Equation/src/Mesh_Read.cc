// Included libraries found in header file
#include "Mesh_Read.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/
// Constructor to initialize object
Mesh::Mesh(int ncell, int nedge, int nbdry, int nvert)
           : iNCell(ncell), iNEdge(nedge), iNBdry(nbdry), iNVert(nvert)
{
	// Initialize cell edge data to -1. They will be updated later
	for(int i = 0; i < iNCell; i++)
	{
		Cell_edge[0][i] = -1;
		Cell_edge[1][i] = -1;
		Cell_edge[2][i] = -1;
	}
}

// Member functions
void Mesh::print_edges()
{
	for(int i = 0; i < iNEdge; i++)
	{
		printf("Edge: %i %i %i %i\n", Edge[0][i], Edge[1][i], Edge[2][i], Edge[3][i]);
	}

	for(int i = 0; i < iNBdry; i++)
	{
		printf("Bdry: %i %i %i\n", Bdry[0][i], Bdry[1][i], Bdry[2][i]);
	}
}

void Mesh::print_vert_coord()
{
	for(int i = 0; i < iNVert; i++)
	{
		printf("Vertex: %i, X: %f, Y: %f\n", i, Vert[0][i], Vert[1][i]);
	}	
}

void Mesh::calc_cell_edge()
{
	// int Bdrycntr = 0;
	// int Edgecntr = 0;
	// for(int cell = 0; cell < iNCell; i++)
	// {
	// 	// Loop through edges to check if edge belongs to cell
	// 	for(int edge = 0; edge < iNEdge; edge++)
	// 	{
	// 		if(Edge[0][edge] == cell || Edge[1][edge])
	// 		{
	// 		}
	// 	}

	// 	for(int bdry = 0; bdry < iNBdry; bdry++)
	// 	{

	// 	}
	// }
}

Mesh read_mesh(std::string meshname)
{
	std::ifstream mesh(meshname); // Read mesh
	
	// Check if mesh is opened, return error otherwise
	if(mesh.is_open())
	{
		// Read number of cells, edges, boundary edges, and vertices
		int iNCell, iNEdge, iNBdry, iNVert;
		mesh >> iNCell >> iNEdge >> iNBdry >> iNVert;
		printf("Number of Cells - %i\nNumber of Edges - %i\nNumber of Boundary Edges - %i\nNumber of Vertices - %i\n", iNCell, iNEdge, iNBdry, iNVert);
		
		// Store X and Y coordinates for each vertex (X is Vert[0], Y is Vert[1])
		std::array<std::vector<double>, 2> Vert = {std::vector<double>(iNVert), std::vector<double>(iNVert)};
		for(int i = 0; i < iNVert; i++)
		{
			mesh >> Vert[0][i] >> Vert[1][i];
		}

		/*
		Store boundary edges and internal edges seperately
		For Bdry - Bdry[0]: Connected cell, Bdry[1]: Edge origin vertex, Bdry[2]: Edge destination vertex
		For Edge - Edge[0]: Left cell Edge[1]: Right cell, Edge[2]: Edge origin vertex, Edge[3]: Edge destination vertex
		*/
		std::array<std::vector<int>, 3> Bdry = {std::vector<int>(iNBdry), std::vector<int>(iNBdry), std::vector<int>(iNBdry)};
		std::array<std::vector<int>, 4> Edge = {std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge)};
		int total_edge = iNBdry + iNEdge; // Total number of edges found in the mesh
		int left, right, origin, destination; // Temporarily store values to allow for selecting which array they belong to
		// Counters to use for indices
		int Bdrycntr = 0;
		int Edgecntr = 0;
		for(int i = 0; i < total_edge; i++)
		{
			mesh >> left >> right >> origin >> destination;
			if(right == -1)
			{
				Bdry[0][Bdrycntr] = left;
				Bdry[1][Bdrycntr] = origin;
				Bdry[2][Bdrycntr] = destination;
				Bdrycntr += 1;
			}
			else
			{
				Edge[0][Edgecntr] = left;
				Edge[1][Edgecntr] = right;
				Edge[2][Edgecntr] = origin;
				Edge[3][Edgecntr] = destination;
				Edgecntr += 1;
			}
		}

		Mesh temp(iNCell, iNEdge, iNBdry, iNVert);
		// Populate array of vectors for vertex coordinate data
		for(int i = 0; i < iNVert; i++)
		{
			temp.Vert[0][i] = Vert[0][i];
			temp.Vert[1][i] = Vert[1][i];
		}

		// Populate array of vectors for boundary index data
		for(int i = 0; i < iNBdry; i++)
		{
			temp.Bdry[0][i] = Bdry[0][i];
			temp.Bdry[1][i] = Bdry[1][i];
			temp.Bdry[2][i] = Bdry[2][i];
		}

		// Populate array of vectors for edge index data
		for(int i = 0; i < iNEdge; i++)
		{
			temp.Edge[0][i] = Edge[0][i];
			temp.Edge[1][i] = Edge[1][i];
			temp.Edge[2][i] = Edge[2][i];
			temp.Edge[3][i] = Edge[3][i];
		}
		return temp;
	}

	else
	{
		printf("Mesh file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

int main()
{
// Read mesh 
std::string verycoarse = "Face-Cell/mech511-square-verycoarse.mesh";
std::string coarse = "Face-Cell/mech511-square-coarse.mesh";
std::string medium = "Face-Cell/mech511-square-medium.mesh";
std::string fine = "Face-Cell/mech511-square-fine.mesh";
std::string veryfine = "Face-Cell/mech511-square-veryfine.mesh";
std::string analytical = "Face-Cell/analytical.mesh";

Mesh mesh = read_mesh(analytical);
mesh.print_edges();

return 0;
}
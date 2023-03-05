// Included libraries
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>

/********************************************************************
// Unstructured Mesh Programming Assignment Function File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

#ifndef MESH_READ_H

#define MESH_READ_H

class Mesh
{
private:

public:
	// Store number of cells, internal edges, boundary edges, vertices
	int iNCell, iNEdge, iNBdry, iNVert;
	// Array of vectors to store vertex coordinate data
	std::array<std::vector<double>, 2> Vert = {std::vector<double>(iNVert), std::vector<double>(iNVert)}; // Store X and Y coordinates for each vertex (X is Vert[0], Y is Vert[1])
	/*
	Store connected cells, and vertex indices for edges and boundary edges seperately
	For Bdry - Bdry[0]: Connected cell, Bdry[1]: Edge origin vertex, Bdry[2]: Edge destination vertex
	For Edge - Edge[0]: Left cell Edge[1]: Right cell, Edge[2]: Edge origin vertex, Edge[3]: Edge destination vertex
	*/
	//int Edge[4][iNEdge];
	//int Bdry[3][iNBdry];
	std::array<std::vector<int>, 3> Bdry = {std::vector<int>(iNBdry), std::vector<int>(iNBdry), std::vector<int>(iNBdry)};
	std::array<std::vector<int>, 4> Edge = {std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge)};
	
	//int cell_neighbor[3][iNCell]; // Store indices of neighbors for each cell

	// Constructors
	Mesh(int ncell, int nedge, int nbdry, int nvert);

	// Member functions
	void print_edges();
	void print_vert_coord();
};

Mesh read_mesh(std::string meshname);

#endif
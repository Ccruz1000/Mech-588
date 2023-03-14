// Included libraries
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <time.h>

/********************************************************************
// Unstructured Mesh Programming Assignment Mesh Reading File Header
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
	std::array<std::vector<double>, 2> Vert = {std::vector<double>(iNVert), std::vector<double>(iNVert)}; // Store X and Y coordinates for each vertex (Xi is Vert[0][i[, Yi is Vert[1][i])
	
	/*
	Store connected cells, and vertex indices for edges and boundary edges seperately
	For Bdry - Bdry[0]: Connected cell, Bdry[1]: Edge origin vertex, Bdry[2]: Edge destination vertex
	For Edge - Edge[0][i]: Left cell, Edge[1][i]: Right cell, Edge[2][i]: Edge origin vertex, Edge[3][i]: Edge destination vertex
	*/
	std::array<std::vector<int>, 3> Bdry = {std::vector<int>(iNBdry), std::vector<int>(iNBdry), std::vector<int>(iNBdry)};
	std::array<std::vector<int>, 4> Edge = {std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge)};
	
	// Store indices of edges, vertices and neighbors for each cell
	// For each of these Cell_edge[0][i] represents the 1st "thing" of the ith cell
	std::array<std::vector<int>, 3> Cell_edge = {std::vector<int>(iNCell), std::vector<int>(iNCell), std::vector<int>(iNCell)};
	std::array<std::vector<int>, 3> Cell_neighbor = {std::vector<int>(iNCell), std::vector<int>(iNCell), std::vector<int>(iNCell)};
	std::array<std::vector<int>, 3> Cell_vert = {std::vector<int>(iNCell), std::vector<int>(iNCell), std::vector<int>(iNCell)};
	
	// Store length of each edge
	std::vector<double> Edge_length = std::vector<double>(iNEdge); 

	// Store area of each cell
	std::vector<double> Cell_area = std::vector<double>(iNCell);

	// Store x and y components for normal vector of each edge. [0][i] is x component, [1][i] is y component for edge i
	std::array<std::vector<double>, 2> Edge_norm = {std::vector<double>(iNEdge), std::vector<double>(iNEdge)};

	// Store coordinates of edge midpoint and cell centroid
	std::array<std::vector<double>, 2> Cell_centroid = {std::vector<double>(iNCell), std::vector<double>(iNCell)};
	std::array<std::vector<double>, 2> Edge_centroid = {std::vector<double>(iNEdge), std::vector<double>(iNEdge)}; 
	
	// Store upwind and downwind cell for each edge. Edge_upwind[0][i] is the upwind cell index, while Edge_upwind[1][i] is 
	// the downwind cell index. 
	std::array<std::vector<int>, 2> Edge_upwind = {std::vector<int>(iNEdge), std::vector<int>(iNEdge)};
	std::vector<double> dot_product = std::vector<double>(iNEdge); // Stores the dot product of vel*norm at each edge
	// Constructors
	Mesh(int ncell, int nedge, int nbdry, int nvert);

	// Member functions
	void print_edges();
	void print_vert_coord();
	void calc_cell_edge();
	void calc_cell_neighbor();
	void calc_cell_vert();
	void calc_edge_length();
	void calc_cell_centroid();
	void calc_edge_centroid();
	void calc_edge_norm();
	void calc_all_param(); // Calculates all of the above parameters in one go
};

Mesh read_mesh(std::string meshname);

#endif
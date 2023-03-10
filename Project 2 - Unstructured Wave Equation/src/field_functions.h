// Included libraries
#include "Mesh_Read.h"

/********************************************************************
// Unstructured Mesh Programming Assignment Field Functions File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

#ifndef FIELD_FUNCTIONS_H

#define FIELD_FUNCTIONS_H

void initial_condition(const Mesh &mesh, std::vector<double> &temp_cent, std::array<std::vector<double>, 2> &vel);
void save_VTK(std::string fileName, const Mesh mesh, const std::vector<double> temp_cent);
void calc_grad(const Mesh mesh, std::vector<double> temp_cent, std::array<std::vector<double>, 2> &Cell_Grad);

#endif
// Included libraries
#include <iostream>
#include <string> 
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <time.h>

/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

#ifndef ROE_H

#define ROE_H

void read_input(std::string input_file, double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond);
void assign_input(int argc, char* argv[], double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond);
void initial_conditions(std::array<std::vector<double>, 4> &Prim, int nx, double *init_cond, double gamma);
void save_results(std::string fileName, std::array<std::vector<double>, 4> &Prim, double dx, int nx);
void matrix_mult3x3(double A[3][3], double B[3][3], double C[3][3]);
void matrix_vector(double A[3][3], double B[3], double C[3]);
void Roe_flux(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, double dx, double gamma, std::vector<double> c, int nx);


#endif
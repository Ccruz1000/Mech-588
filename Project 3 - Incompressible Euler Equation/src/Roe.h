// Included libraries
#include <iostream>
#include <string> 
#include <fstream>
#include <vector>
#include <array>

/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

#ifndef ROE_H

#define ROE_H

void read_input(std::string input_file, double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond);
void assign_input(int argc, char* argv[], double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond);
void initial_conditions(std::vector<double> &rho, std::vector<double> &u, std::vector<double> &P, int nx, double *init_cond);
void save_results(std::string fileName, std::vector<double> rho, std::vector<double> u, std::vector<double> P, double dx, int nx);

#endif
// Included libraries
#include <iostream>
#include <string> 
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <time.h>
#include <algorithm>

/**************************************************************************************
// Compressible Flow Programming Assignment General Function File Header
// This file will contain all general purpose functions and classes. This includes 
// things such as saving files, matrix multiplication, and reading input data.
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

#ifndef FUNCTIONS_H

#define FUNCTIONS_H
// Included header files
#include "Roe.h"
#include "Steger.h"

void read_input(std::string input_file, double &CFL, int &nx, double &xbegin, double &xend, double &gamma, double *init_cond, double & t_final, int &numsave, int &timescheme);
void assign_input(int argc, char* argv[], double &CFL, int &nx, double &xbegin, double &xend, double &gamma, double *init_cond, double &t_final, int &numsave, int &timescheme);
void initial_conditions(std::array<std::vector<double>, 4> &Prim, int nx, double *init_cond, double gamma);
void run_sim(int argc, char* argv[], int scheme);
void save_results(std::string fileName, std::array<std::vector<double>, 4> &Prim, double dx, int nx);
void matrix_mult3x3(double A[3][3], double B[3][3], double C[3][3]);
void matrix_vector(double A[3][3], double B[3], double C[3]);
void calc_Prim(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, double gamma, double nx);
void calc_bound(std::array<std::vector<double>, 3> &F, const std::array<std::vector<double>, 4> Prim, double nx);
void assemble_matrices(double rho, double u, double c, double gamma, double du_dv[3][3], double XR[3][3], double XL[3][3], double dv_du[3][3]);
void assemble_lambda(double lam1, double lam2, double lam3, double lambda[3][3]);
void calc_Jacobian(double du_dv[3][3], double XR[3][3], double lambda[3][3], double XL[3][3], double dv_du[3][3], double A[3][3]);
void calc_timestep(std::array<std::vector<double>, 4> Prim, double &dt, double gamma, int nx, double CFL, double dx);
void forward_euler(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 3> &F, double dt, double dx, int nx);
void RK2(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 3> F, double dt, double dx, int nx, double gamma);
void midpoint(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 3> F, double dt, double dx, int nx, double gamma);


#endif
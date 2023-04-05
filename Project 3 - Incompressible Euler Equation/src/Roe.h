// Included libraries and headers
#include "functions.h"

/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

#ifndef ROE_H

#define ROE_H

void Roe_flux(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, int nx);

#endif
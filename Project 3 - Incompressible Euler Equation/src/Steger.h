/**************************************************************************************
// Compressible Flow Programming Assignment Steger Warming Flux Vector Splitting File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

#ifndef STEGER_H

#define STEGER_H
// Included libraries and headers
#include "functions.h"

void Steger_warming(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, int nx);

#endif
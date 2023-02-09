#include "Fcurv.h"
#include "Flux.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Flux Integral File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/



int main()
{
	// double PI = 3.14159265;
	// std::string meshname = ("testmesh-40x40.vel");
	// Solution test = readmesh(meshname);

	Vector d_etaix(test.nx, test.ny); // Vector to calculate d_eta w.r.t x on the i + 1/2 face
	Vector d_etaiy(test.nx, test.ny); // Vector to calculate d_eta w.r.t y on the i + 1/2 face
	Vector d_etajx(test.nx, test.ny); // Vector to calculate d_eta w.r.t x on the j + 1/2 face
	Vector d_etajy(test.nx, test.ny); // Vector to calculate d_eta w.r.t y on the j + 1/2 face

	Vector d_xiix(test.nx, test.ny); // Vector to calculate d_xi w.r.t x on the i + 1/2 face
	Vector d_xiiy(test.nx, test.ny); // Vector to calculate d_xi w.r.t y on the i + 1/2 face
	Vector d_xijx(test.nx, test.ny); // Vector to calculate d_xi w.r.t x on the j + 1/2 face
	Vector d_xijy(test.nx, test.ny); // Vector to calculate d_xi w.r.t j on the j + 1/2 face

	// Interior points for i + 1/2 face
	for(int i = 1; i < test.nx - 1; i++)
		for(int j = 1; j < test.ny - 1; j++)
		{
			d_etaix(i, j) = test.X(i, j) - test.X(i, j - 1);
			d_etaiy(i, j) = test.Y(i, j) - test.Y(i, j - 1);
			d_xiix(i, j) = (test.X(i + 1, j) + test.X(i + 1, j - 1) - test.X(i - 1, j) - test.X(i - 1, j - 1)) / 4;
			d_xiiy(i, j) = (test.Y(i + 1, j) + test.Y(i + 1, j - 1) - test.Y(i - 1, j) - test.Y(i - 1, j - 1)) / 4;
		}

	return 0;
}
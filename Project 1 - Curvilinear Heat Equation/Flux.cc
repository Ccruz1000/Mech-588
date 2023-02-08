#include "Fcurv.h"
#include "Flux.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Flux Integral File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/



int main()
{
	double PI = 3.14159265;
	std::string meshname = ("testmesh-20x20.vel");
	Solution test = readmesh(meshname);

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

	Vector deta_dx_exact(test.nx, test.ny);
	Vector dxi_dx_exact(test.nx, test.ny);
	Vector deta_dy_exact(test.nx, test.ny);
	Vector dxi_dy_exact(test.nx, test.ny);

	// Calculate exact values for derivatives
	for(int i = 1; i < test.nx - 1; i++)
		for(int j = 1; j < test.ny - 1; j++)
		{
			deta_dx_exact(i, j) = PI / ((test.nx - 1) * ((test.nx - 1)) * cos(test.X(i, j) * PI / (test.nx - 1)) * sin(test.Y(i, j)*PI/(test.ny - 1)));
		}
	double error = d_etaix.L2Norm() - deta_dx_exact.L2Norm();
	std::cout << "d_etaix " << d_etaix.L2Norm() << std::endl;
	std::cout << "deta_dx_exact " << deta_dx_exact.L2Norm() << std::endl;
	std::cout << "Error " << fabs(error) << std::endl;
	return 0;
}
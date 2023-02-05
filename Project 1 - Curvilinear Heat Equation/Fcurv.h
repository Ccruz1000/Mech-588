#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#ifndef FCURV_H

#define FCURV_H

// Taken from base code provided in Mech 587
class Operations{
public:
	void alloc1D(double **V, const unsigned long &N, bool isInit = true, double def = 0.0) const;
	void deAlloc1D(double **V, const unsigned long &N) const;
	void copyArray(double *S, double *T, const unsigned long &N) const;
	void addVec(double *R, const double *V1, const double *V2, const unsigned long &N) const;
	void subtractVec(double *R, const double *V1, const double *V2, const unsigned long &N) const;
	void scaleVec(double *R, const double &a, const double *V1, const unsigned long &N) const;
	void dotVec(double &R, const double *V1, const double *V2, const unsigned long &N) const;
	void findAbsMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const;
	void findMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const;
	void findMin(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const;
};

class Vector{
	unsigned long N, Nx, Ny;
	double *b;
	bool isInit;
	Operations O;

	inline bool isValid(unsigned long idxNode) const {return idxNode < N ? true:false; }

public:
	Vector();
	Vector(unsigned long Nx, unsigned long Ny, bool isInit = true, double initVal = 0.0);
	Vector(const Vector &V1);
	~Vector();
	inline size_t size() const{return N;}
	inline size_t sizeNx() const{return Nx;}
	inline size_t sizeNy() const{return Ny;}
	void setSize(unsigned long Nx, unsigned long Ny);
	double& operator()(unsigned long i, unsigned long j);
	double operator()(unsigned long i, unsigned long j) const;

	double& operator()(unsigned long i);
	double operator()(unsigned long i) const;

	Vector operator -();
	Vector operator + (const Vector &V1) const;
	Vector operator - (const Vector &V1) const;
	Vector operator * (const double &S) const;
	Vector operator + (const Vector &V1);
	Vector operator - (const Vector &V1);
	Vector operator * (const double &S);
	Vector operator = (const Vector &V1);
	friend Vector operator * (double const &S, Vector const &V);

	const double L2Norm();
	const double LinfNorm(unsigned long &i, unsigned long &j);
	const double GetMax();
	const double GetMin();

	void storeV(char filename[50]);
};


class Solution
{
private:

public:
	// Member Variables
	Vector X, Y, u, v, T; // Variables to store solution vectors
	unsigned long nx, ny; // Value to store number of X and Y points
	unsigned long N; // Store total number of points

	// Constructors
	// Null constructor, initializes vectors as empty
	Solution(); 
	// Void constructor, initializes with all 0's
	Solution(unsigned long iNx, unsigned long iNy); 
	// Constructor to initialize solution object with pre defined vectors
	Solution(unsigned long iNx, unsigned long iNy, Vector iX, Vector iY, Vector iu, Vector iv, Vector iT);
	
	// Member Functions
	void printmesh();
};

Solution readmesh(std::string meshname);
void storeVTKSolution(Solution &s, std::string fileName);

#endif
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>


/********************************************************************
// Curvillinear Mesh Programming Assignment Function File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

#ifndef FCURV_H

#define FCURV_H

// Taken from base code provided in Mech 587
/********************************************************************
// Operations from Mech 587 - Used in operator overloading
********************************************************************/
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

/********************************************************************
// Vector class from Mech 587 - Wrapper for 1D array to 2D array
********************************************************************/
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

/********************************************************************
// Matrix class from Mech 587 
********************************************************************/
class Matrix
{
	unsigned long N, Nx, Ny;
	double *A[5];
	bool isInit;
	Operations O;

	bool isValid(unsigned short pos, unsigned long idxNode) const;
public:
	Matrix();
	Matrix(unsigned long iNx, unsigned long iNy);
	~Matrix();
	inline size_t size() const {return N;}
	inline size_t sizeNx() const {return Nx;}
	inline size_t sizeNy() const {return Ny;}
	void setSize(unsigned long iNx, unsigned long iNy);
	double& operator()(unsigned long i, unsigned long j, unsigned short pos);
	double operator()(unsigned long i, unsigned long j, unsigned short pos) const;

	double& operator()(unsigned long i, unsigned short pos);
	double operator()(unsigned long i, unsigned short pos) const;

	void storeA(char filename[50]);
};
// Gauss-Seidel Solver from mech 587
void solveGS(Vector &u, const Matrix &A, const Vector &b);

/********************************************************************
// Solution class constructors and member functions
// Used to store solution, and perform any operations in solving 
********************************************************************/
class Solution
{
private:

public:
	// Member Variables
	Vector X, Y, u, v, T; // Variables to store solution vectors
	unsigned long nx, ny; // Value to store number of X and Y points
	unsigned long N; // Store total number of points
	Vector idx_deta, idy_deta, idx_dxi, idy_dxi; // Vectors for mesh metrics on i + 1/2, j faces
	Vector jdx_deta, jdy_deta, jdx_dxi, jdy_dxi; // Vectors for mesh metrics on i, j + 1/2 faces
	Vector Ji; // J component on i + 1/2, j faces
	Vector Jj; // J component on i, j + 1/2 faces
	Vector Si; // S component on i + 1/2, j faces
	Vector Sj; // S component on i, j + 1/2 faces
	Vector Di; // D component onn i + 1/2, j faces
	Vector Dj; // D component on i, j + 1/2 faces

	// Constructors
	// Null constructor, initializes vectors as empty
	Solution(); 
	// Void constructor, initializes with all 0's
	Solution(unsigned long iNx, unsigned long iNy); 
	// Constructor to initialize solution object with pre defined vectors
	Solution(unsigned long iNx, unsigned long iNy, Vector iX, Vector iY, Vector iu, Vector iv, Vector iT);
	
	// Member Functions
	void printmesh();
	void iprintmetrics(); // Print metrics for i + 1/2 faces
	void jprintmetrics(); // Print metrics for j + 1/2 faces
	void calc_mesh_metric(); // Calculate mesh metrics
};

/********************************************************************
// Useful functions - Various functions to be used in the main file
********************************************************************/
Solution readmesh(std::string meshname);  // Read mesh and create a solution object
void storeVTKSolution(Solution &s, std::string fileName); // Save VTK file from solution object



#endif
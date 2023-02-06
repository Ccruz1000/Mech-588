#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>

#ifndef BASE_H
#define BASE_H
#define INITVEC true
#define TOL 1e-8
#define PI 4.0*atan(1.0)

/**
 *	The Operations class contains the underlying operations such as memory allocation, and deallocation. 
 *  It also contains vector operations such as copying, addition, scalar multiplication, dot products, and Norms of vectors.
 */


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



class Grid;
class Vector;
class Matrix;

void solveGS(Vector &u, const Matrix &A, const Vector &b);
void solve(Vector&u, const Matrix &A, const Vector&b);
void storeVTKStructured(const Vector &u, const Grid &G, const char *fileName);
void storeVTKSolution(const Vector &u, const Vector &v, const Vector &p, const Grid &G, const char *fileName);


// Wrap all these functions in a utils class, and call them through the Utils object



/**
 * The Grid class initializes the 2-Dimensional Cartesian Grid, necessary for handling the 
 * problem. 
 * 
 * It contains the following data members:
 * X, Y 	-> Co-ordinates of Grid in X and Y directions
 * nx, ny 	-> No. of grid-points in X and Y directions.
 * dX, dY	-> Grid spacing in X and Y directions
 * xlim[2], ylim[2] -> Domain Limits in X and Y directions
 * N 		-> Total no. of grid-points in the Grid, which will be the same as the size of the Solution vector.
 * 
 * 
 * The methods for the Grid class are:
 * Grid()			-> Default constructor
 * Grid(Nx, Ny)  	-> Constructor which takes in no. of Grid Points as input, and initializes a unit square as default domain.
 * Grid(Nx, Ny, xlim, ylim) -> Constructor which takes in no. of Grid Points and domain limits as input 
 * x(i)			-> Returns the x-co-ordinate of ith point along x-direction
 * y(j)			-> Returns the y-co-ordinate of jth point along y-direction
 *  
 */

class Grid
{
	private:
	double *X, *Y;
	unsigned long nx, ny;
	double dX, dY;
	unsigned long N;
	double xlim[2], ylim[2];
	Operations O;
	
	public:
	Grid();
	Grid(unsigned long Nx, unsigned long Ny);
	Grid(unsigned long Nx, unsigned long Ny, double xlim[2], double ylim[2]);
	~Grid();
	void setGridSize(const unsigned long Nx, const unsigned long Ny);
	inline double x(const unsigned long i) const {return X[i];}
	inline double y(const unsigned long j) const {return Y[j];}
	inline double dx() const {return dX;}
	inline double dy() const {return dY;}
	inline size_t Nx() const {return nx;}
	inline size_t Ny() const {return ny;}
	inline size_t size() const {return N;}
	void storeVTKStructured(const Vector &u, const char *fname) const;
};

// class Solver
// {

// };

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

#endif
// Included libraries found in header file
#include "Mesh_Read.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

/********************************************************************
// Operations from Mech 587 - Used in operator overloading
********************************************************************/
//Allocate memory for a 1D array
void Operations::alloc1D(double **V, const unsigned long &N, bool isInit, double def) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}

	*V = new double[N];

	if(isInit)
		for(unsigned long i = 0; i < N; i++)
			(*V)[i] = 0.0;
}

// Deallocate memory for a 1D array
void Operations::deAlloc1D(double **V, const unsigned long &N) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}
}

// Copy array to new space in memory
void Operations::copyArray(double *S, double *T, const unsigned long &N) const
{
	for(int i = 0; i < N; i++)
		T[i] = S[i];
}

// Addition used for vector + operator overloading
void Operations::addVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] + V2[i];
}

// Subtraction used for vector - operator overloading
void Operations::subtractVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] - V2[i];
}

// Scale vector by scalar
void Operations::scaleVec(double *R, const double &a, const double *V1, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = a*V1[i];
}

// Take dot product of two vectors
void Operations::dotVec(double &R, const double *V1, const double *V2, const unsigned long &N) const
{
	R = 0.0;
	for(unsigned long i = 0; i < N; i++)
		R = R + (V1[i]*V2[i]);
}

// Find maximum value of the absolute value of the vector
void Operations::findAbsMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(fabs(V2[i]) > R){
			idx = i;
			R = fabs(V2[i]);
		}
}

// Find maximum value in vector
void Operations::findMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] > R){
			idx = i;
			R = V2[i];
		}
}

// Find minimum value in vector
void Operations::findMin(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = 1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] < R){
			idx = i;
			R = V2[i];
		}
}

/********************************************************************
// Vector class from Mech 587 - Wrapper for 1D array to 2D array
********************************************************************/
// Constructors
Vector::Vector() : Nx(0), Ny(0), N(0), b(nullptr), isInit(false){}

Vector::Vector(unsigned long iNx, unsigned long iNy, bool initflag, double initVal): Nx(iNx), Ny(iNy), N(iNx*iNy), isInit(false), b(nullptr)
{
	O.alloc1D(&b,N,initflag,initVal);
	isInit = true;
}

Vector::Vector(const Vector &V1) : Nx(V1.Nx), Ny(V1.Ny), N(V1.N), isInit(false), b(nullptr)
{
	O.alloc1D(&b,N,false);
	isInit = true;
	O.copyArray(V1.b, b, N);
}

// Destructor
Vector::~Vector()
{
	O.deAlloc1D(&b,N);
	isInit = false;
	Nx = 0, Ny = 0; N = 0;
}

// Set size of vector and allocate memory 
void Vector::setSize(unsigned long iNx, unsigned long iNy)
{
	Nx = iNx, Ny = iNy;
	N = Nx*Ny;
	b = nullptr;
	O.alloc1D(&b,N);
	isInit = true;
}

// Operator overloads
double& Vector::operator()(unsigned long i, unsigned long j)
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	unsigned long idxNode = i*Ny+j;
	if(isInit && isValid(idxNode))
		return b[idxNode];
	else{
		printf("attempted to access invalid memory location. Please check. Exiting\n");
		exit(0);
	}
}

double Vector::operator()(unsigned long i, unsigned long j) const
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}
	unsigned long idxNode = i*Ny+j;
	if(isInit && isValid(idxNode))
		return b[idxNode];
	else
		return 0.0;
}

double& Vector::operator()(unsigned long i)
{
	if(!isInit){
		printf("Vector has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(i))
		return b[i];
	else{
		printf("attempted to access invalid memory location. Please check. Exiting\n");
		exit(0);
	}
}

double Vector::operator()(unsigned long i) const
{
	if(!isInit){
		printf("Vector has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(i))
		return b[i];
	else
		return 0.0;
}

Vector Vector::operator = (const Vector &V1)
{
	Nx = V1.Nx;
	Ny = V1.Ny;

	if(N == V1.N)
		O.copyArray(V1.b, b, V1.N);
	else{
		O.deAlloc1D(&b, N);
		N = V1.N;
		O.alloc1D(&b, N);
		O.copyArray(V1.b, b, N);
	}

	return *this;
}

Vector Vector::operator - ()
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b, -1.0, this->b, this->N);

	return R;
}

Vector Vector::operator + (Vector const &V)
{
	if(V.N != this->N){
		printf("Invalid vector addition operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny, false);
	O.addVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator - (const Vector &V)
{
	if(V.N != this->N){
		printf("Invalid vector subtraction operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny,false);
	O.subtractVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator * (const double &a)
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b,a,this->b,this->N);
	return R;
}

Vector Vector::operator + (Vector const &V) const
{
	if(V.N != this->N){
		printf("Invalid vector addition operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny, false);
	O.addVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator - (const Vector &V) const
{
	if(V.N != this->N){
		printf("Invalid vector subtraction operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny,false);
	O.subtractVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator * (const double &a) const
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b,a,this->b,this->N);
	return R;
}

Vector operator * (const double &a, const Vector &V1)
{
	Operations O;
	Vector R(V1.Nx, V1.Ny, false);
	O.scaleVec(R.b,a,V1.b,V1.N);
	return R;
}

// Computer L2Norm of vector
const double Vector::L2Norm()
{
	double L2;
	O.dotVec(L2,this->b, this->b, this->N);
	L2 = sqrt(L2/(this->N));
	return L2;
}

// Computer L infinity norm of vector
const double Vector::LinfNorm(unsigned long &ix, unsigned long &iy)
{
	double Linf = -1e20;
	unsigned long idx;
	O.findAbsMax(Linf, idx, this->b, this->N);
	ix = idx/Ny;
	iy = idx%Ny;
	return Linf;
}

// Find maximum value in vector
const double Vector::GetMax()
{
	double max = -1e20;
	unsigned long idx;
	O.findMax(max,idx,this->b, this->N);
	return max;
}

// Find minimum value in vector
const double Vector::GetMin()
{
	double min = 1e20;
	unsigned long idx;
	O.findMin(min, idx, this->b, this->N);
	return min;
}

// Store vector in file
void Vector::storeV(char filename[50])
{
	FILE *f;
	f = fopen(filename,"w");
	for(unsigned long i = 0; i < Nx; i++){
		for(unsigned long j = 0; j < Ny; j++)
			fprintf(f,"%lu %lu %10e\n",i, j, (*this)(i,j));
		fprintf(f,"\n");
	}
}

/********************************************************************
// Matrix class from Mech 587 
********************************************************************/
Matrix::Matrix(): isInit(false)
{
	N = 0, Nx = 0, Ny = 0;
	for(int i = 0; i < 5; i++)
		A[i] = nullptr;
}

Matrix::Matrix(unsigned long iNx, unsigned long iNy) : Nx(iNx), Ny(iNy), N(iNx*iNy), isInit(false)
{
	for(int i = 0; i < 5; i++){
		A[i] = nullptr;
		O.alloc1D(&A[i], N);
	}
	isInit = true;

}

Matrix::~Matrix()
{
	for(int i = 0; i < 5; i++)
		O.deAlloc1D(&A[i], N);
	Nx = 0, Ny = 0, N = 0;
}

bool Matrix::isValid(unsigned short pos, unsigned long idxNode) const
{
	if(idxNode < N && (pos >= 0 && pos <= 4))
		return true;
	return false;
}

double& Matrix::operator()(unsigned long i, unsigned long j, unsigned short pos)
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	unsigned long idxNode = i*Ny+j;
	if(isValid(pos, idxNode))
		return A[pos][idxNode];
	else
	{
		printf("Attempt to access invalid memory. Please check. Exiting \n");
		exit(0);
		//Trigger exit code
	}
}

double Matrix::operator()(unsigned long i, unsigned long j, unsigned short pos) const
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	unsigned long idxNode = i*Ny+j;
	if(isValid(pos, idxNode))
		return A[pos][idxNode];
	else
		return 0.0;
}

double& Matrix::operator()(unsigned long i, unsigned short pos)
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(pos, i))
		return A[pos][i];
	else
	{
		printf("Attempt to access invalid memory. Please check. Exiting \n");
		exit(0);
		//Trigger exit code
	}
}

double Matrix::operator()(unsigned long i, unsigned short pos) const
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(pos, i))
		return A[pos][i];
	else
		return 0.0;
}

void Matrix::setSize(unsigned long iNx, unsigned long iNy){
	Nx = iNx, Ny = iNy;
	N = Nx*Ny;
	for(int i = 0; i < 5; i++){
		A[i] = nullptr;
		O.alloc1D(&A[i], N);
	}
	isInit = true;
}

// Gauss-Seidel Solver from Mech 587 
void solveGS(Vector &u, const Matrix &A, const Vector &b)
{
	size_t Ny = b.sizeNy();
	unsigned long i;
	const Vector* const u0 = &u;

	int iter = 0;
	while(1){
		double L2Norm = 0.0;
		for(i = 0; i < A.size(); i++){
			double Res = 0.0;
			Res = Res + A(i,0) * (*u0)(i-Ny);
			Res = Res + A(i,1) * (*u0)(i-1);
			Res = Res + A(i,2) * (*u0)(i);
			Res = Res + A(i,3) * (*u0)(i+1);
			Res = Res + A(i,4) * (*u0)(i+Ny);

			Res = b(i) - Res;

			L2Norm += (Res*Res);
			u(i) = u(i) + 1.8*(1.0/A(i,2))*Res;
		}
		L2Norm = sqrt(L2Norm/b.size());
		++iter;
		// printf("%d %14.12e\n",iter,L2Norm);
		if(L2Norm < 1e-8 || iter > 50000)
		{
			printf("Gauss Seidel Iterations converged\n");
			printf("iterations = %d, Error Residuals = %14.12e\n",iter,L2Norm);
			break;
		}
	}
}

/********************************************************************
Class used to define mesh
The following data is contained within this object
- Number of cells, vertices, boundary edges, and edges
- Utilizes .mesh file type

********************************************************************/
// class Mesh
// {
// private:


// public:
// 	 Mesh size data 
// 	int iNVerts, iNEdges, iNBdryEdges, iNCells;
// }

// /********************************************************************
// Class used to define individual cells when reading meshes
// The following data is required for each cell
// - Cell centroid location
// - Cell vertices
// - Cell face center locations
// ********************************************************************/
// class Cell
// {
// private:

// public:
// double (*cell_cent)[NDIM];
// };

// /********************************************************************
// Class used to define individual cells when reading meshes
// The following data is required for each vertex
// - Cell centroid location
// - Cell vertices
// - Cell face center locations
// ********************************************************************/
int main()
{
// Read mesh 
std::string meshname = "../Face-Cell/mech511-square-verycoarse.mesh";


std::ifstream mesh(meshname);

int iNCell, iNEdge, iNBdry, iNVert;

mesh >> iNCell >> iNEdge >> iNBdry >> iNVert;

printf("%i %i %i %i \n", iNCell, iNEdge, iNBdry, iNVert);
}
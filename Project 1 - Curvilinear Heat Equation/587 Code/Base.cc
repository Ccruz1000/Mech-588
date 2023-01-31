#include "Base.h"

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

void Operations::deAlloc1D(double **V, const unsigned long &N) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}
}

void Operations::copyArray(double *S, double *T, const unsigned long &N) const
{
	for(int i = 0; i < N; i++)
		T[i] = S[i];
}

void Operations::addVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] + V2[i];
}

void Operations::subtractVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] - V2[i];
}

void Operations::scaleVec(double *R, const double &a, const double *V1, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = a*V1[i];
}

void Operations::dotVec(double &R, const double *V1, const double *V2, const unsigned long &N) const
{
	R = 0.0;
	for(unsigned long i = 0; i < N; i++)
		R = R + (V1[i]*V2[i]);
}

void Operations::findAbsMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(fabs(V2[i]) > R){
			idx = i;
			R = fabs(V2[i]);
		}
}

void Operations::findMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] > R){
			idx = i;
			R = V2[i];
		}
}

void Operations::findMin(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = 1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] < R){
			idx = i;
			R = V2[i];
		}
}

Grid::Grid() : X(nullptr), Y(nullptr), nx(0), ny(0), N(0)
{}

Grid::Grid(unsigned long iNx, unsigned long iNy) : nx(iNx), ny(iNy), X(nullptr), Y(nullptr)
{
	xlim[0] = 0.0, xlim[1] = 1.0;
	ylim[0] = 0.0, ylim[1] = 1.0;
	N = nx*ny;
	dX = 1.0/(nx-1);
	dY = 1.0/(ny-1);

	O.alloc1D(&X, nx);
	O.alloc1D(&Y, ny);

	unsigned long i;
	for(i = 0; i < nx; i++)
		X[i] = i*dX;

	for(i = 0; i < ny; i++)
		Y[i] = i*dY;
}

Grid::Grid(unsigned long iNx, unsigned long iNy, double ixlim[2], double iylim[2]) : nx(iNx), ny(iNy), X(nullptr), Y(nullptr)
{
	xlim[0] = ixlim[0], xlim[1] = ixlim[1];
	ylim[0] = iylim[0], ylim[1] = iylim[1];
	N = nx*ny;
	dX = (xlim[1]-xlim[0])/(nx-1);
	dY = (ylim[1]-ylim[0])/(ny-1);

	O.alloc1D(&X, nx);
	O.alloc1D(&Y, ny);

	unsigned long i;
	for(i = 0; i < nx; i++)
		X[i] = xlim[0] + i*dX;

	for(i = 0; i < ny; i++)
		Y[i] = ylim[0] + i*dY;
}

Grid::~Grid()
{
	O.deAlloc1D(&X, nx);
	O.deAlloc1D(&Y, ny);
	nx = 0, ny = 0;
	N = 0;
	dX = 0.0, dY = 0.0;
}

void Grid::setGridSize(const unsigned long iNx, const unsigned long iNy)
{
	if(X != nullptr && Y != nullptr)
		return;

	nx = iNx, ny = iNy;
	dX = 1.0/(nx-1); dY = 1.0/(ny-1);
	N = nx*ny;
	O.alloc1D(&X,nx);
	O.alloc1D(&Y,ny);

	unsigned long i;
	for(i = 0; i < nx; i++)
		X[i] = i*dX;

	for(i = 0; i < ny; i++)
		Y[i] = i*dY;
}

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

Vector::~Vector()
{
	O.deAlloc1D(&b,N);
	isInit = false;
	Nx = 0, Ny = 0; N = 0;
}

void Vector::setSize(unsigned long iNx, unsigned long iNy)
{
	Nx = iNx, Ny = iNy;
	N = Nx*Ny;
	b = nullptr;
	O.alloc1D(&b,N);
	isInit = true;
}

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



const double Vector::L2Norm()
{
	double L2;
	O.dotVec(L2,this->b, this->b, this->N);
	L2 = sqrt(L2/(this->N));
	return L2;
}

const double Vector::LinfNorm(unsigned long &ix, unsigned long &iy)
{
	double Linf = -1e20;
	unsigned long idx;
	O.findAbsMax(Linf, idx, this->b, this->N);
	ix = idx/Ny;
	iy = idx%Ny;
	return Linf;
}

const double Vector::GetMax()
{
	double max = -1e20;
	unsigned long idx;
	O.findMax(max,idx,this->b, this->N);
	return max;
}

const double Vector::GetMin()
{
	double min = 1e20;
	unsigned long idx;
	O.findMin(min, idx, this->b, this->N);
	return min;
}


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

void Matrix::storeA(char filename[50])
{
	FILE *f;
	f = fopen(filename,"w");
	for(unsigned long i = 0; i < Nx; i++){
		for(unsigned long j = 0; j < Ny; j++)
			fprintf(f,"%lu %lu %10e %10e %10e %10e %10e \n",i, j, (*this)(i,j,0), (*this)(i,j,1), (*this)(i,j,2), (*this)(i,j,3), (*this)(i,j,4));
		fprintf(f,"\n");
	}
}

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

void solve(Vector &u, const Matrix &A, const Vector &b)
{
	solveGS(u,A,b);
}

void storeVTKStructured(const Vector &u, const Grid &G, const char *fileName) 
{
	FILE *vtkFile;
	vtkFile = fopen(fileName,"w");
	unsigned long i,j;
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	fprintf(vtkFile,"# vtk DataFile Version 2.0\n");
	fprintf(vtkFile,"TITLE = \"Quad data\"\n");
	fprintf(vtkFile,"ASCII\n");
	fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
	fprintf(vtkFile,"DIMENSIONS %lu %lu %d\n", Nx, Ny, 1);
	fprintf(vtkFile,"POINTS %lu FLOAT\n",Nx*Ny);
	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){
			fprintf(vtkFile, "%lf %lf %lf\n",G.x(i), G.y(j), 0.0);
		}

	fprintf(vtkFile,"POINT_DATA %lu\n",Nx*Ny);
	fprintf(vtkFile,"SCALARS phi FLOAT 1\n");
	fprintf(vtkFile,"LOOKUP_TABLE default\n");
	for(i= 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
			fprintf(vtkFile,"%14.12lf\n",u(i,j));

	fclose(vtkFile);
}

void storeVTKSolution(const Vector &u, const Vector &v, const Vector &p, const Grid &G, const char *fileName) 
{;
	FILE *vtkFile;
	vtkFile = fopen(fileName,"w");
	unsigned long i,j;
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	fprintf(vtkFile,"# vtk DataFile Version 2.0\n");
	fprintf(vtkFile,"TITLE = \"Quad data\"\n");
	fprintf(vtkFile,"ASCII\n");
	fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
	fprintf(vtkFile,"DIMENSIONS %lu %lu %d\n", Nx, Ny, 1);
	fprintf(vtkFile,"POINTS %lu FLOAT\n",Nx*Ny);
	for(j = 0; j < Ny; j++)
		for(i = 0; i < Nx; i++){
			fprintf(vtkFile, "%f %f %f\n",G.x(i), G.y(j), 0.0);
		}

	fprintf(vtkFile,"POINT_DATA %lu\n",Nx*Ny);
	fprintf(vtkFile,"SCALARS pressure FLOAT 1\n");
	fprintf(vtkFile,"LOOKUP_TABLE default\n");
	for(j= 0; j < Ny; j++)
		for(i = 0; i < Nx; i++)
			fprintf(vtkFile,"%lf\n",p(i,j));

	fprintf(vtkFile,"VECTORS velocity FLOAT\n");
	for(j= 0; j < Ny; j++)
		for(i = 0; i < Nx; i++)
			fprintf(vtkFile,"%lf %lf %lf\n",u(i,j),v(i,j),0.0);

	fclose(vtkFile);
}

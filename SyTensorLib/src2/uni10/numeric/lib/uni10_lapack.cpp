#include "mkl.h"
#include "mkl_lapack.h"
#include <string.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
namespace uni10{
void myDgemm(double* A, double* B, int M, int N, int K, double* C){
	double alpha = 1, beta = 0;
	dgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
}

void vecAdd(double* X, double* Y, int64_t N){
	double a = 1.0;
	int inc = 1;
	int64_t left = N;
	int64_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		daxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
	}
}
void vecScal(double a, double* X, int64_t N){
	int inc = 1;
	int64_t left = N;
	int64_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		dscal(&chunk, &a, X + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
	}
}

/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N){
	int eleNum = M*N;
	double *random = (double*)malloc(eleNum * sizeof(double));
	randomNums(random, M * N, 0);
	assert(M <= N);
	int min = M; //min = min(M,N)
	int ldA = M, ldu = M, ldvT = min;
	double *S = (double*)malloc(min*sizeof(double));
	double *u = (double*)malloc(ldu*M*sizeof(double));
	int lwork = 16*N;
	double *work = (double*)malloc(lwork*sizeof(double));
	int info;
	//tmpT = u*S*vT
	dgesvd((char*)"N", (char*)"S", &M, &N, random, &ldA, S, u, &ldu, elem, &ldvT, work, &lwork, &info);
	//reshape from Fortran format to C format
	memcpy(random, elem, eleNum * sizeof(double));
	myTranspose(random, N, M, elem, 0);
	free(random);
	free(S);
	free(u);
	free(work);
}
void syDiag(double* Kij, int N, double* Eig, double* EigVec){
	memcpy(EigVec, Kij, N * N * sizeof(double));
	int ldA = N;
	int lwork = 4*N;
	double* work= (double*)malloc(sizeof(double)*lwork);
	int info;
	dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, &info);
	assert(info == 0);
	free(work);
}

void myDgesvd(double* Mij_ori, int M, int N, double* U, double* S, double* vT){ //not tested yet
	//Mij = U * S * VT
	double* Mij = (double*)malloc(M * N * sizeof(double));
	memcpy(Mij, Mij_ori, M * N * sizeof(double));
	int min = M < N ? M : N;	//min = min(M,N)
	int ldA = N, ldu = N, ldvT = min;
	int lwork = 12*N;
	double *work = (double*)malloc(lwork*sizeof(double));
	int info;
	dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, &info);
	assert(info == 0);
	free(work);
	free(Mij);
}
void myTranspose(double* A, int M, int N, double* AT, int status){
	for(int i = 0; i < M; i++)
		for(int j = 0; j < N; j++)
			AT[j * M + i] = A[i * N + j];
}

void myEye(double* elem, int M, int N, int status){
	int min;
	if(M < N)	min = M;
	else		min = N;
	memset(elem, 0, M * N * sizeof(double));
	for(int i = 0; i < min; i++)
		elem[i * N + i] = 1;
}

};	/* namespace uni10 */	

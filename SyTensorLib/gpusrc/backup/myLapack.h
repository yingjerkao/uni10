#pragma once
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include "mkl.h"
#include "mkl_lapack.h"
#include "myTools.h"
#include <iostream>
using namespace std;
#include <stdio.h>


void myDgemm(double* A, double* B, int M, int N, int K, double* C);
void vecAdd(double* X, double* Y, int64_t N);	// Y = X + Y
void vecScal(double a, double* X, int64_t N);	// X = a * X

/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N, int status);
void syDiag(double* Kij, int N, double* Eig, double* EigVec);
void myDgesvd(double* Mij_ori, int M, int N, double* U, double* S, double* vT);

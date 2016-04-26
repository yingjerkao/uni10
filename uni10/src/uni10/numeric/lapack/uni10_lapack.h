/****************************************************************************
*  @file uni10_lapack.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Header file for Lapack wrapping functions
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef UNI10_LAPACK_H
#define UNI10_LAPACK_H
#include <cstdint>
#include <limits.h>
#include <assert.h>
#include <sstream>
#include <stdexcept>
#include <cmath>
//#include <Accelerate/Accelerate.h>
#include <complex>
namespace uni10{
enum mmtype{
	MM_DDD = 0,
	MM_DDH = 1,
	MM_DHD = 2,
	MM_DHH = 3,
	MM_HDD = 4,
	MM_HDH = 5,
	MM_HHD = 6,
	MM_HHH = 7
};
void uni10Dgemm(int p, int q, int M, int N, int K, double* A, double* B, double* C, mmtype how);
void matrixMul(double* A, double* B, int M, int N, int K, double* C, bool ongpuA, bool ongpuB, bool ongpuC);
void vectorAdd(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu);// Y = Y + X
void vectorScal(double a, double* X, size_t N, bool ongpu);	// X = a * X
void vectorMul(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu); // Y = Y * X, element-wise multiplication;
double vectorSum(double* X, size_t N, int inc, bool ongpu);
double vectorNorm(double* X, size_t N, int inc, bool ongpu);
void vectorExp(double a, double* X, size_t N, bool ongpu);
void diagRowMul(double* mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu);
void diagColMul(double* mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu);
/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N, bool ongpu);
void eigDecompose(double* Kij, int N, std::complex<double>* Eig, std::complex<double> *EigVec, bool ongpu);
void eigSyDecompose(double* Kij, int N, double* Eig, double* EigVec, bool ongpu);
void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT, bool ongpu);
void matrixInv(double* A, int N, bool diag, bool ongpu);
void setTranspose(double* A, size_t M, size_t N, double* AT, bool ongpu, bool ongpuT);
void setTranspose(double* A, size_t M, size_t N, bool ongpu);
void setCTranspose(double* A, size_t M, size_t N, double *AT, bool ongpu, bool ongpuT);
void setCTranspose(double* A, size_t M, size_t N, bool ongpu);
void setIdentity(double* elem, size_t M, size_t N, bool ongpu);
void reshapeElem(double* elem, size_t* transOffset);
bool lanczosEV(double* A, double* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, double* eigVec, bool ongpu);
//====== real qr rq ql lq ======//
void matrixQR(double* Mij_ori, int M, int N, double* Q, double* R, bool ongpu);
void matrixRQ(double* Mij_ori, int M, int N, double* Q, double* R, bool ongpu);
void matrixQL(double* Mij_ori, int M, int N, double* Q, double* L, bool ongpu);
void matrixLQ(double* Mij_ori, int M, int N, double* Q, double* L, bool ongpu);
//==============================//
/***** Complex version *****/
void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, double *S, std::complex<double>* vT, bool ongpu);
void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, std::complex<double>* S, std::complex<double>* vT, bool ongpu);
void matrixInv(std::complex<double>* A, int N, bool diag, bool ongpu);
std::complex<double> vectorSum(std::complex<double>* X, size_t N, int inc, bool ongpu);
double vectorNorm(std::complex<double>* X, size_t N, int inc, bool ongpu);
void matrixMul(std::complex<double>* A, std::complex<double>* B, int M, int N, int K, std::complex<double>* C, bool ongpuA, bool ongpuB, bool ongpuC);
void vectorAdd(std::complex<double>* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu);// Y = Y + X
void vectorAdd(std::complex<double>* Y, std::complex<double>* X, size_t N, bool y_ongpu, bool x_ongpu);// Y = Y + X
void vectorScal(double a, std::complex<double>* X, size_t N, bool ongpu);	// X = a * X
void vectorScal(const std::complex<double>& a, std::complex<double>* X, size_t N, bool ongpu);	// X = a * X
void vectorMul(std::complex<double>* Y, std::complex<double>* X, size_t N, bool y_ongpu, bool x_ongpu); // Y = Y * X, element-wise multiplication;
void diagRowMul(std::complex<double>* mat, std::complex<double>* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu);
void diagColMul(std::complex<double>* mat, std::complex<double>* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu);
void vectorExp(double a, std::complex<double>* X, size_t N, bool ongpu);
void vectorExp(const std::complex<double>& a, std::complex<double>* X, size_t N, bool ongpu);
void orthoRandomize(std::complex<double>* elem, int M, int N, bool ongpu);
void setTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double>* AT, bool ongpu, bool ongpuT);
void setTranspose(std::complex<double>* A, size_t M, size_t N, bool ongpu);
void setCTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double>* AT, bool ongpu, bool ongpuT);
void setCTranspose(std::complex<double>* A, size_t M, size_t N, bool ongpu);
void eigDecompose(std::complex<double>* Kij, int N, std::complex<double>* Eig, std::complex<double> *EigVec, bool ongpu);
void eigSyDecompose(std::complex<double>* Kij, int N, double* Eig, std::complex<double>* EigVec, bool ongpu);
void setConjugate(std::complex<double> *A, size_t N, bool ongpu);
void setIdentity(std::complex<double>* elem, size_t M, size_t N, bool ongpu);
bool lanczosEV(std::complex<double>* A, std::complex<double>* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, std::complex<double>* eigVec, bool ongpu);
bool lanczosEVL(std::complex<double>* A, std::complex<double>* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, std::complex<double>* eigVec, bool ongpu);
void matrixQR(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R, bool ongpu);
void matrixRQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R, bool ongpu);
void matrixQL(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L, bool ongpu);
void matrixLQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L, bool ongpu);

};	/* namespace uni10 */
#endif /* UNI10_LAPACK_H */

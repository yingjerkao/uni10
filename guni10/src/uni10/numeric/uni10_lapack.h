/****************************************************************************
*  @file uni10_lapack.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao
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
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
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
double vectorSum(double* X, size_t N, int inc, bool ongpu);
double vectorNorm(double* X, size_t N, int inc, bool ongpu);
void vectorExp(double a, double* X, size_t N, bool ongpu);

/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N, bool ongpu);

void syDiag(double* Kij, int N, double* Eig, double* EigVec, bool ongpu);
void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT, bool ongpu);
void setTranspose(double* A, size_t M, size_t N, double* AT, bool ongpu);
void setIdentity(double* elem, size_t M, size_t N, bool ongpu);
void reshapeElem(double* elem, size_t* transOffset);

};	/* namespace uni10 */
#endif /* UNI10_LAPACK_H */

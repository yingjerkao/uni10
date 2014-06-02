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
namespace uni10{
void myDgemm(double* A, double* B, int M, int N, int K, double* C);
void vecAdd(double* X, double* Y, int64_t N);	// Y = X + Y
void vecScal(double a, double* X, int64_t N);	// X = a * X

/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N);

void syDiag(double* Kij, int N, double* Eig, double* EigVec);
void myDgesvd(double* Mij_ori, int M, int N, double* U, double* S, double* vT);
void myTranspose(double* A, int M, int N, double* AT, int status);
void myEye(double* elem, int M, int N, int status);

};	/* namespace uni10 */	
#endif /* UNI10_LAPACK_H */

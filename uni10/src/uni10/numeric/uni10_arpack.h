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
*  @brief Header file for Arapack wrapping functions
*  @author Chen-Yen Lai
*  @date 2016-01-22
*  @since 0.1.0
*
*****************************************************************************/
#ifndef UNI10_ARPACK_H
#define UNI10_ARPACK_H
#include <complex>

bool arpackEigh(double* A, double* psi, size_t n, size_t& max_iter,
    double& eigVal, double* eigVec, bool ongpu, double err_tol=0.0e0, int nev=1);

bool arpackEigh(std::complex<double>* A, std::complex<double>* psi, size_t n,
    size_t& max_iter, double& eigVal, std::complex<double>* eigVec, bool ongpu,
    double err_tol=0.0e0, int nev=1);

#endif /* end of include guard: UNI10_ARPACK_H */

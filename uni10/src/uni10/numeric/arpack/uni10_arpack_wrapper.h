/****************************************************************************
*  @file uni10_lapack_wrapper.h
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
*  @brief C wrapper functions for Arpack libraries
*  @author Chen-Yen Lai
*  @date 2016-01-22
*  @since 0.9.2
*
*****************************************************************************/
#ifndef UNI10_ARPACK_WRAPPER_H
#define UNI10_ARPACK_WRAPPER_H
#include <complex>

extern "C" {
void znaupd_(int *ido, char *bmat, int *n, char *which,
             int *nev, double *tol, std::complex<double> *resid, int *ncv,
             std::complex<double> *v, int *ldv, int *iparam, int *ipntr,
             std::complex<double> *workd, std::complex<double> *workl,
             int *lworkl, double *rwork, int *info);

void zneupd_(int *rvec, char *All, int *select, std::complex<double> *d,
             std::complex<double> *z, int *ldz, std::complex<double> *sigma,
             std::complex<double> *workev, char *bmat, int *n, char *which, int *nev,
             double *tol, std::complex<double> *resid, int *ncv,
             std::complex<double> *v,
             int *ldv, int *iparam, int *ipntr, std::complex<double> *workd,
             std::complex<double> *workl, int *lworkl, double *rwork, int *info);

void dsaupd_(int *ido, char *bmat, int *n, char *which,
            int *nev, double *tol, double *resid, int *ncv,
            double *v, int *ldv, int *iparam, int *ipntr,
            double *workd, double *workl, int *lworkl, int *info);

void dseupd_(int *rvec, char *All, int *select, double *d,
            double *z, int *ldz, double *sigma,
            char *bmat, int *n, char *which, int *nev,
            double *tol, double *resid, int *ncv, double *v,
            int *ldv, int *iparam, int *ipntr, double *workd,
            double *workl, int *lworkl, int *info);
}

#endif /* end of include guard: UNI10_ARPACK_WRAPPER_H */

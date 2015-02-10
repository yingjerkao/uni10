/****************************************************************************
*  @file Matrix.h
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
*  @brief Header file for Matrix class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <string.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <uni10/data-structure/Block.h>
#include <uni10/data-structure/CBlock.h>
//Type of Matrix
namespace uni10{

class Matrix: public Block {
	public:
		Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=true);
		Matrix(size_t _Rnum, size_t _Cnum, const double* _elem, bool _diag=false, bool src_ongpu=false);
		Matrix(size_t _Rnum, size_t _Cnum, const std::vector<double>& _elem, bool _diag=false, bool src_ongpu=false);
		Matrix(const Matrix& _m);
    Matrix(const Block& _b);
#ifndef UNI10_PURE_REAL
    Matrix(const CBlock& _b);
#endif
		Matrix();
		~Matrix();
		Matrix& operator=(const Matrix& _m);
		Matrix& operator=(const Block& _m);
		double& operator[](size_t idx);
		double& at(size_t i, size_t j);
		double* getHostElem();
		void setElem(const double* elem, bool _ongpu = false);
		void setElem(const std::vector<double>& elem, bool _ongpu = false);
		Matrix& resize(size_t row, size_t col);
		void load(const std::string& fname);
		void identity();
		void set_zero();
		void randomize();
		void orthoRand();
		bool toGPU();
		Matrix& transpose();
		Matrix& cTranspose();
    Matrix& conj(){return *this;};
		Matrix& operator*= (double a);
		Matrix& operator*= (const Block& Mb);
		Matrix& operator+= (const Block& Mb);
	private:
		void init(bool togpu);
		void init(const double* elem, bool _ongpu);
};
Matrix takeExp(double a, const Block& mat);
#ifndef UNI10_PURE_REAL
Matrix exp(double a, const Block& mat);
CMatrix exp(const std::complex<double>& a, const Block& mat);
#endif
Matrix exph(double a, const Block& mat);
#ifndef UNI10_PURE_REAL
Matrix exp(const Block& mat);
#endif
Matrix exph(const Block& mat);
Matrix otimes(const Block& Ma, const Block& Mb);
};	/* namespace uni10 */
#endif /* MATRIX_H */

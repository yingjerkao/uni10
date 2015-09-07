/****************************************************************************
*  @file CMatrix.h
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
#ifndef CMATRIX_H
#define CMATRIX_H
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
//Type of CMatrix
namespace uni10{
class CMatrix: public CBlock {
	public:
		CMatrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=true);
		CMatrix(size_t _Rnum, size_t _Cnum, const std::complex<double>* _elem, bool _diag=false, bool src_ongpu=false);
		CMatrix(size_t _Rnum, size_t _Cnum, const std::vector<std::complex<double> >& _elem, bool _diag=false, bool src_ongpu=false);
		CMatrix(const CMatrix& _m);
    CMatrix(const CBlock& _b);  
    CMatrix(const Block& _b);   //////
		CMatrix();
		~CMatrix();
		CMatrix& operator=(const CMatrix& _m);
		CMatrix& operator=(const CBlock& _m);
		std::complex<double>& operator[](size_t idx);
		std::complex<double>& at(size_t i, size_t j);
		std::complex<double>* getHostElem();
		void setElem(const std::complex<double>* elem, bool _ongpu = false);
		void setElem(const std::vector<std::complex<double> >& elem, bool _ongpu = false);
		CMatrix& resize(size_t row, size_t col);
		void load(const std::string& fname);
		void identity();
		void set_zero();
		void randomize();
		void orthoRand();
		bool toGPU();
		CMatrix& transpose();
		CMatrix& cTranspose();
    CMatrix& conj();
		CMatrix& operator*= (double a);
		CMatrix& operator*= (const std::complex<double>& a);
		CMatrix& operator*= (const Block& Mb);
		CMatrix& operator*= (const CBlock& Mb);
		CMatrix& operator+= (const Block& Mb);
		CMatrix& operator+= (const CBlock& Mb);
	private:
		void init(bool togpu);
		void init(const std::complex<double>* elem, bool _ongpu);
};
CMatrix exp(double a, const CBlock& mat);
CMatrix exp(const std::complex<double>& a, const CBlock& mat);
CMatrix exph(double a, const CBlock& mat);
CMatrix exp(const CBlock& mat);
CMatrix exph(const CBlock& mat);
CMatrix otimes(const CBlock& Ma, const CBlock& Mb);
};	/* namespace uni10 */
#endif /* CMATRIX_H */

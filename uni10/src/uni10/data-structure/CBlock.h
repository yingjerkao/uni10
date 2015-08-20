/****************************************************************************
*  @file CBlock.h
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
*  @brief Header file for CBlock class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef CBLOCK_H
#define CBLOCK_H
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <cstdint>
#include <vector>
#include <uni10/datatype.hpp>
#include <stdexcept>
#include <complex>

namespace uni10{
class UniTensor;
class Matrix;
class CMatrix;
class CBlock{
	public:
		CBlock();
		CBlock(const CBlock& _b);
    CBlock(size_t _Rnum, size_t _Cnum, bool _diag = false);
		virtual ~CBlock();
		size_t row()const;
		size_t col()const;
		bool isDiag()const;
		bool isOngpu()const;
		size_t elemNum()const;
		std::complex<double> operator[](size_t idx)const;
		std::complex<double> at(size_t i, size_t j)const;
		std::complex<double>* getElem()const;
		CMatrix getDiag()const;
		void save(const std::string& fname)const;
		std::vector<CMatrix> eig()const;
		std::vector<CMatrix> eigh()const;
		std::vector<CMatrix> svd()const;
		/*** qr rq lq ql ***/
		std::vector<CMatrix> qr()const;
		std::vector<CMatrix> rq()const;
		std::vector<CMatrix> ql()const;
		std::vector<CMatrix> lq()const;
		/*******************/
    size_t lanczosEigh(double& E0, CMatrix& psi, size_t max_iter=200, double err_tol = 5E-15)const;
    CMatrix inverse()const;
		std::complex<double> trace()const;
		double norm()const;
		std::complex<double> sum()const;
		friend CMatrix operator*(const CBlock& Ma, const CBlock& Mb);
		friend CMatrix operator*(const CBlock& Ma, const Block& Mb);
		friend CMatrix operator*(const Block& Ma, const CBlock& Mb);
		friend CMatrix operator*(double a, const CBlock& Ma);
		friend CMatrix operator*(const CBlock& Ma, double a);
		friend CMatrix operator*(const std::complex<double>& a, const CBlock& Ma);
		friend CMatrix operator*(const CBlock& Ma, const std::complex<double>& a);
		friend CMatrix operator+(const CBlock& Ma, const CBlock& Mb);
		friend CMatrix operator+(const Block& Ma, const CBlock& Mb);
		friend CMatrix operator+(const CBlock& Ma, const Block& Mb);
		friend bool operator== (const CBlock& m1, const CBlock& m2);
		friend bool operator== (const Block& m1, const CBlock& m2);
		friend bool operator== (const CBlock& m1, const Block& m2);
		friend class UniTensor;
		friend class CUniTensor;
		friend class Block;
		friend class Matrix;
		friend class CMatrix;
		friend std::ostream& operator<< (std::ostream& os, const CBlock& b);
		//friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
	protected:
		std::complex<double>* m_elem;
		size_t Rnum;		//number of rows of the block
		size_t Cnum;		//number of columns of the block
		bool diag;
		bool ongpu;
};
};
#endif /* CBLOCK_H */

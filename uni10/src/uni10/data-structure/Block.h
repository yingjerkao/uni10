/****************************************************************************
*  @file Block.h
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
*  @brief Header file for Block class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef BLOCK_H
#define BLOCK_H
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <uni10/datatype.hpp>
#include <stdexcept>
#include <complex>

namespace uni10{
class UniTensor;
class Matrix;
class CMatrix;
class Block{
	public:
		Block();
		Block(const Block& _b);
    Block(size_t _Rnum, size_t _Cnum, bool _diag = false);
		virtual ~Block();
		size_t row()const;
		size_t col()const;
		bool isDiag()const;
		bool isOngpu()const;
		size_t elemNum()const;
		double operator[](size_t idx)const;
		double at(size_t i, size_t j)const;
		double* getElem()const;
    Matrix getDiag()const;
		void save(const std::string& fname)const;
		std::vector<CMatrix> eig()const;
		std::vector<Matrix> eigh()const;
		std::vector<Matrix> svd()const;
    size_t lanczosEigh(double& E0, Matrix& psi, size_t max_iter=200, double err_tol = 5E-15)const;
    Matrix inverse()const;
		double trace()const;
		double norm()const;
		double sum()const;
		friend Matrix operator*(const Block& Ma, const Block& Mb);
		friend Matrix operator*(double a, const Block& Ma);
		friend Matrix operator*(const Block& Ma, double a);
		friend CMatrix operator*(const std::complex<double>& a, const Block& Ma);
		friend CMatrix operator*(const Block& Ma, const std::complex<double>& a);
		friend Matrix operator+(const Block& Ma, const Block& Mb);
		friend bool operator==(const Block& m1, const Block& m2);
		friend class UniTensor;
		friend class CUniTensor;
		friend class CBlock;
		friend class Matrix;
		friend class CMatrix;
		friend std::ostream& operator<< (std::ostream& os, const Block& b);
		//friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
	protected:
		double* m_elem;
		size_t Rnum;		//number of rows of the block
		size_t Cnum;		//number of columns of the block
		bool diag;
		bool ongpu;
};
};
#endif /* BLOCK_H */

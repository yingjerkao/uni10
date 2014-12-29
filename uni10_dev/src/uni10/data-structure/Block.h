/****************************************************************************
*  @file Block.h
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
*  @brief Header file for Block class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef BLOCK_H
#define BLOCK_H
#ifndef UNI10_DTYPE
#define UNI10_DTYPE double
#endif
#ifndef UNI10_BLOCK
#define UNI10_BLOCK Block
#endif

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
class UNI10_BLOCK{
	public:
		UNI10_BLOCK();
		UNI10_BLOCK(const Block& _b);
    UNI10_BLOCK(size_t _Rnum, size_t _Cnum, bool _diag = false);
		virtual ~UNI10_BLOCK();
		size_t row()const;
		size_t col()const;
		bool isDiag()const;
		bool isOngpu()const;
		size_t elemNum()const;
		UNI10_DTYPE operator[](size_t idx)const;
		UNI10_DTYPE at(size_t i, size_t j)const;
		UNI10_DTYPE* getElem()const;
		void save(const std::string& fname)const;
		std::vector<Matrix> eigh()const;
		std::vector<Matrix> svd()const;
    size_t lanczosEigh(UNI10_DTYPE& E0, Matrix& psi, size_t max_iter=200, UNI10_DTYPE err_tol = 5E-15)const;
    Matrix inverse()const;
		UNI10_DTYPE trace()const;
		UNI10_DTYPE norm()const;
		UNI10_DTYPE sum()const;
		friend Matrix operator* (const Block& Ma, const Block& Mb);
		friend Matrix operator*(const Block& Ma, UNI10_DTYPE a);
		friend Matrix operator*(UNI10_DTYPE a, const Block& Ma);
		friend Matrix operator+(const Block& Ma, const Block& Mb);
		friend bool operator== (const Block& m1, const Block& m2);
		friend class UniTensor;
		friend class Matrix;
		friend std::ostream& operator<< (std::ostream& os, const UNI10_BLOCK& b);
		friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
	protected:
		UNI10_DTYPE* m_elem;
		size_t Rnum;		//number of rows of the block
		size_t Cnum;		//number of columns of the block
		bool diag;
		bool ongpu;
};
};

#ifdef UNI10_BLOCK
#undef UNI10_BLOCK
#endif
#ifdef UNI10_DTYPE
#undef UNI10_DTYPE
#endif
#endif /* BLOCK_H */

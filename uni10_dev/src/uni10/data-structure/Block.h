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
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <cstdint>
#include <uni10/datatype.hpp>
#include <stdexcept>
namespace uni10{
class UniTensor;
class Block{
	public:
		Block();
    Block(size_t _Rnum, size_t _Cnum);
		Block(const Block& _b);
		~Block();
		friend class UniTensor;
		friend std::ostream& operator<< (std::ostream& os, const Block& b);
		friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
		friend bool operator== (const Block& b1, const Block& b2);
	private:
		//Qnum qnum;
		double* m_elem;
		size_t Rnum;		//number of rows of the block
		size_t Cnum;		//number of columns of the block
		size_t m_elemNum;
		bool diag;
		bool ongpu;
		//size_t offset;	//index of the first element of a block element in Tensor
};
};
#endif /* BLOCK_H */


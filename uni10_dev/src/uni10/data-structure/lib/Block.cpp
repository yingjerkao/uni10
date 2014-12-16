/****************************************************************************
*  @file Block.cpp
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
*  @brief Implementation file of Block class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/data-structure/Block.h>
//using namespace uni10::datatype;
namespace uni10{
Block::Block(): Rnum(0), Cnum(0), m_elemNum(0), diag(false), ongpu(false), m_elem(NULL){}
Block::Block(size_t _Rnum, size_t _Cnum): Rnum(_Rnum), Cnum(_Cnum), m_elemNum(_Rnum * _Cnum), diag(false), ongpu(false), m_elem(NULL){}
Block::Block(const Block& _b): Rnum(_b.Rnum), Cnum(_b.Cnum), m_elemNum(_b.m_elemNum), diag(_b.diag), ongpu(_b.ongpu), m_elem(_b.m_elem){}
Block::~Block(){}
std::ostream& operator<< (std::ostream& os, const Block& b){
	os << ": " << b.Rnum << " x " << b.Cnum << " = " << b.Rnum * b.Cnum << " ---\n\n";
	for(size_t r = 0; r < b.Rnum; r++){
		for(size_t c = 0; c < b.Cnum; c++)
			os<< std::setw(7) << std::setprecision(3) << b.m_elem[r * b.Cnum + c];
		os << "\n\n";
	}
	return os;
}
/*
bool operator== (const Block& b1, const Block& b2){
	return (b1.qnum == b2.qnum) && (b1.Rnum == b2.Rnum) && (b1.Cnum == b2.Cnum);
}
*/
};	/* namespace uni10 */

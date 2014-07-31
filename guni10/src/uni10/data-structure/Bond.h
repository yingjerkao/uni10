/****************************************************************************
*  @file Bond.h
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
*  @brief Declaration file for Bond Class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef BOND_H
#define BOND_H

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
#include <uni10/datatype.hpp>

namespace uni10{
enum bondType{
	BD_IN = 1,
	BD_OUT = -1
};
class UniTensor;
class Bond {
	public:
		Bond(){};
		Bond(bondType _type, int dim);
		Bond(bondType, const std::vector<Qnum>& qnums);
		Bond(const Bond& _b);
		void assign(bondType, int dim);
		void assign(bondType, const std::vector<Qnum>& qnums);
		bondType type()const;
		int dim()const;
		friend class UniTensor;
		friend class Node;
		friend std::ostream& operator<< (std::ostream& os, const Bond& b);
		friend bool operator== (const Bond& b1, const Bond& b2);
		void change(bondType tp);
		Bond& combine(const Bond bd);
		friend Bond combine(bondType tp, const std::vector<Bond>& bds);
		friend Bond combine(const std::vector<Bond>& bds);
		std::map<Qnum, int> degeneracy()const;
		std::vector<Qnum> Qlist()const;
		~Bond();
	private:
		void setting(const std::vector<Qnum>& qnums);
		bondType m_type;
		int m_dim;
		std::vector<Qnum>Qnums;	//Quantum numbers
		std::vector<int>Qdegs;	//Degeneracy in each quantum sector
		std::vector<int>offsets;
};
Bond combine(bondType tp, const std::vector<Bond>& bds);
Bond combine(const std::vector<Bond>& bds);
};
#endif /* BOND_H */

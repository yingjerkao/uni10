/****************************************************************************
*  @file Qnum.h
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
*  @brief Header file for Qnum (Quantum Number) class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef QNUM_H
#define QNUM_H
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <stdexcept>
#include <sstream>

namespace uni10{
enum parityType{
	PRT_EVEN = 0,
	PRT_ODD = 1
};
enum parityFType{
	PRTF_EVEN = 0,
	PRTF_ODD = 1
};

class Qnum {
	public:
		Qnum();
		Qnum(int _U1);
		Qnum(int _U1, parityType _prt);
		Qnum(parityFType _prtF);
		Qnum(parityFType _prtF, int _U1);
		Qnum(parityFType _prtF, int _U1, parityType _prt);
		Qnum(const Qnum& _q);
		~Qnum(){};
		int U1()const;
		parityType prt()const;
		parityFType prtF()const;
		void assign(int _U1 = 0, parityType _prt = PRT_EVEN);
		void assign(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);
		//static bool isFermionic();
		static bool isFermionic(){return Fermionic;}
		friend bool operator< (const Qnum& q1, const Qnum& q2);
		friend bool operator<= (const Qnum& q1, const Qnum& q2);
		friend bool operator== (const Qnum& q1, const Qnum& q2);
		friend Qnum operator- (const Qnum& q1);
		friend Qnum operator* (const Qnum& q1, const Qnum& q2);
		friend std::ostream& operator<< (std::ostream& os, const Qnum& q);
		static const int U1_UPB = 1000;	//Upper bound of U1
		static const int U1_LOB = -1000;//Lower bound of U1
	private:
		static bool Fermionic;
		int m_U1;
		parityType m_prt;
		parityFType m_prtF;
};

//};
};
#endif /* QNUMF_H */

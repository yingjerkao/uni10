/****************************************************************************
*  @file CMakeLists.txt
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
*  @brief Main specification file for CMake
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/datatype/Qnum.h>
#include <uni10/tools/uni10_tools.h>

namespace uni10{
bool Qnum::Fermionic = false;
Qnum::Qnum(int _U1, parityType _prt): m_U1(_U1), m_prt(_prt), m_prtF(PRTF_EVEN){
  try{
    if(!(m_U1 < U1_UPB && m_U1 > U1_LOB)){
      std::ostringstream err;
      err<<"U1 is out of range. "<<U1_LOB<<" < U1 < "<<U1_UPB<<".";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Qnum::Qnum(int, parityType):");
  }
}
Qnum::Qnum(parityFType _prtF, int _U1, parityType _prt): m_U1(_U1), m_prt(_prt), m_prtF(_prtF){
  try{
	if(!(m_U1 < U1_UPB && m_U1 > U1_LOB)){
    std::ostringstream err;
    err<<"U1 is out of range. "<<U1_LOB<<" < U1 < "<<U1_UPB<<".";
    throw std::runtime_error(exception_msg(err.str()));
  }
	if(_prtF == PRTF_ODD)
		Fermionic = true;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Qnum::Qnum(parityFType, int, parityType):");
  }
}

Qnum::Qnum(const Qnum& _q):m_U1(_q.m_U1), m_prt(_q.m_prt), m_prtF(_q.m_prtF){}
bool operator< (const Qnum& q1, const Qnum& q2){
	return q1.hash() < q2.hash();
}
bool operator<= (const Qnum& q1, const Qnum& q2){
	return q1.hash() <= q2.hash();
}
bool operator== (const Qnum& q1, const Qnum& q2){
	return (q1.m_U1 == q2.m_U1) && (q1.m_prt == q2.m_prt) && (q1.m_prtF == q2.m_prtF);
}
Qnum operator- (const Qnum& q1){
	Qnum q2(q1.m_prtF, -q1.m_U1, q1.m_prt);
	return q2;
}

Qnum operator* (const Qnum& q1, const Qnum& q2){
	Qnum q3((parityFType)(q1.m_prtF ^ q2.m_prtF), q1.m_U1 + q2.m_U1, (parityType)(q1.m_prt ^ q2.m_prt));
	return q3;
}
std::ostream& operator<< (std::ostream& os, const Qnum& q){
	os << "(U1 = " << std::setprecision(2) << q.m_U1 << ", P = " << std::setprecision(1) << (parityType)q.m_prt << ", " << (parityType)q.m_prtF << ")";
	return os;
}
void Qnum::assign(int _U1, parityType _prt){
	m_U1 = _U1;
	m_prt = _prt;
	m_prtF = PRTF_EVEN;
}
void Qnum::assign(parityFType _prtF, int _U1, parityType _prt){
	m_U1 = _U1;
	m_prt = _prt;
	m_prtF = _prtF;
	if(_prtF == PRTF_ODD)
		Fermionic = true;
}

int Qnum::U1()const{return m_U1;}
parityType Qnum::prt()const{return m_prt;}
parityFType Qnum::prtF()const{return m_prtF;}
long int Qnum::hash()const{ return m_U1 * 4 + m_prt * 2 + m_prtF; }
};

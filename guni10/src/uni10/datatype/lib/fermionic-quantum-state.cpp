/****************************************************************************
*  @file CMakeLists.txt
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
*  @brief Main specification file for CMake
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
/*****************************************************************************
*
* Universal Tensor Network Library
*
* Copyright (C) 2013-2014 
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
* Please direct any enquiry to <development@uni10.org>
*
*****************************************************************************/

#ifndef UNI10_DATATYPE_FERMIONIC_QUANTUM_STATE_CPP
#define UNI10_DATATYPE_FERMIONIC_QUANTUM_STATE_CPP

#include <uni10/datatype/fermionic-quantum-state.h>

namespace uni10 {
namespace datatype {

template<>
fermionic_quantum_state<int,char>
::fermionic_quantum_state
  ( int U1_
  , int prt_
  , int prtF_
  )
  : quantum_state(U1_,prt_)
  , _prtF (prtF_)
{
  assert(_prtF > prt_lbound && _prtF < prt_ubound);
}

template<>
std::ostream&
operator<<
  ( std::ostream & os
  , const fermionic_quantum_state<int,char> & obj
  )
{
  os << "U1 = " << std::setprecision(2) << obj._U1 << "\t"
     << "prt = " << std::setprecision(1) << static_cast<int>(obj._prt) << "\t"
     << "prtF = " << static_cast<int>(obj._prtF) << "\t" 
  ;
  return os;
}

bool operator<
  ( fermionic_quantum_state<int,char> const & q1
  , fermionic_quantum_state<int,char> const & q2
  )
{
  return (q1.U1() * 10 + q1.prt() * 2 + q1.prtF()) < (q2.U1() * 10 + q2.prt() * 2 + q2.prtF());
}

template<>
bool operator<=
  ( fermionic_quantum_state<int,char> const & q1
  , fermionic_quantum_state<int,char> const & q2
  )
{
  return (q1.U1() * 10 + q1.prt() * 2 + q1.prtF()) <= (q2.U1() * 10 + q2.prt() * 2 + q2.prtF());
}

template<>
bool operator==
  ( fermionic_quantum_state<int,char> const & q1
  , fermionic_quantum_state<int,char> const & q2
  )
{
  return (q1.U1() == q2.U1()) && (q1.prt() == q2.prt()) && (q1.prtF() == q2.prtF());
}

template<>
fermionic_quantum_state<int,char>
operator-
  ( fermionic_quantum_state<int,char> const & q
  )
{
  return fermionic_quantum_state<int,char>(-q.U1(), q.prt(), q.prtF());
}

template<>
fermionic_quantum_state<int,char>
operator*
  ( fermionic_quantum_state<int,char> const & q1
  , fermionic_quantum_state<int,char> const & q2
  )
{
  return fermionic_quantum_state<int,char>(q1.U1() + q2.U1(), q1.prt() ^ q2.prt(), q1.prtF() ^ q2.prtF());
}

} // ending namespace datatype
} // ending namespace uni10

#endif

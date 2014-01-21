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

#ifndef UNI10_DATATYPE_QUANTUM_STATE_CPP
#define UNI10_DATATYPE_QUANTUM_STATE_CPP

#include <uni10/datatype/quantum-state.h>

namespace uni10 {
namespace datatype {

template<>
quantum_state<int,char>
::quantum_state
  ( int U1_ 
  , int prt_
  )
  : _U1  (U1_)
  , _prt (prt_)
{
  assert(  _U1 > U1_lbound && _U1 < U1_ubound
        && _prt > prt_lbound && _prt < prt_ubound
        );
}

template<>
std::ostream&
operator<<
  ( std::ostream & os
  , const quantum_state<int,char> & obj
  )
{
  os << "U1 = " << std::setprecision(2) << obj._U1 << "\t"
     << "prt = " << std::setprecision(1) << static_cast<int>(obj._prt) << "\t"
  ;
  return os;
}

template<>
bool operator<
  ( quantum_state<int,char> const & q1
  , quantum_state<int,char> const & q2
  ) 
{
	return (q1.U1() * 10 + q1.prt()) < (q2.U1() * 10 + q2.prt());
}

template<>
bool operator<= 
  ( quantum_state<int,char> const & q1
  , quantum_state<int,char> const & q2
  )
{
  return (q1.U1() * 10 + q1.prt()) <= (q2.U1() * 10 + q2.prt());
}

template<>
bool operator== 
  ( quantum_state<int,char> const & q1
  , quantum_state<int,char> const & q2
  )
{
	return (q1.U1() == q2.U1()) && (q1.prt() == q2.prt());
}

template<>
quantum_state<int,char>
operator- 
  ( quantum_state<int,char> const & q
  )
{
	return quantum_state<int,char>(-q.U1(), q.prt());
}

template<>
quantum_state<int,char> 
operator* 
  ( quantum_state<int,char> const & q1
  , quantum_state<int,char> const & q2
  )
{
	return quantum_state<int,char>(q1.U1() + q2.U1(), q1.prt() ^ q2.prt());
}


} // ending namespace datatype
} // ending namespace uni10

#endif

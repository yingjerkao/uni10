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

#ifndef UNI10_DATATYPE_QNUM_CPP
#define UNI10_DATATYPE_QNUM_CPP

#include <uni10/datatype/Qnum.h>

namespace uni10 {
namespace datatype {

template<>
Qnum<int,char>
::Qnum
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
void 
Qnum<int,char>
::set
  ( int U1_
  , int prt_
  )
{
  Qnum tmp(U1_,prt_);  
  std::swap(*this,tmp);
}

template<>
std::ostream&
operator<<
  ( std::ostream & os
  , const Qnum<int,char> & obj
  )
{
  os << "(U1 = " 
     << std::setprecision(2) << obj._U1 
     << ", P = " 
     << std::setprecision(1) << static_cast<int>(obj._prt) 
     << ")"
  ;
  return os;
}

bool operator<
  ( Qnum<int,char> const & q1
  , Qnum<int,char> const & q2
  )
{
	return (q1.U1() * 10 + q1.prt()) < (q2.U1() * 10 + q2.prt());
}

bool operator<= 
  ( Qnum<int,char> const & q1
  , Qnum<int,char> const & q2
  )
{
  return (q1.U1() * 10 + q1.prt()) <= (q2.U1() * 10 + q2.prt());
}

bool operator== 
  ( Qnum<int,char> const & q1
  , Qnum<int,char> const & q2
  )
{
	return (q1.U1() == q2.U1()) && (q1.prt() == q2.prt());
}

Qnum<int,char>
operator- 
  ( Qnum<int,char> const & q
  )
{
	return Qnum<int,char>(-q.U1(), q.prt());
}

Qnum<int,char> 
operator* 
  ( Qnum<int,char> const & q1
  , Qnum<int,char> const & q2
  )
{
	return Qnum<int,char>(q1.U1() + q2.U1(), q1.prt() ^ q2.prt());
}


} // ending namespace datatype
} // ending namespace uni10

#endif

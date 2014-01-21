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

#ifndef UNI10_DATATYPE_QUANTUM_STATE_H
#define UNI10_DATATYPE_QUANTUM_STATE_H

#include <iostream>
#include <iomanip>
#include <cassert>

namespace uni10 {
namespace datatype {

template <class I, class S>
class quantum_state 
{
public:
  // typedefs
  typedef I  int_type;    
  typedef S  short_type;

  // Non-standard constructor
  quantum_state
    ( int_type U1_   = int_type()
    , int_type prt_  = int_type()
    );

  // I/O
  inline int_type U1()  const  { return _U1;  }
  inline int_type prt() const  { return _prt; } 

  inline int_type get_U1()  const  { return U1();  }   // to be depreciated...
  inline int_type get_prt() const  { return prt(); }   // to be depreciated...

  inline void set_U1  (int_type U1_)    { _U1  = U1_;  }
  inline void set_prt (int_type prt_)   { _prt = prt_; }
  inline void set     (int_type U1_, int_type prt_)   { set_U1(U1_); set_prt(prt_); } 

  template <class I1, class S1>
  friend bool operator<  (quantum_state<I1,S1> const & q1, quantum_state<I1,S1> const & q2);
  template <class I1, class S1>
  friend bool operator<= (quantum_state<I1,S1> const & q1, quantum_state<I1,S1> const & q2);
  template <class I1, class S1>
  friend bool operator== (quantum_state<I1,S1> const & q1, quantum_state<I1,S1> const & q2);

  template <class I1, class S1>
  friend quantum_state<I1,S1> operator- (quantum_state<I1,S1> const & q);
  template <class I1, class S1>
  friend quantum_state<I1,S1> operator* (quantum_state<I1,S1> const & q1, quantum_state<I1,S1> const & q2);

  template <class I1, class S1>
  friend std::ostream& operator<< (std::ostream& os, const quantum_state<I1,S1> & obj);

protected:
  int_type    _U1;
  short_type  _prt;

  static const int U1_ubound  = 100; 
  static const int U1_lbound  = -100;
  static const int prt_ubound = 2;
  static const int prt_lbound = -1;
};

} // ending namespace datatype
} // ending namespace uni10

#endif

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

#ifndef UNI10_DATATYPE_FERMIONIC_QUANTUM_STATE_H
#define UNI10_DATATYPE_FERMIONIC_QUANTUM_STATE_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <uni10/datatype/quantum-state.h>

namespace uni10 {
namespace datatype {

#define FERMIONIC 1  // hmm... we should not do this...

template <class I, class S>
class fermionic_quantum_state 
  : public quantum_state<I,S>
{
public:
  typedef typename quantum_state<I,S>::int_type    int_type;
  typedef typename quantum_state<I,S>::short_type  short_type;

  // Non-standard constructor
  fermionic_quantum_state
    ( int_type U1_   = int_type()
    , int_type prt_  = int_type()
    , int_type prtF_ = int_type()
    );

  inline int_type prtF()     const  { return _prtF; }
  inline int_type get_prtF() const  { return prtF(); }   // to be depreciated...

  inline void set_prtF (int_type prtF_)   { _prtF = prtF_; }
  inline void set      (int_type U1_, int_type prt_)                   { this->set_U1(U1_); this->set_prt(prt_); } 
  inline void set      (int_type U1_, int_type prt_, int_type prtF_)   { this->set_U1(U1_); this->set_prt(prt_); set_prtF(prtF_); }

  template <class I1, class S1>
  friend bool operator<  (fermionic_quantum_state<I1,S1> const & q1, fermionic_quantum_state<I1,S1> const & q2);
  template <class I1, class S1>
  friend bool operator<= (fermionic_quantum_state<I1,S1> const & q1, fermionic_quantum_state<I1,S1> const & q2);
  template <class I1, class S1>
  friend bool operator== (fermionic_quantum_state<I1,S1> const & q1, fermionic_quantum_state<I1,S1> const & q2);

  template <class I1, class S1>
  friend fermionic_quantum_state<I1,S1> operator- (fermionic_quantum_state<I1,S1> const & q);
  template <class I1, class S1>
  friend fermionic_quantum_state<I1,S1> operator* (fermionic_quantum_state<I1,S1> const & q1, fermionic_quantum_state<I1,S1> const & q2);

  template <class I1, class S1>
  friend std::ostream& operator<< (std::ostream& os, const fermionic_quantum_state<I1,S1> & obj);

private:
  short_type  _prtF;

};

} // ending namespace datatype
} // ending namespace uni10

#endif

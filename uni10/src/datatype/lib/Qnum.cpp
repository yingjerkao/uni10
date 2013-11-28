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

const int Qnum::U1_UPB = 100;      // this should be made nicer...
const int Qnum::U1_LOB = -100;      // ....
const int Qnum::ptr_UPB = 2;
const int Qnum::ptr_LOB = -1;


bool operator< (const Qnum_t& q1, const Qnum_t& q2){
	return (q1.U1 * 10 + q1.prt) < (q2.U1 * 10 + q2.prt);
}
bool operator<= (const Qnum_t& q1, const Qnum_t& q2){
	return (q1.U1 * 10 + q1.prt) <= (q2.U1 * 10 + q2.prt);
}
bool operator== (const Qnum_t& q1, const Qnum_t& q2){
	return (q1.U1 == q2.U1) && (q1.prt == q2.prt);
}

Qnum_t operator- (const Qnum_t& q1){
	Qnum_t q2(-q1.U1, q1.prt);
	return q2;
}

Qnum_t operator* (const Qnum_t& q1, const Qnum_t& q2){
	Qnum_t q3(q1.U1 + q2.U1, q1.prt ^ q2.prt);
	return q3;
}

ostream& operator<< (ostream& os, const Qnum_t& q){
	os << "(U1 = " << setprecision(2) << q.U1 << ", P = " << setprecision(1) << (int)q.prt << ")";
	return os;
}

void Qnum_t::set(int _U1, int _prt){
	U1 = _U1;	
	prt = _prt;
}

} // ending namespace datatype
} // ending namespace uni10

#endif

#include "Qnum.h"
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
	Qnum_t q2;
	q2.U1 = -q1.U1;
	q2.prt = q1.prt;
	return q2;
}

Qnum_t operator* (const Qnum_t& q1, const Qnum_t& q2){
	Qnum_t q3;
	q3.U1 = q1.U1 + q2.U1;
	q3.prt = q1.prt ^ q2.prt;
	return q3;
}

ostream& operator<< (ostream& os, const Qnum_t& q){
	os << "(U1 = " << setprecision(2) << q.U1 << ", P = " << setprecision(1) << q.prt << ")";
	return os;
}


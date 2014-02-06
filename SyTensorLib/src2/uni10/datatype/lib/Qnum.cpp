#include <uni10/datatype/Qnum.h>
//using namespace uni10::datatype;
//namespace uni10{
//namespace datatype{
bool operator< (const Qnum& q1, const Qnum& q2){
	return ((q1.U1 * 10) + (q1.prt * 2) + q1.prtF) < ((q2.U1 * 10) + (q2.prt * 2) + q2.prtF);
}
bool operator<= (const Qnum& q1, const Qnum& q2){
	return ((q1.U1 * 10) + (q1.prt * 2) + q1.prtF) <= ((q2.U1 * 10) + (q2.prt * 2) + q2.prtF);
}
bool operator== (const Qnum& q1, const Qnum& q2){
	return (q1.U1 == q2.U1) && (q1.prt == q2.prt) && (q1.prtF == q2.prtF);
}

Qnum operator- (const Qnum& q1){
	Qnum q2(-q1.U1, q1.prt, q1.prtF);
	return q2;
}

Qnum operator* (const Qnum& q1, const Qnum& q2){
	Qnum q3(q1.U1 + q2.U1, q1.prt ^ q2.prt, q1.prtF ^ q2.prtF);
	return q3;
}

std::ostream& operator<< (std::ostream& os, const Qnum& q){
	os << "(U1 = " << std::setprecision(2) << q.U1 << ", P = " << std::setprecision(1) << (int)q.prt << ", " << (int)q.prtF << ")";
	return os;
}

void Qnum::set(int _U1, int _prt, int _prtF){
	U1 = _U1;	
	prt = _prt;
	prtF = _prtF;
}
//};
//};

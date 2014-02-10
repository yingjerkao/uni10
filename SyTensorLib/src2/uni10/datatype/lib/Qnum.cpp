#include <uni10/datatype/Qnum.h>
//using namespace uni10::datatype;
namespace uni10{
Qnum::Qnum(): m_U1(0), m_prt(PRT_EVEN), m_prtF(PRTF_EVEN){}
Qnum::Qnum(int _U1): m_U1(_U1), m_prt(PRT_EVEN), m_prtF(PRTF_EVEN){
	assert(m_U1 < U1_UPB && m_U1 > U1_LOB);
}
Qnum::Qnum(int _U1, parityType _prt): m_U1(_U1), m_prt(_prt), m_prtF(PRTF_EVEN){
	assert(m_U1 < U1_UPB && m_U1 > U1_LOB && m_prt < PRT_UPB && m_prt > PRT_LOB);
}
Qnum::Qnum(parityFType _prtF): m_prtF(_prtF), m_U1(0), m_prt(PRT_EVEN){
	assert(m_prtF < PRT_UPB && m_prtF > PRT_LOB);
}
Qnum::Qnum(parityFType _prtF, int _U1): m_prtF(_prtF), m_U1(_U1), m_prt(PRT_EVEN){
	assert(m_prtF < PRT_UPB && m_prtF > PRT_LOB && m_U1 < U1_UPB && m_U1 > U1_LOB);
}
Qnum::Qnum(parityFType _prtF, int _U1, parityType _prt): m_U1(_U1), m_prt(_prt), m_prtF(_prtF){
	assert(m_prtF < PRT_UPB && m_prtF > PRT_LOB && m_U1 < U1_UPB && m_U1 > U1_LOB && m_prt < PRT_UPB && m_prt > PRT_LOB);
}
Qnum::Qnum(const Qnum& _q):m_U1(_q.m_U1), m_prt(_q.m_prt), m_prtF(_q.m_prtF){}
bool operator< (const Qnum& q1, const Qnum& q2){
	return ((q1.m_U1 * 10) + (q1.m_prt * 2) + q1.m_prtF) < ((q2.m_U1 * 10) + (q2.m_prt * 2) + q2.m_prtF);
}
bool operator<= (const Qnum& q1, const Qnum& q2){
	return ((q1.m_U1 * 10) + (q1.m_prt * 2) + q1.m_prtF) <= ((q2.m_U1 * 10) + (q2.m_prt * 2) + q2.m_prtF);
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
}

int Qnum::U1()const{return m_U1;}
parityType Qnum::prt()const{return m_prt;}
parityFType Qnum::prtF()const{return m_prtF;}
};

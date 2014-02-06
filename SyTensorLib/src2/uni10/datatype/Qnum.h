#ifndef QNUM_H
#define QNUM_H
#include <iostream>
#include <iomanip>
#include <assert.h>

#define FERMIONIC 1
//namespace uni10{
//namespace datatype{

class Qnum {
	public:
		Qnum(): U1(0), prt(0), prtF(0){}
		Qnum(int _U1): U1(_U1), prt(0), prtF(0){
			assert(U1 < U1_UPB && U1 > U1_LOB);
		}
		Qnum(int _U1, int _prt): U1(_U1), prt(_prt), prtF(0){
			assert(U1 < U1_UPB && U1 > U1_LOB && prt < prt_UPB && prt > prt_LOB);
		}
		Qnum(int _U1, int _prt, int _prtF): U1(_U1), prt(_prt), prtF(_prtF){
			assert(U1 < U1_UPB && U1 > U1_LOB && prt < prt_UPB && prt > prt_LOB && prtF < prt_UPB && prtF > prt_LOB);
		}
		Qnum(const Qnum& _q):U1(_q.U1), prt(_q.prt), prtF(_q.prtF){}
		~Qnum(){};
		int getU1()const{return U1;}
		int getPrt()const{return prt;}
		int getPrtF()const{return prtF;}
		void set(int _U1 = 0, int _prt = 0, int _prtF = 0);
		friend bool operator< (const Qnum& q1, const Qnum& q2);
		friend bool operator<= (const Qnum& q1, const Qnum& q2);
		friend bool operator== (const Qnum& q1, const Qnum& q2);
		friend Qnum operator- (const Qnum& q1);
		friend Qnum operator* (const Qnum& q1, const Qnum& q2);
		friend std::ostream& operator<< (std::ostream& os, const Qnum& q);
	private:
		static const int U1_UPB = 100;	//Upper bound of U1
		static const int U1_LOB = -100;//Lower bound of U1
		static const int prt_UPB = 2;  //Upper bound of prt
		static const int prt_LOB = -1; //Lower bound of prt
		int U1;
		char prt;
		char prtF;
};

//};
//};
#endif /* QNUMF_H */

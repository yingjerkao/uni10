#pragma once
#include <iostream>
#include <iomanip>
#include <assert.h>
using namespace std;
#define FERMIONIC 1

const int U1_UPB = 100;	//Upper bound of U1
const int U1_LOB = -100;//Lower bound of U1
const int prt_UPB = 2;  //Upper bound of prt
const int prt_LOB = -1; //Lower bound of prt

class Qnum_t {
	public:
		Qnum_t(): U1(0), prt(0), prtF(0){}
		Qnum_t(int _U1): U1(_U1), prt(0), prtF(0){
			assert(U1 < U1_UPB && U1 > U1_LOB);
		}
		Qnum_t(int _U1, int _prt): U1(_U1), prt(_prt), prtF(0){
			assert(U1 < U1_UPB && U1 > U1_LOB && prt < prt_UPB && prt > prt_LOB);
		}
		Qnum_t(int _U1, int _prt, int _prtF): U1(_U1), prt(_prt), prtF(_prtF){
			assert(U1 < U1_UPB && U1 > U1_LOB && prt < prt_UPB && prt > prt_LOB && prtF < prt_UPB && prtF > prt_LOB);
		}
		Qnum_t(const Qnum_t& _q):U1(_q.U1), prt(_q.prt), prtF(_q.prtF){}
		~Qnum_t(){};
		int getU1()const{return U1;}
		int getPrt()const{return prt;}
		int getPrtF()const{return prtF;}
		void set(int _U1 = 0, int _prt = 0, int _prtF = 0);
		friend bool operator< (const Qnum_t& q1, const Qnum_t& q2);
		friend bool operator<= (const Qnum_t& q1, const Qnum_t& q2);
		friend bool operator== (const Qnum_t& q1, const Qnum_t& q2);
		friend Qnum_t operator- (const Qnum_t& q1);
		friend Qnum_t operator* (const Qnum_t& q1, const Qnum_t& q2);
		friend ostream& operator<< (ostream& os, const Qnum_t& q);
	private:
		int U1;
		char prt;
		char prtF;
};


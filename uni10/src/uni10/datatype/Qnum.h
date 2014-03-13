#ifndef QNUM_H
#define QNUM_H
#include <iostream>
#include <iomanip>
#include <assert.h>

namespace uni10{

enum parityType{
	PRT_EVEN = 0,
	PRT_ODD = 1
};
enum parityFType{
	PRTF_EVEN = 0,
	PRTF_ODD = 1
};

class Qnum {
	public:
		Qnum();
		Qnum(int _U1);
		Qnum(int _U1, parityType _prt);
		Qnum(parityFType _prtF);
		Qnum(parityFType _prtF, int _U1);
		Qnum(parityFType _prtF, int _U1, parityType _prt);
		Qnum(const Qnum& _q);
		~Qnum(){};
		int U1()const;
		parityType prt()const;
		parityFType prtF()const;
		void assign(int _U1 = 0, parityType _prt = PRT_EVEN);
		void assign(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);
		//static bool isFermionic();
		static bool isFermionic(){return Fermionic;}
		friend bool operator< (const Qnum& q1, const Qnum& q2);
		friend bool operator<= (const Qnum& q1, const Qnum& q2);
		friend bool operator== (const Qnum& q1, const Qnum& q2);
		friend Qnum operator- (const Qnum& q1);
		friend Qnum operator* (const Qnum& q1, const Qnum& q2);
		friend std::ostream& operator<< (std::ostream& os, const Qnum& q);
	private:
		static const int U1_UPB = 100;	//Upper bound of U1
		static const int U1_LOB = -100;//Lower bound of U1
		static const int PRT_UPB = 2;  //Upper bound of prt
		static const int PRT_LOB = -1; //Lower bound of prt
		static bool Fermionic;
		int m_U1;
		parityType m_prt;
		parityFType m_prtF;
};

//};
};
#endif /* QNUMF_H */

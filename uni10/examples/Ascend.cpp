#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
//#include "elemCmp.cpp"
#define FERMION 1


int main(){
	Qnum q10(1, PRT_EVEN);
	Qnum q_10(-1, PRT_EVEN);
	Qnum q30(3, PRT_EVEN);
#ifdef FERMION
	Qnum q_11(PRTF_ODD, -1, PRT_EVEN);
	Qnum q11(PRTF_ODD, 1, PRT_EVEN);
	Qnum q_31(PRTF_ODD, -3, PRT_EVEN);
#else
	Qnum q_11(-1, PRT_ODD);
	Qnum q11(1, PRT_ODD);
	Qnum q_31(-3, PRT_ODD);
#endif
	vector<Bond> bonds;
	vector<Qnum> qnums;
	qnums.push_back(q10);qnums.push_back(q_11);
	Bond bdr(BD_IN, qnums);
	Bond bdc(BD_OUT, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

	UniTensor H0(bonds, "Ob");
	H0.addRawElem(H_elem);

	UniTensor U(bonds, "U");
	U.orthoRand();
	UniTensor UT = U;
	UT.transpose();
	
	bonds.clear();
	vector<Qnum> qnums1;
	qnums1.push_back(q30); qnums1.push_back(q11); qnums1.push_back(q11); qnums1.push_back(q11);
	qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_31);
	Bond bdr1(BD_IN, qnums1);
	Bond bdc1(BD_OUT, qnums1);
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	UniTensor W1(bonds, "W1");
	W1.orthoRand();
	UniTensor W1T = W1;
	W1T.transpose();
	UniTensor W2(bonds, "W2");
	W2.orthoRand();
	UniTensor W2T = W2;
	W2T.transpose();
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdr1);
	bonds.push_back(bdc1);
	bonds.push_back(bdc1);
	UniTensor Rho(bonds, "Rho");
	Rho.orthoRand();

	vector<UniTensor*> tens;
	tens.push_back(&W1);
	tens.push_back(&W2);
	tens.push_back(&U);
	tens.push_back(&H0);
	tens.push_back(&UT);
	tens.push_back(&W1T);
	tens.push_back(&W2T);
	tens.push_back(&Rho);
	UniTensor H1, H2;
	Network asdL("AscendL", tens);
	H1 = asdL.launch();
	asdL.putTensor(0, &W1);
	asdL.putTensor(1, &W2);
	asdL.putTensor(2, &U);
	asdL.putTensor(3, &H0);
	asdL.putTensor(4, &UT);
	asdL.putTensor(5, &W1T);
	asdL.putTensor(6, &W2T);
	asdL.putTensor(7, &Rho);
	H2 = asdL.launch();
	cout<<asdL;
	cout<<H1;
	cout<<H2;
	//H1.printRawElem();
}


#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include "MERA_Operator2.cpp"
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
	cout<<W1;


	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdr1);
	bonds.push_back(bdc1);
	bonds.push_back(bdc1);
	UniTensor Rho(bonds, "Rho");
	Rho.orthoRand();
	//cout<<Rho;
	
	H0.rawElem().save("H0_elem1");
	W1.rawElem().save("W1_elem1");
	U.rawElem().save("U_elem1");
	Rho.rawElem().save("Rho_elem1");
	
	
	// write out tensors before contraction.
	
	int label_H0[] = {1, 4, 2, 5};
	H0.addLabel(label_H0);
	H0.save("tenH0");
	U.save("tenU");
	UT.save("tenUT");
	W1.save("tenW1");
	W2.save("tenW2");
	W1T.save("tenW1T");
	W2T.save("tenW2T");
	Rho.save("tenRho");
	// Network replaceWith()
	UniTensor H1;
	Network asdL("../Diagrams/AscendL");
	Network asdC("../Diagrams/AscendC");
	Network asdR("../Diagrams/AscendR");
	H1 = Ascend(H0, W1, U, asdL, asdC, asdR);
	//cout<<asdL;
	//cout<<H1;
	//H1.save("tenH1");
	H1.printRawElem();
	H0.check();
	//Check whether tensors are modified or not after contractions
	UniTensor H0_tmp("tenH0");
	UniTensor U_tmp("tenU");
	UniTensor UT_tmp("tenUT");
	UniTensor W1_tmp("tenW1");
	UniTensor W1T_tmp("tenW1T");
	UniTensor W2_tmp("tenW2");
	UniTensor W2T_tmp("tenW2T");
	UniTensor Rho_tmp("tenRho");
	UniTensor H1_tmp("tenH1");
	
	assert(H0.elemCmp(H0_tmp));
	assert(U.elemCmp(U_tmp));
	assert(UT.elemCmp(UT_tmp));
	assert(W1.elemCmp(W1_tmp));
	assert(W1T.elemCmp(W1T_tmp));
	assert(W2.elemCmp(W2_tmp));
	assert(W2T.elemCmp(W2T_tmp));
	assert(Rho.elemCmp(Rho_tmp));
	assert(H1.elemCmp(H1_tmp));
}


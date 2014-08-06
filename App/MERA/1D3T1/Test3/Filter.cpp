#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#include "MERA_Operator.cpp"
#define FERMION 1



int main(){
	Qnum_t q10(1, 0);
	Qnum_t q_10(-1, 0);
	Qnum_t q30(3, 0);
#ifdef FERMION
	Qnum_t q_11(-1, 0, 1);
	Qnum_t q11(1, 0, 1);
	Qnum_t q_31(-3, 0, 1);
#else
	Qnum_t q_11(-1, 1);
	Qnum_t q11(1, 1);
	Qnum_t q_31(-3, 1);
#endif
	vector<Bond_t> bonds;
	vector<Qnum_t> qnums;
	qnums.push_back(q10);qnums.push_back(q_11);
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

	SyTensor_t H0(bonds, "Ob");
	H0.addRawElem(H_elem);

	//printRawElem(H0, "H0_elem");
	
	// write out tensors before contraction.
	
	int label_H0[] = {1, 4, 2, 5};
	H0.addLabel(label_H0);
	H0.save("tenH0");

	// Network replaceWith()
	SyTensor_t newRho;

	vector<Qnum_t> Qnums = H0.qnums();
	/*
	for(int i = 0;i < Qnums.size(); i++)
		cout << Qnums[i] << endl;
	*/
	Qnum_t Qnum = Qnums[1];
	newRho = Filter(H0,Qnum);
	
	//cout << newRho;
	newRho.save("tennewRho");
	newRho.printRawElem();
	H0.check();
	//Check whether tensors are modified or not after contractions
	SyTensor_t H0_tmp("tenH0");
	SyTensor_t newRho_tmp("tennewRho");
	
	assert(H0.elemCmp(H0_tmp));
	assert(newRho.elemCmp(newRho_tmp));
}


#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#define FERMION 1

bool elemCmp(SyTensor_t& ten1, SyTensor_t& ten2);
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

	SyTensor_t U(bonds, "U");
	U.orthoRand();
	SyTensor_t UT = U;
	UT.transpose();
	
	bonds.clear();
	vector<Qnum_t> qnums1;
	qnums1.push_back(q30); qnums1.push_back(q11); qnums1.push_back(q11); qnums1.push_back(q11);
	qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_31);
	Bond_t bdr1(BD_ROW, qnums1);
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t W1(bonds, "W1");
	W1.orthoRand();
	SyTensor_t W1T = W1;
	W1T.transpose();

	SyTensor_t W2(bonds, "W2");
	W2.orthoRand();
	SyTensor_t W2T = W2;
	W2T.transpose();
	SyTensor_t W2F = W2;
	W2F.orthoRand();
	
	SyTensor_t H1;
	
	// raw contractions
	/*
	int label_W1[] = {-1, 0, 1, 3};
	W1.addLabel(label_W1);
	int label_W1T[] = {0, 2, 6, -3};
	W1T.addLabel(label_W1T);
	int label_U[] = {3, 7, 4, 8};
	U.addLabel(label_U);
	int label_H0[] = {1, 4, 2, 5};
	H0.addLabel(label_H0);
	int label_UT[] = {5, 8, 6, 9};
	UT.addLabel(label_UT);
	int label_W2[] = {-2, 7, 10, 11};
	W2.addLabel(label_W2);
	int label_W2T[] = {9, 10, 11, -4};
	W2T.addLabel(label_W2T);

	H0.save("tenH0");
	U.save("tenU");
	UT.save("tenUT");
	W1.save("tenW1");
	W2.save("tenW2");
	W1T.save("tenW1T");
	W2T.save("tenW2T");

	int label_tmp[] = {-2, 10, 11, 3, 4, 8};
	H1 = W1 * W2;
	H1 *= U;
	H1 *= H0;
	H1 *= UT;
	H1 *= W1T;
	H1 *= W2T;
	*/
	vector<SyTensor_t*> tens;
	tens.push_back(&W1);
	tens.push_back(&W2);
	tens.push_back(&U);
	tens.push_back(&H0);
	tens.push_back(&UT);
	tens.push_back(&W1T);
	tens.push_back(&W2T);
	Network_t net3("AscendL", tens);
	cout<< net3;
	H1 = net3.launch();
	int label_out[] = {-1, -2, -3, -4};
	H1.reshape(label_out, 2);
	//H1.save("tenH1");
	
	cout<<H1;

	//Sanity Check
	SyTensor_t H0_tmp("tenH0");
	SyTensor_t U_tmp("tenU");
	SyTensor_t UT_tmp("tenUT");
	SyTensor_t W1_tmp("tenW1");
	SyTensor_t W1T_tmp("tenW1T");
	SyTensor_t W2_tmp("tenW2");
	SyTensor_t W2T_tmp("tenW2T");
	SyTensor_t H1_tmp("tenH1");
	
	assert(elemCmp(H0, H0_tmp));
	assert(elemCmp(U, U_tmp));
	assert(elemCmp(UT, UT_tmp));
	assert(elemCmp(W1, W1_tmp));
	assert(elemCmp(W1T, W1T_tmp));
	assert(elemCmp(W2, W2_tmp));
	assert(elemCmp(W2T, W2T_tmp));
	assert(elemCmp(H1, H1_tmp));

	H0.check();
}

bool elemCmp(SyTensor_t& ten1, SyTensor_t& ten2){
	int n1 = ten1.getElemNum();
	int n2 = ten2.getElemNum();
	double diff;
	if(n1 == n2){
		for(int i = 0; i < n1; i++){
			diff = fabs(ten1.elem[i] - ten2.elem[i]);
			if(diff > 1E-6)
				return false;
		}
	}
	else
		return false;
	return true;
}


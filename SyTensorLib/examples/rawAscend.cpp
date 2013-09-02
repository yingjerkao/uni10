#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"

int main(){
	Qnum_t q10(1, 0);
	Qnum_t q_10(-1, 0);
	Qnum_t q30(3, 0);
	Qnum_t q_30(-3, 0);
	vector<Bond_t> bonds;		
	vector<Qnum_t> qnums;
	qnums.push_back(q10);qnums.push_back(q_10);
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
	int label_H0[] = {2, 5, 3, 6};
	SyTensor_t H0(bonds, label_H0, "H0");
	H0.addRawElem(H_elem);

	SyTensor_t U(bonds, "U");
	U.orthoRand();
	SyTensor_t UT = U;
	UT.transpose();
	int label_U[] = {4, 8, 5, 9};
	U.addLabel(label_U);
	int label_UT[] = {6, 9, 7, 10};
	UT.addLabel(label_UT);
	
	bonds.clear();
	vector<Qnum_t> qnums1;
	//qnums1.push_back(q30);qnums1.push_back(q10);
	//qnums1.push_back(q_10);qnums1.push_back(q_30);
	qnums1.push_back(q30);qnums1.push_back(q10);qnums1.push_back(q10);qnums1.push_back(q10);
	qnums1.push_back(q_10);qnums1.push_back(q_10);qnums1.push_back(q_10);qnums1.push_back(q_30);
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
	int label_W1[] = {-1, 1, 2, 4};
	W1.addLabel(label_W1);
	int label_W1T[] = {1, 3, 7, -3};
	W1T.addLabel(label_W1T);

	SyTensor_t W2(bonds, "W2");
	W2.orthoRand();
	SyTensor_t W2T = W2;
	W2T.transpose();
	int label_W2[] = {-2, 8, 11, 12};
	W2.addLabel(label_W2);
	int label_W2T[] = {10, 11, 12, -4};
	W2T.addLabel(label_W2T);
	U.transpose();

	
	SyTensor_t H1 = W1T * W1;		
	SyTensor_t tmp = UT * H0;
	tmp *= U;
	tmp *= W2;
	tmp *= W2T;
	H1 *= tmp;
	
	int label_out[] = {-1, -2, -3, -4};
	H1.reshape(label_out, 2);
	cout<<H1;
	printRawElem(H1);
	H0.check();
}

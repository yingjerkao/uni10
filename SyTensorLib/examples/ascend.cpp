#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
//#define FERMION 1

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

	H1 = W1 * W2;		
	H1 *= U;
	H1 *= H0;
	H1 *= UT;
	H1 *= W1T;
	H1 *= W2T;
	*/

	// END raw contraction
	
	// Network replaceWith()
	/*
	Network_t net2("AscendL");
	net2.replaceWith(0, &W1);
	net2.replaceWith(1, &W2);
	net2.replaceWith(2, &U);
	net2.replaceWith(3, &H0);
	net2.replaceWith(4, &UT);
	net2.replaceWith(5, &W1T);
	net2.replaceWith(6, &W2T);
	H1 = net2.launch();
	cout<<net2;
	*/
	//  END replaceWith()
	// Network Construct()
	vector<SyTensor_t*> tens;
	tens.push_back(&W1);
	tens.push_back(&W1T);
	tens.push_back(&U);
	tens.push_back(&H0);
	tens.push_back(&UT);
	tens.push_back(&W2);
	tens.push_back(&W2T);
	Network_t net3("AscendL", tens);
	H1 = net3.launch();
	// END Network Construct()
	//int label_out[] = {-1, -2, -3, -4};
	int label_out[] = {-1, -2, -3, -4};
	H1.reshape(label_out, 2);
	cout<<H1;
	//printRawElem(H1);
	//cout<<W1T;

	H0.check();
}


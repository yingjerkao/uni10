#include <iostream>
#include <assert.h>
#include <map>
#include <uni10.hpp>
#define FERMION 1

using namespace std;

int main(){
	uni10::Qnum q10(1,uni10::PRT_EVEN);
	uni10::Qnum q_10(-1, uni10::PRT_EVEN);
	uni10::Qnum q30(3, uni10::PRT_EVEN);
	
#ifdef FERMION
	uni10::Qnum q_11(uni10::PRTF_EVEN,-1, uni10::PRT_ODD);
	uni10::Qnum q11(uni10::PRTF_EVEN,1, uni10::PRT_ODD);
	uni10::Qnum q_31(uni10::PRTF_EVEN,-3, uni10::PRT_ODD);
#else
	uni10::Qnum q_11(-1, uni10::PRT_ODD);
	uni10::Qnum q11(1, uni10::PRT_ODD);
	uni10::Qnum q_31(-3, uni10::PRT_ODD);
#endif
	vector<uni10::Bond> bonds;
	vector<uni10::Qnum> qnums;
	qnums.push_back(q10);qnums.push_back(q_11);
	uni10::Bond bdr(uni10::BD_IN, qnums);
	uni10::Bond bdc(uni10::BD_OUT, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

	uni10::UniTensor H0(bonds, "Ob");
	H0.addRawElem(H_elem);

	uni10::UniTensor U(bonds, "U");
	U.orthoRand();
	uni10::UniTensor UT = U;
	UT.transpose();
	
	bonds.clear();
	vector<uni10::Qnum> qnums1;
	qnums1.push_back(q30); qnums1.push_back(q11); qnums1.push_back(q11); qnums1.push_back(q11);
	qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_31);
	uni10::Bond bdr1(uni10::BD_IN, qnums1);
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	uni10::UniTensor W1(bonds, "W1");
	W1.orthoRand();
	uni10::UniTensor W1T = W1;
	W1T.transpose();

	uni10::UniTensor W2(bonds, "W2");
	W2.orthoRand();
	uni10::UniTensor W2T = W2;
	W2T.transpose();
	uni10::UniTensor W2F = W2;
	W2F.orthoRand();
	
	uni10::UniTensor H1;
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
	
	// Network putTensor()
	uni10::Network net2("AscendL");
	net2.putTensor(0, &W1);
	net2.putTensor(1, &W2);
	net2.putTensor(2, &U);
	net2.putTensor(3, &H0);
	net2.putTensor(4, &UT);
	net2.putTensor(5, &W1T);
	net2.putTensor(6, &W2T);
	H1 = net2.launch();
	net2.putTensor(0, &W1);
	net2.putTensor(0, &W1);
	net2.putTensor(0, &W1);
	H1 = net2.launch();
	cout<<net2;
	//  END putTensor()
	// Network Construct()
	/*
	vector<uni10::UniTensor*> tens;
	tens.push_back(&W1);
	tens.push_back(&W2);
	tens.push_back(&U);
	tens.push_back(&H0);
	tens.push_back(&UT);
	tens.push_back(&W1T);
	tens.push_back(&W2T);
	uni10::Network net3("AscendL", tens);
	cout<< net3;
	H1 = net3.launch();
	*/
	// END Network Construct()
	int label_out[] = {-1, -2, -3, -4};
	H1.permute(label_out, 2);
	cout<<H1;
	//printRawElem(H1);
	//cout<<W1T;
	uni10::UniTensor H1_tmp("tenH1");
	H0.check();
	assert(H1.elemCmp(H1_tmp));
}


#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"

int main(){
	Qnum_t q_20(-2, 0);
	Qnum_t q_21(-2, 1);
	Qnum_t q_10(-1, 0, 0);
	Qnum_t q_11(-1, 0, 1);
	Qnum_t q00(0, 0, 0);
	Qnum_t q01(0, 0, 1);
	Qnum_t q10(1, 0, 0);
	Qnum_t q11(1, 0, 1);
	Qnum_t q20(2, 0);
	Qnum_t q21(2, 1);
	vector<Bond_t> bonds;
	vector<Qnum_t> qnums;
	qnums.push_back(q_10);
	qnums.push_back(q_11);
	qnums.push_back(q00);
	qnums.push_back(q00);
	qnums.push_back(q01);
	qnums.push_back(q01);
	qnums.push_back(q10);
	qnums.push_back(q11);
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	
	/*
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4
					  };
	*/

	SyTensor_t H(bonds);
	//H.addRawElem(H_elem);
	H.orthoRand();
	int labels[] = {-1, -2, -3, -4};
	int rsp_labels[] = {-4, -1, -2, -3};
	int ord[] = {3, 0, 1, 2};
	vector<_Swap> swaps = _recSwap(ord, 4);
	SyTensor_t H1 = H;
	H.addLabel(labels);
	H1.addLabel(labels);
	H.reshape(rsp_labels, 1);
	H1.addGate(swaps);
	H1.reshape(rsp_labels, 1, 0);
	H.setName("H");
	H1.setName("H1");
	//cout << H;
	cout << H1;
	H.save("H_out");
	H1.save("H1_out");
}

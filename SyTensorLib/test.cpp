#include <iostream>
#include <assert.h>
#include <vector>
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
	H0.save("dada");
	cout<<H0;
	SyTensor_t H1("H1");
	cout<<H1;
}


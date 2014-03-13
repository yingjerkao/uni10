#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#define FERMION 1


int main(){
	Qnum_t q10(1, 0);
	Qnum_t q_10(-1, 0);
	vector<Bond_t> bonds;
	vector<Qnum_t> qnums;
	qnums.push_back(q10);qnums.push_back(q_10);
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	int labels[] = {1, 2, 3, 4, 5, 6};
	int rsp_labels[] = {1, 4, 3, 2, 5, 6};
	int labels1[] = {1, 2, 3, -4, -5, -6};
	int rsp_labels1[] = {1, -4, 3, 2, -5, -6};
	SyTensor_t H1(bonds);
	H1.addLabel(labels);
	H1.orthoRand();
	SyTensor_t H2(bonds);
	H2.addLabel(labels1);
	H2.orthoRand();
	cout<<H1;
	
	cout<<H1.trace(H2)<<endl;

	cout<<(H1 * H2).trace()<<endl;
	H1.reshape(rsp_labels, 3);
	H2.reshape(rsp_labels1, 3);
	cout<<(H1 * H2).trace()<<endl;
	//cout<<"Trace Hf1: "<<Hf1.trace(Hf)<<endl;
	
}


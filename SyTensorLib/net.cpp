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
	int label_H[] = {1, 5, 3, 6};
	SyTensor_t H(bonds, label_H, "H");
	int label_U[] = {7, 8, 2, 5};
	SyTensor_t U(bonds, label_U, "U");
	H.addRawElem(H_elem);
	U.orthoRand();
	SyTensor_t UT = U;
	UT.transpose();
	//cout << UT;
	int label_UT[] = {6, 9, 7, 10};
	UT.addLabel(label_UT);
	//cout<<H;
	//cout<<U;

	Node_t ndH(&H);
	Node_t ndU(&U);
	cout << ndH.metric(&ndU) <<endl;

	//SyTensor_t tmp = H * U;
	//cout<<tmp;
	//H.check();

	//cout<<ndH;
	//cout<<ndU;
	//cout << ndH.contract(&ndU)<<endl;

	Network_t net;
	net.add(&H);
	net.add(&U);
	net.add(&UT);
	net.optimize();
	cout<<net;
	
	/*
	vector<SyTensor_t*> tens;
	tens.push_back(&H);
	tens.push_back(&U);
	tens.push_back(&UT);
	Network_t net1(tens);
	*/
}

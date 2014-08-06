#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#include "MERA_Operator.cpp"

int main(){
	Qnum_t q_30(-3, 0);
	Qnum_t q_20(-2, 0);
	Qnum_t q_10(-1, 0);
	Qnum_t q00(0, 0);
	Qnum_t q10(1, 0);
	Qnum_t q20(2, 0);
	Qnum_t q30(3, 0);
	vector<Bond_t> bonds;

	vector<Qnum_t> qnum0s;
	qnum0s.push_back(q_10);
	qnum0s.push_back(q00);
	qnum0s.push_back(q10);
	Bond_t bdr(BD_ROW, qnum0s);
	Bond_t bdc(BD_COL, qnum0s);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t H(bonds, "Ob");
	SyTensor_t U(bonds, "U");

	vector<Qnum_t> qnums;
	qnums.push_back(q_30);
	for(int i= 0; i < 3; i++)
		qnums.push_back(q20);
	for(int i= 0; i < 6; i++)
		qnums.push_back(q10);
	for(int i= 0; i < 7; i++)
		qnums.push_back(q00);
	for(int i= 0; i < 6; i++)
		qnums.push_back(q_10);
	for(int i= 0; i < 3; i++)
		qnums.push_back(q_20);
	qnums.push_back(q30);
	Bond_t bdr1(BD_ROW, qnums);
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t W1(bonds, "W1");
	/*
	cout<<H;
	cout<<U;
	cout<<W1;
	*/


	H.orthoRand();
	U.orthoRand();
	W1.orthoRand();
	
	//printRawElem(H, "H_elemS");
	//printRawElem(U, "U_elemS");
	//printRawElem(W1, "W1_elemS");

	Network_t asdL("../Diagrams/AscendL");
	Network_t asdC("../Diagrams/AscendC");
	Network_t asdR("../Diagrams/AscendR");

	time_t start, end;
	start = clock();		
	SyTensor_t Tout = Ascend(H, W1, U, asdL, asdC, asdR);
	end = clock();
	printf("time = %f(sec)\n", float(end - start) / CLOCKS_PER_SEC);
	printRawElem(Tout);
}

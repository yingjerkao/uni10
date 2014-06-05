#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#include "MERA_Operator.cpp"

int main(){
	Qnum_t q_20(-2, 0);
	Qnum_t q_10(-1, 0);
	Qnum_t q00(0, 0, 0);
	Qnum_t q01(0, 0, 1);
	Qnum_t q10(1, 0);
	Qnum_t q20(2, 0);
	vector<Bond_t> bonds;
	vector<Qnum_t> qnums;
	qnums.push_back(q_20);
	for(int i = 0; i < 5; i++)
		qnums.push_back(q_10);
	for(int i = 0; i < 5; i++)
		qnums.push_back(q00);
	for(int i = 0; i < 5; i++)
		qnums.push_back(q01);
	for(int i = 0; i < 5; i++)
		qnums.push_back(q10);
	qnums.push_back(q20);
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t H(bonds, "Rho");
	SyTensor_t Rho(bonds, "Rho");
	SyTensor_t U(bonds, "U");
	bonds.clear();
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t W1(bonds, "W1");

	Rho.orthoRand();
	H.orthoRand();
	U.orthoRand();
	W1.orthoRand();
	/*	
	printRawElem(H, "H_elemS");
	printRawElem(U, "U_elemS");
	printRawElem(W1, "W1_elemS");
	printRawElem(Rho, "Rho_elemS");
	*/

	Network_t desL("../Diagrams/DescendL");
	Network_t desC("../Diagrams/DescendC");
	Network_t desR("../Diagrams/DescendR");

	time_t start, end;
	start = clock();		
	SyTensor_t Tout = Descend(Rho, W1, U, desL, desC, desR);
	end = clock();
	printf("time = %f(sec)\n", float(end - start) / CLOCKS_PER_SEC);
	
	cout<<desL;
	cout<<desC;
	cout<<desR;
	/*
	printRawElem(H);
	printRawElem(U);
	printRawElem(W1);
	*/
	printRawElem(Tout);

	/*
	cout << H;
	cout << U;
	cout << W1;
	cout << Tout;
	*/
}

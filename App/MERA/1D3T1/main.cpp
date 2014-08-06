#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#include "MERA_Operator.cpp"

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
	SyTensor_t H0(bonds, "Ob");
	H0.addRawElem(H_elem);

	SyTensor_t U(bonds, "U");
	U.orthoRand();
	
	bonds.clear();
	vector<Qnum_t> qnums1;
	qnums1.push_back(q30);qnums1.push_back(q10);qnums1.push_back(q10);qnums1.push_back(q10);
	qnums1.push_back(q_10);qnums1.push_back(q_10);qnums1.push_back(q_10);qnums1.push_back(q_30);
	Bond_t bdr1(BD_ROW, qnums1);
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t W(bonds);
	W.orthoRand();

	Network_t asdL("Diagrams/AscendL");
	Network_t asdC("Diagrams/AscendC");
	Network_t asdR("Diagrams/AscendR");
	Network_t desL("Diagrams/DescendL");
	Network_t desC("Diagrams/DescendC");
	Network_t desR("Diagrams/DescendR");
	SyTensor_t H = Ascend(H0, W, U, asdL, asdC, asdR);
	H.save("H1");
	cout<<asdL;
	/*
	H.setName("Rho");
	SyTensor_t Rho = Descend(H, W, U, desL, desC, desR);
	cout<<desC;
	cout<< Rho;
	*/

	H.check();
}

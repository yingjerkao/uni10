#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#include "MERA_Operator.cpp"

int main(){
	Qnum_t q00(0, 0);
	vector<Bond_t> bonds;
	vector<Qnum_t> qnums;
	qnums.push_back(q00);
	qnums.push_back(q00);
	qnums.push_back(q00);
	qnums.push_back(q00);
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t H(bonds, "Ob");
	SyTensor_t U(bonds, "U");
	bonds.clear();
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t W1(bonds, "W1");
	
	double *H_elem = (double*)malloc(sizeof(double) * H.getElemNum());
	double *U_elem = (double*)malloc(sizeof(double) * U.getElemNum());
	double *W1_elem = (double*)malloc(sizeof(double) * W1.getElemNum());
	FILE *fp = fopen("H_elem", "r");
	fread(H_elem, sizeof(double), H.getElemNum(), fp);
	fclose(fp);
	fp = fopen("U_elem", "r");
	fread(U_elem, sizeof(double), U.getElemNum(), fp);
	fclose(fp);
	fp = fopen("W1_elem", "r");
	fread(W1_elem, sizeof(double), W1.getElemNum(), fp);
	fclose(fp);	
	
	H.addRawElem(H_elem);	
	U.addRawElem(U_elem);	
	W1.addRawElem(W1_elem);	


	Network_t asdL("../Diagrams/AscendL");
	Network_t asdC("../Diagrams/AscendC");
	Network_t asdR("../Diagrams/AscendR");



	SyTensor_t Tout = Ascend(H, W1, U, asdL, asdC, asdR);
	cout << H;
	cout << U;
	cout << W1;
	cout << Tout;

}

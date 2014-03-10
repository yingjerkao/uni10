#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "SyTensor.h"

int main(){
	Qnum_t q10(0, 0);
	Qnum_t q00(0, 0);
	Qnum_t q_10(0, 0);
	vector<Bond_t> bonds;		
	vector<Qnum_t> qnums;
	qnums.push_back(q10);qnums.push_back(q00);qnums.push_back(q_10);;
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdc);

	SyTensor_t Sz(bonds, "Sz");
	int label_tmp[] = {-1, -3};
	vector<int> labels(label_tmp, label_tmp + sizeof(label_tmp) / sizeof(int));
	Sz.addLabel(labels);
	double Sz_elem[] = {1, 0, 0,\
					    0, 0, 0,\
					    0, 0,-1\
					   };
	Sz.addRawElem(Sz_elem);
	SyTensor_t Sz2 = Sz;
	int label2_tmp[] = {-2, -4};
	vector<int> labels2(label2_tmp, label2_tmp + sizeof(label2_tmp) / sizeof(int));
	Sz2.addLabel(labels2);
	//rSyT.reshape(labels2, 2);
	SyTensor_t SzSz = Sz * Sz2;

	int label3_tmp[] = {-1, -2, -3, -4};
	vector<int> labels3(label3_tmp, label3_tmp + sizeof(label3_tmp) / sizeof(int));
	SzSz.reshape(labels3, 2);
	cout<<SzSz;

	SyTensor_t Sp(bonds, "Sp");
	Sp.addLabel(labels);
	double Sp_elem[] = {0, 1, 0,\
					    0, 0, 1,\
					    0, 0, 0\
					   };
	Sp.addRawElem(Sp_elem);
	Sp *= sqrt(2);
	SyTensor_t Sp2 = Sp;
	Sp2.addLabel(labels2);
	SyTensor_t SpSp = Sp * Sp2;
	SpSp.reshape(labels3, 2);
	cout<<SpSp;

	SyTensor_t Sm(bonds, "Sm");
	Sm.addLabel(labels);
	double Sm_elem[] = {0, 0, 0,\
					    1, 0, 0,\
					    0, 1, 0\
					   };
	Sm.addRawElem(Sm_elem);
	Sm *= sqrt(2);
	SyTensor_t Sm2 = Sm;
	Sm2.addLabel(labels2);
	SyTensor_t SmSm = Sm * Sm2;
	SmSm *= sqrt(2.0);	
	SmSm.reshape(labels3, 2);
	cout<<SmSm;

	SyTensor_t SS = 0.5 * (Sp * Sm2 + Sm * Sp2) + Sz * Sz2;
	SS.reshape(labels3, 2);
	cout<<SS;
	printRawElem(SS);

	Sz.check();
}


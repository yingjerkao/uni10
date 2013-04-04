#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "SyTensor.h"

int main(){
	Qnum_t q10(1, 0);
	Qnum_t q00(0, 0);
	Qnum_t q_10(-1, 0);
	vector<Bond_t> bonds;		
	vector<Qnum_t> qnums;
	qnums.push_back(q10);qnums.push_back(q00);qnums.push_back(q_10);;
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdc);

	SyTensor_t Sz(bonds, "Sz");
	int label_tmp[] = {-1, -2, -3, -4};
	vector<int> labels(label_tmp, label_tmp + sizeof(label_tmp) / sizeof(int));
	//SyT.addLabel(labels);
	double Sz_elem[] = {1, 0, 0,\
					 0, 0, 0,\
					 0, 0,-1\
					};
	Sz.addRawElem(Sz_elem);
	int label2_tmp[] = {-2, -3, -1, -4};
	vector<int> labels2(label2_tmp, label2_tmp + sizeof(label2_tmp) / sizeof(int));
	//rSyT.reshape(labels2, 2);
	cout << Sz;		            
	Sz.check();
}


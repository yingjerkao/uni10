#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"

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
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);

	SyTensor_t SyT(bonds);
	int label_tmp[] = {-1, -2, -3, -4};
	vector<int> labels(label_tmp, label_tmp + sizeof(label_tmp) / sizeof(int));
	SyT.addLabel(labels);
	double elem[] = {1, 0, 0, 0, 0, 0, 0, 0, 0,\
					 0, 0, 0, 1, 0, 0, 0, 0, 0,\
					 0, 0,-1, 0, 1, 0, 0, 0, 0,\
					 0, 1, 0, 0, 0, 0, 0, 0, 0,\
					 0, 0, 1, 0, 0, 0, 1, 0, 0,\
					 0, 0, 0, 0, 0, 0, 0, 1, 0,\
					 0, 0, 0, 0, 1, 0,-1, 0, 0,\
					 0, 0, 0, 0, 0, 1, 0, 0, 0,\
					 0, 0, 0, 0, 0, 0, 0, 0, 1\
					};

	SyT.addRawElem(elem);
	cout << SyT;
}

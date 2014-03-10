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
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);

	SyTensor_t SyT(bonds, "HelloD");
	cout<<"LALA"<<endl;
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
	/*
	double elem[] = {1, 0, 0, 0, 0, 0, 0, 0, 0,\
					 0, 0, 0, 2, 0, 0, 0, 0, 0,\
					 0, 0, 3, 0, 4, 0, 0, 0, 0,\
					 0, 5, 0, 0, 0, 0, 0, 0, 0,\
					 0, 0, 6, 0, 0, 0, 7, 0, 0,\
					 0, 0, 0, 0, 0, 0, 0, 8, 0,\
					 0, 0, 0, 0, 9, 0,10, 0, 0,\
					 0, 0, 0, 0, 0,11, 0, 0, 0,\
					 0, 0, 0, 0, 0, 0, 0, 0,12\
					};
	*/
	SyT.addRawElem(elem);
	SyTensor_t rSyT = SyT;
	//rSyT.reshape(labels2, 2);
	cout << SyT;		            
	//cout<< rSyT;
	//cout<<SyT;
	printRawElem(SyT);
	SyT.check();
}

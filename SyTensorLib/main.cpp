#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "SyTensor.h"

int main(){
	Qnum_t q11(1, 1);
	/*
	Qnum_t q10(1, 0);
	Qnum_t q00(0, 0);
	Qnum_t q01(0, 1);
	Qnum_t q_10(-1, 0);
	Qnum_t q_11(-1, 1);
	*/
	vector<Qnum_t> qnums;
	vector<Qnum_t> qnums2;
	qnums.push_back(q11);
	qnums2.push_back(-q11);
	vector<Bond_t> bonds;
		
	Bond_t bd(BD_ROW, qnums);
	printf("main line 25\n");
	Bond_t bd2(BD_ROW, qnums2);
	printf("main line 27\n");
	bonds.push_back(bd);
	printf("main line 29\n");
	bonds.push_back(bd2);
	printf("main line 31\n");
	SyTensor_t SyT(bonds, "HelloD");
	cout<< SyT;
}

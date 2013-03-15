#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "Qnum.h"
#include "Bond.h"
#include "Block.h"

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
	qnums.push_back(q11);
	
	Bond_t bd(BD_ROW, qnums);
	Bond_t bd2(BD_ROW, qnums);
	cout<< bd;
	Block_t blk(12, 15, 6);
	cout<< blk << endl;
}

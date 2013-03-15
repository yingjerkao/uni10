#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
using namespace std;
#include "Qnum.h"
#include "Bond.h"

Bond_t::Bond_t(bondType _type, vector<Qnum_t>& qnums) : type(_type), dim(0){
	map<Qnum_t, bool> mark;
	int cnt = 0;
	for(int i = 0; i < qnums.size(); i++){
		if(mark.find(qnums[i]) == mark.end()){
			cout<<"line32\n";
			mark[ qnums[i] ] = true;
			cout<<"line34\n";
			Qnums.push_back(qnums[i]);
			cout<<"line36\n";
			Qdegs.push_back(1);
			offsets.push_back(dim);
			cnt++;
		}
		else{
			assert(qnums[i - 1] == qnums[i]);
			Qdegs[cnt - 1]++;
		}
		dim++;
	}
}
Bond_t::~Bond_t(){
	cout<<"Destructing Bond...\n";
}

ostream& operator<< (ostream& os, const Bond_t& b){
	if(b.type == BD_ROW)
		os<<"ROW: ";
	else
		os<<"COL: ";
	for(int i = 0; i < b.Qnums.size(); i++)
		os << b.Qnums[i] << "|" << b.Qdegs[i] << ", ";
	os<<"Dim = "<< b.dim << endl;
	return os;
}

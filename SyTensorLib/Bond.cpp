#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
using namespace std;
#include "Qnum.h"
#include "Bond.h"

Bond_t::Bond_t(bondType _type, vector<Qnum_t>& qnums) : type(_type){
	setting(qnums);
}

void Bond_t::assign(bondType _type, vector<Qnum_t>& qnums){
	type = _type;
	Qnums.clear();
	Qdegs.clear();
	offsets.clear();
	setting(qnums);
}

void Bond_t::setting(vector<Qnum_t>& qnums){
	map<Qnum_t, bool> mark;
	int cnt = 0;
	dim = 0;
	//cout<<"Constructing Bond "<< this << endl;
	for(int i = 0; i < qnums.size(); i++){
		if(mark.find(qnums[i]) == mark.end()){
			mark[ qnums[i] ] = true;
			Qnums.push_back(qnums[i]);
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
	//cout<<"Destructing Bond "<< this << endl;
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

bool operator== (const Bond_t& b1, const Bond_t& b2){
	return (b1.type == b2.type) && (b1.Qnums == b2.Qnums) && (b1.Qdegs == b2.Qdegs);
}

#include "../QnumF.h"
using namespace uni10::datatype;
#include "../Bond.h"
Bond_t::Bond_t(bondType _type, std::vector<Qnum_t>& qnums) : type(_type){
	setting(qnums);
}

void Bond_t::assign(bondType _type, std::vector<Qnum_t>& qnums){
	type = _type;
	Qnums.clear();
	Qdegs.clear();
	offsets.clear();
	setting(qnums);
}

void Bond_t::setting(std::vector<Qnum_t>& qnums){
	assert(qnums.size() > 0);
	std::map<Qnum_t, bool> mark;
	int cnt = 0;
	dim = 0;
	//cout<<"Constructing Bond "<< this << endl;
	for(int i = 0; i < qnums.size(); i++){
		if(i == 0 || !(qnums[i] == qnums[i - 1])){
			Qnums.push_back(qnums[i]);
			Qdegs.push_back(1);
			offsets.push_back(dim);
			cnt++;
		}
		else
			Qdegs[cnt - 1]++;
		dim++;
		/*
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
		*/
	}
}

Bond_t::~Bond_t(){
	//cout<<"Destructing Bond "<< this << endl;
}

std::ostream& operator<< (std::ostream& os, const Bond_t& b){
	if(b.type == BD_ROW)
		os<<"ROW: ";
	else
		os<<"COL: ";
	for(int i = 0; i < b.Qnums.size(); i++)
		os << b.Qnums[i] << "|" << b.offsets[i] << "-"<<b.Qdegs[i]<<", ";
	os<<"Dim = "<< b.dim << std::endl;
	return os;
}

bool operator== (const Bond_t& b1, const Bond_t& b2){
	return (b1.type == b2.type) && (b1.Qnums == b2.Qnums) && (b1.Qdegs == b2.Qdegs);
}
void Bond_t::change(bondType tp){
	if(type != tp){
		for(int q = 0; q < Qnums.size(); q++)
			Qnums[q] = -Qnums[q];
		type = tp;
	}
}

void Bond_t::combine(Bond_t bd){
	bd.change(type);
	std::vector<Qnum_t> qnums;
	std::vector<int> qdegs;
	offsets.clear();
	dim = 0;
	Qnum_t qnum;
	int qdim;
	int cnt = 0;
	for(int q = 0; q < Qnums.size(); q++)
		//for(int d = 0; d < Qdegs[q]; d++){
			for(int qq = 0; qq < bd.Qnums.size(); qq++){
				qnum = Qnums[q] * bd.Qnums[qq];
				qdim = Qdegs[q] * bd.Qdegs[qq];
				if(qnums.size() == 0 || !(qnum == qnums[cnt - 1])){
					qnums.push_back(qnum);
					qdegs.push_back(qdim);
					offsets.push_back(dim);
					cnt++;
				}
				else{
					qdegs[cnt - 1] += qdim;
				}
				dim += qdim;
			}
		//}
	Qnums = qnums;
	Qdegs = qdegs;
}
Bond_t combine(bondType tp, const std::vector<Bond_t>& bds){
	assert(bds.size() > 1);
	int bd_num = bds.size();
	Bond_t outBond1 = bds[bd_num - 1];
	Bond_t outBond2 = bds[bd_num - 2];
	int b = 0;
	outBond2.change(tp);
	outBond2.combine(outBond1);
	for(b = 0; b < bd_num - 2; b++){
		if(b % 2 == 0){
			outBond1 = bds[bd_num - 3 - b];
			outBond1.change(tp);
			outBond1.combine(outBond2);
		}
		else{
			outBond2 = bds[bd_num - 3 - b];
			outBond2.change(tp);
			outBond2.combine(outBond1);
		}
	}
	if(b % 2 == 0)
		return outBond2;
	else
		return outBond1;	
}
Bond_t combine(const std::vector<Bond_t>& bds){
	return combine(bds[0].type, bds);
}

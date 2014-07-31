#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){
	Qnum q0;
	Bond bdi(BD_IN, 5);
	Bond bdo(BD_OUT, 5);
	vector<Bond> bonds(2, bdi);
	bonds.push_back(bdo);
	bonds.push_back(bdo);
	//UniTensor T(bonds);
	//T.orthoRand();
	int per_label[] = {1, 3, 2, 0};
	UniTensor T("Telem");
	cout<<T;
	//T.save("Telem");
	T.permute(per_label, 1);
	cout<<T;

	uni10::Qnum q1(1);
	uni10::Qnum q_1(-1);
	uni10::Qnum q2(2);
	uni10::Qnum q_2(-2);
	std::vector<uni10::Qnum> sy_qnums;
	sy_qnums.push_back(q2);
	sy_qnums.push_back(q1);
	sy_qnums.push_back(q0);
	sy_qnums.push_back(q0);
	sy_qnums.push_back(q_1);
	sy_qnums.push_back(q_2);
	uni10::Bond bd_in(uni10::BD_IN, sy_qnums);
	uni10::Bond bd_out(uni10::BD_OUT, sy_qnums);
	std::vector<uni10::Bond> sy_bonds;
	sy_bonds.push_back(bd_in);
	sy_bonds.push_back(bd_in);
	sy_bonds.push_back(bd_out);
	sy_bonds.push_back(bd_out);
	UniTensor ST(sy_bonds);
	ST.orthoRand();
	//cout<<ST;
}

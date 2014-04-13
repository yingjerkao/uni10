#include <iostream>
#include <uni10.hpp>

int main(){
	uni10::Qnum q1(1);	
	uni10::Qnum q0(0);	
	uni10::Qnum q_1(-1);	
	std::vector<uni10::Qnum> qnums;
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);

	uni10::Bond bd1(uni10::BD_IN, qnums);
	uni10::Bond bd2(uni10::BD_IN, qnums);

	return 0;
}



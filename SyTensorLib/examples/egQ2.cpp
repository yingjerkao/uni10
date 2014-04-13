#include <iostream>
#include <uni10.hpp>

int main(){
	// U1 = 1, parity even.
	uni10::Qnum q10(1, uni10::PRT_EVEN);	
	// U1 = -1, parity odd.
	uni10::Qnum q_11(-1, uni10::PRT_ODD);	

	std::cout<<"q10: "<< q10 << std::endl;
	std::cout<<"q_11: "<< q_11 << std::endl;

	std::cout<<"----- Operations -----\n";
	std::cout<<"-q_11 = " << -q_11 << std::endl;
	std::cout<<"q10 * q_11 = " << q10 * q_11 <<std::endl;
	std::cout<<"q10 * (-q_11) = " << q10 * (-q_11) <<std::endl;

	return 0;
}


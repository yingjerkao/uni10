#include <iostream>
#include <uni10.hpp>

int main(){
	// U1 = 1, parity even.
	uni10::Qnum q10(1, uni10::PRT_EVEN);	
	// U1 = -1, parity odd.
	uni10::Qnum q_11(-1, uni10::PRT_ODD);	

	std::cout<<"q10: "<<q10<<std::endl;
	std::cout<<"q_11: "<<q_11<<std::endl;
	std::cout<<"q_11: U1 = "<<q_11.U1()<<", parity = "<<q_11.prt()<<std::endl; 
	q_11.assign(-2, uni10::PRT_EVEN);
	std::cout<<"q_11(after assign): "<<q_11<<std::endl;
	// check the for fermionic
	std::cout<<"isFermioinc: "<<uni10::Qnum::isFermionic()<<std::endl ;

	// Fermionic system
	std::cout<<"----- Fermionic -----\n";
	// fermionic parity even, U1 = 1, parity even.
	uni10::Qnum f0_q10(uni10::PRTF_EVEN, 1, uni10::PRT_EVEN);	
	// fermionic parity odd, U1 = 1, parity even.
	uni10::Qnum f1_q10(uni10::PRTF_ODD, 1, uni10::PRT_EVEN);	

	std::cout<<"f0_q10: "<<f0_q10<<std::endl;
	std::cout<<"f1_q10: "<<f1_q10<<std::endl;
	std::cout<<"f1_q10: fermionic parity = " <<f1_q10.prtF()<<std::endl;
	std::cout<<"isFermioinc: "<<uni10::Qnum::isFermionic()<<std::endl;

	return 0;
}


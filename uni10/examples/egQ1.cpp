/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egQ1.cpp
* 
*  @license
*  Copyright (C) 2013-2014 
*  This file is part of Uni10
*  
*  Uni10 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/
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
/* OUTPUT

q10: (U1 = 1, P = 0, 0)
q_11: (U1 = -1, P = 1, 0)
q_11: U1 = -1, parity = 1
q_11(after assign): (U1 = -2, P = 0, 0)
isFermioinc: 0
----- Fermionic -----
f0_q10: (U1 = 1, P = 0, 0)
f1_q10: (U1 = 1, P = 0, 1)
f1_q10: fermionic parity = 1
isFermioinc: 1

*/

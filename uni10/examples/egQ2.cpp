/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egQ2.cpp
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

	std::cout<<"q10: "<< q10 << std::endl;
	std::cout<<"q_11: "<< q_11 << std::endl;

	std::cout<<"----- Operations -----\n";
	std::cout<<"-q_11 = " << -q_11 << std::endl;
	std::cout<<"q10 * q_11 = " << q10 * q_11 <<std::endl;
	std::cout<<"q10 * (-q_11) = " << q10 * (-q_11) <<std::endl;

	return 0;
}
/* Output

q10: (U1 = 1, P = 0, 0)
q_11: (U1 = -1, P = 1, 0)
----- Operations -----
-q_11 = (U1 = 1, P = 1, 0)
q10 * q_11 = (U1 = 0, P = 1, 0)
q10 * q_11 = (U1 = 2, P = 1, 0)

*/

/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egM1.cpp
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
	// Spin 1/2 Heisenberg hamiltonian
	double elem[] = {1.0/4,      0,      0,     0,
						 0, -1.0/4,  1.0/2,     0,
						 0,  1.0/2, -1.0/4,     0,
						 0,      0,      0, 1.0/4};
	uni10::Matrix H(4, 4, elem);
	std::cout<<H;
	// Diagonlize H
	std::vector<uni10::Matrix> results = H.eigh();
	std::cout<<"The eigen values: \n\n"<<results[0];
	std::cout<<"The eigen vectors: \n\n"<<results[1];

	// Access element in a diagonal matrix
	uni10::Matrix D = results[0];
	std::cout<<"D.at(1, 1) = "<<D.at(1, 1)<<std::endl;;
	std::cout<<"D[2] = " << D[2]<<std::endl;
	// Assign element
	std::cout<<"\nAssign D.at(3, 3) = 7.0\n\n";
	D.at(3, 3) = 7.0;
	std::cout<<D;

	// Access element
	std::cout<<"H.at(1, 2) = "<<H.at(1, 2)<<std::endl;;
	std::cout<<"H[5] = " << H[5]<<std::endl;

	// Make a pure density matrix from ground state
	uni10::Matrix U = results[1];
	// Generate ground state by taking the first H.rol() elements from U.
	uni10::Matrix GS(1, H.col(), U.getElem());
	// Transposed GS
	uni10::Matrix GST = GS;
	GST.transpose();

	std::cout<<"\nThe ground state: \n\n";
	std::cout<< GS;
	std::cout<< GST;

	// Compose a pure density matrix from ground state
	uni10::Matrix Rho = GST * GS;
	std::cout<<"\nPure density matrix of ground state: \n\n";
	std::cout<< Rho;

	// Measure ground state energy
	std::cout<<"\nThe ground state energy: " << (Rho * H).trace() << std::endl;
}

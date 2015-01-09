/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egM2.cpp
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
	uni10::Matrix M(4, 5);
	M.randomize();
	std::cout<<M;
	// carry out SVD
	std::vector<uni10::Matrix> rets = M.svd();
	
	// write matrice out to file
	rets[0].save("mat_U");
	rets[1].save("mat_Sigma");
	rets[2].save("mat_VT");

	uni10::Matrix U(rets[0].row(), rets[0].col(), rets[0].isDiag());
	uni10::Matrix S(rets[1].row(), rets[1].col(), rets[1].isDiag());
	uni10::Matrix VT(rets[2].row(), rets[2].col(), rets[2].isDiag());

	// read in the matrice we just write out
	U.load("mat_U");
	S.load("mat_Sigma");
	VT.load("mat_VT");
	std::cout<< S;
	std::cout<< U * S * VT;
}
/* Output:
4 x 5 = 20

  0.840  0.394  0.783  0.798  0.912

  0.198  0.335  0.768  0.278  0.554

  0.477  0.629  0.365  0.513  0.952

  0.916  0.636  0.717  0.142  0.607

4 x 4 = 4, Diagonal

  2.736  0.000  0.000  0.000

  0.000  0.555  0.000  0.000

  0.000  0.000  0.449  0.000

  0.000  0.000  0.000  0.382

4 x 5 = 20

  0.840  0.394  0.783  0.798  0.912

  0.198  0.335  0.768  0.278  0.554

  0.477  0.629  0.365  0.513  0.952

  0.916  0.636  0.717  0.142  0.607
*/

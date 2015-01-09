/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egU3.cpp
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
	// Construct spin 1 Heisenberg model by reading in the tensor which is written out in example egU1
	uni10::UniTensor H_U1("egU1_H_U1");

	// Randomly create an isometry tensor
	uni10::Qnum q0(0);
	uni10::Qnum q1(1);
	uni10::Qnum q_1(-1);
	uni10::Qnum q2(2);
	uni10::Qnum q_2(-2);
	std::vector<uni10::Qnum> in_qnums;
	in_qnums.push_back(q2);
	in_qnums.push_back(q1);
	in_qnums.push_back(q0);
	in_qnums.push_back(q0);
	in_qnums.push_back(q_1);
	in_qnums.push_back(q_2);
	std::vector<uni10::Qnum> out_qnums;
	out_qnums.push_back(q1);
	out_qnums.push_back(q0);
	out_qnums.push_back(q_1);
	uni10::Bond bd_in(uni10::BD_IN, in_qnums);
	uni10::Bond bd_out(uni10::BD_OUT, out_qnums);
	std::vector<uni10::Bond> bonds;
	bonds.push_back(bd_in);
	bonds.push_back(bd_out);
	bonds.push_back(bd_out);
	// Create isometry tensor W and transposed WT
	uni10::UniTensor W(bonds, "W");
	W.orthoRand();
	uni10::UniTensor WT = W;
	WT.transpose();

	// Operate W and WT on H_U1, see the contraction labels in the documentation.
	int label_H[] = {1, 2, 3, 4};
	int label_W[] = {-1, 1, 2};
	int label_WT[] = {3, 4, -2};
	H_U1.setLabel(label_H);
	W.setLabel(label_W);
	WT.setLabel(label_WT);
	//std::cout<<W;
	std::cout<<W * H_U1 * WT;

	// Write the tensors W and WT out to file
	W.save("egU3_W");
	WT.save("egU3_WT");

	// Check the memory usage.
	uni10::UniTensor::profile();

	return 0;
}



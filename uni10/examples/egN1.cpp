/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egN1.cpp
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
	// Read in the tensor H_U1 which is written out in example egU1 and W, WT in example egU3
	uni10::UniTensor H_U1("egU1_H_U1");
	uni10::UniTensor W("egU3_W");
	uni10::UniTensor WT("egU3_WT");

	// Create network by reading in network file "egN1_network"
	uni10::Network net("egN1_network");
	// Put tensors to the Network net
	net.putTensor("H", H_U1);
	net.putTensor("W", W);
	net.putTensor("WT", WT);

	// Perform contractions inside the tensor network
	std::cout<<net.launch();
	// Print out the network
	std::cout<<net;
  // Print out the memory usage
  net.profile();

	return 0;
}



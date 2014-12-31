/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egB2.cpp
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



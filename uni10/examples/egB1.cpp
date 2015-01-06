/*
*
*  Universal Tensor Network Library (Uni10)
*  @file
*  egB1.cpp
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
	// Create an array of Qnums for the states of a bond.
	std::vector<uni10::Qnum> qnums;
	qnums.push_back(q1);
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);

	// Constrcut Bond with Qnum array
	uni10::Bond bd(uni10::BD_IN, qnums);
	// Print out a Bond
	std::cout<<"Bond bd: \n"<<bd<<std::endl;
	std::cout<<"Bond type: "<<bd.type()<<"(IN)"<<", Bond dimension: "<<bd.dim()<<std::endl<<std::endl;

	// List the degeneracy of states
	std::cout<<"Degeneracies: "<<std::endl;
	std::map<uni10::Qnum, int> degs = bd.degeneracy();
	for(std::map<uni10::Qnum,int>::const_iterator it=degs.begin(); it!=degs.end(); ++it)
		std::cout<<it->first<<": "<<it->second<<std::endl;
	std::cout<<std::endl;

	std::vector<uni10::Qnum> qlist = bd.Qlist();
	std::cout<<"Qnum list: "<<std::endl;
	for(int i = 0; i < qlist.size(); i++)
		std::cout<<qlist[i]<<", ";
	std::cout<<std::endl<<std::endl;

	// Change bond type
	bd.change(uni10::BD_OUT);
	std::cout<<"bd changes to BD_OUT:\n"<<bd<<std::endl;

	bd.change(uni10::BD_IN);

	// Combine bond
	qnums.clear();
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);
	uni10::Bond bd2(uni10::BD_IN, qnums);
	std::cout<<"Bond bd2: \n"<<bd2<<std::endl;

	// bd.combine(bd2);
	std::cout<<"bd2.combine(bd): \n"<<bd2.combine(bd)<<std::endl;

	std::cout<<"Degeneracies of bd2 after combining bd: "<<std::endl;
	degs = bd2.degeneracy();
	for(std::map<uni10::Qnum,int>::const_iterator it=degs.begin(); it!=degs.end(); ++it)
		std::cout<<it->first<<": "<<it->second<<std::endl;
	std::cout<<std::endl;

	return 0;
}

/* Output

Bond bd: 
IN : (U1 = 1, P = 0, 0)|2, (U1 = 0, P = 0, 0)|3, (U1 = -1, P = 0, 0)|1, Dim = 6

Bond type: 1(IN), Bond dimension: 6

Degeneracies: 
(U1 = -1, P = 0, 0): 1
(U1 = 0, P = 0, 0): 3
(U1 = 1, P = 0, 0): 2

Qnum list: 
(U1 = 1, P = 0, 0), (U1 = 1, P = 0, 0), (U1 = 0, P = 0, 0), (U1 = 0, P = 0, 0), (U1 = 0, P = 0, 0), (U1 = -1, P = 0, 0), 

bd changes to BD_OUT:
OUT: (U1 = -1, P = 0, 0)|2, (U1 = 0, P = 0, 0)|3, (U1 = 1, P = 0, 0)|1, Dim = 6

Bond bd2: 
IN : (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|2, (U1 = -1, P = 0, 0)|1, Dim = 4

bd2.combine(bd): 
IN : (U1 = 2, P = 0, 0)|2, (U1 = 1, P = 0, 0)|3, (U1 = 0, P = 0, 0)|1, (U1 = 1, P = 0, 0)|2, (U1 = 0, P = 0, 0)|3, (U1 = -1, P = 0, 0)|1, (U1 = 1, P = 0, 0)|2, (U1 = 0, P = 0, 0)|3, (U1 = -1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|2, (U1 = -1, P = 0, 0)|3, (U1 = -2, P = 0, 0)|1, Dim = 24

Degeneracies of bd2 after combining bd: 
(U1 = -2, P = 0, 0): 1
(U1 = -1, P = 0, 0): 5
(U1 = 0, P = 0, 0): 9
(U1 = 1, P = 0, 0): 7
(U1 = 2, P = 0, 0): 2


*/

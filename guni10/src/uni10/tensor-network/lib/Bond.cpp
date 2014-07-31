/****************************************************************************
*  @file CMakeLists.txt
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao 
*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Main specification file for CMake
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
using namespace std;
#include "QnumF.h"
#include "Bond.h"

Bond_t::Bond_t(bondType _type, vector<Qnum_t>& qnums) : type(_type){
	setting(qnums);
}

void Bond_t::assign(bondType _type, vector<Qnum_t>& qnums){
	type = _type;
	Qnums.clear();
	Qdegs.clear();
	offsets.clear();
	setting(qnums);
}

void Bond_t::setting(vector<Qnum_t>& qnums){
	assert(qnums.size() > 0);
	map<Qnum_t, bool> mark;
	int cnt = 0;
	dim = 0;
	//cout<<"Constructing Bond "<< this << endl;
	for(int i = 0; i < qnums.size(); i++){
		if(mark.find(qnums[i]) == mark.end()){
			mark[ qnums[i] ] = true;
			Qnums.push_back(qnums[i]);
			Qdegs.push_back(1);
			offsets.push_back(dim);
			cnt++;
		}
		else{
			assert(qnums[i - 1] == qnums[i]);
			Qdegs[cnt - 1]++;
		}
		dim++;
	}
}

Bond_t::~Bond_t(){
	//cout<<"Destructing Bond "<< this << endl;
}

ostream& operator<< (ostream& os, const Bond_t& b){
	if(b.type == BD_ROW)
		os<<"ROW: ";
	else
		os<<"COL: ";
	for(int i = 0; i < b.Qnums.size(); i++)
		os << b.Qnums[i] << "|" << b.Qdegs[i] << ", ";
	os<<"Dim = "<< b.dim << endl;
	return os;
}

bool operator== (const Bond_t& b1, const Bond_t& b2){
	return (b1.type == b2.type) && (b1.Qnums == b2.Qnums) && (b1.Qdegs == b2.Qdegs);
}
void Bond_t::change(bondType tp){
	if(type != tp)
		for(int q = 0; q < Qnums.size(); q++)
			Qnums[q] = -Qnums[q];
	type = tp;
}

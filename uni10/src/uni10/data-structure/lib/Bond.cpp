/****************************************************************************
*  @file Bond.cpp
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

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
*  @brief Implementation file for Matrix class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/data-structure/Bond.h>
#include <uni10/tools/uni10_tools.h>

namespace uni10{

Bond::Bond(bondType _type, size_t dim) : m_type(_type){
  try{
    Qnum q0(0);
    std::vector<Qnum> qnums(dim, q0);
    setting(qnums);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Bond::Bond(bondType, size_t):");
  }
}
Bond::Bond(bondType _type, const std::vector<Qnum>& qnums) : m_type(_type){
  try{
    setting(qnums);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Bond::Bond(bondType, std::vector<Qnum>&):");
  }
}

Bond::Bond(const Bond& _b):m_type(_b.m_type), m_dim(_b.m_dim), Qnums(_b.Qnums), Qdegs(_b.Qdegs), offsets(_b.offsets){
}
bondType Bond::type()const{
	return m_type;
}
int Bond::dim()const{
	return m_dim;
}

void Bond::assign(bondType _type, size_t dim){
  try{
    m_type = _type;
    Qnums.clear();
    Qdegs.clear();
    offsets.clear();
    Qnum q0(0);
    std::vector<Qnum> qnums(dim, q0);
    setting(qnums);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Bond::assign(bondType, size_t):");
  }
}

void Bond::assign(bondType _type, const std::vector<Qnum>& qnums){
  try{
    m_type = _type;
    Qnums.clear();
    Qdegs.clear();
    offsets.clear();
    setting(qnums);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Bond::assign(bondType, std::vector<Qnum>&):");
  }
}

void Bond::setting(const std::vector<Qnum>& qnums){
  if(!(qnums.size() > 0)){
    std::ostringstream err;
    err<<"Cannot create a bond of dimension 0.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	std::map<Qnum, bool> mark;
	int cnt = 0;
	m_dim = 0;
	for(int i = 0; i < qnums.size(); i++){
		if(i == 0 || !(qnums[i] == qnums[i - 1])){
			Qnums.push_back(qnums[i]);
			Qdegs.push_back(1);
			offsets.push_back(m_dim);
			cnt++;
		}
		else
			Qdegs[cnt - 1]++;
		m_dim++;
	}
}

Bond::~Bond(){}

std::ostream& operator<< (std::ostream& os, const Bond& b){
	if(b.m_type == BD_IN)
		os<<"IN : ";
	else
		os<<"OUT: ";
	for(int i = 0; i < b.Qnums.size(); i++)
		os << b.Qnums[i] << "|" << b.Qdegs[i]<<", ";
	os<<"Dim = "<< b.m_dim << std::endl;
	return os;
}

std::map<Qnum, int> Bond::degeneracy()const{
	std::map<Qnum, int>hst;
	for(int i = 0; i < Qnums.size(); i++){
		if(hst.find(Qnums[i]) == hst.end())
			hst[Qnums[i]] = Qdegs[i];
		else
			hst[Qnums[i]] += Qdegs[i];
	}
	return hst;
}

std::vector<Qnum> Bond::Qlist()const{
	std::vector<Qnum>list(m_dim);
	int cnt = 0;
	for(int q = 0; q < Qnums.size(); q++)
		for(int d = 0; d < Qdegs[q]; d++){
			list[cnt] = Qnums[q];
			cnt++;
		}
	return list;
}

bool operator== (const Bond& b1, const Bond& b2){
	return (b1.m_type == b2.m_type) && (b1.Qnums == b2.Qnums) && (b1.Qdegs == b2.Qdegs);
}
Bond& Bond::change(bondType tp){
	if(m_type != tp){
		for(int q = 0; q < Qnums.size(); q++)
			Qnums[q] = -Qnums[q];
		m_type = tp;
	}
  return *this;
}
Bond& Bond::dummy_change(bondType tp){
	if(m_type != tp)
		m_type = tp;
  return *this;
}
Bond& Bond::combine(Bond bd){
  try{
    bd.change(m_type);
    std::vector<Qnum> qnums;
    std::vector<int> qdegs;
    offsets.clear();
    m_dim = 0;
    Qnum qnum;
    int qdim;
    int cnt = 0;
    for(int q = 0; q < Qnums.size(); q++)
      for(int d = 0; d < Qdegs[q]; d++){
        for(int qq = 0; qq < bd.Qnums.size(); qq++){
          qnum = Qnums[q] * bd.Qnums[qq];
          qdim = bd.Qdegs[qq];
          if(qnums.size() == 0 || !(qnum == qnums[cnt - 1])){
            qnums.push_back(qnum);
            qdegs.push_back(qdim);
            offsets.push_back(m_dim);
            cnt++;
          }
          else
            qdegs[cnt - 1] += qdim;
          m_dim += qdim;
        }
      }
    Qnums = qnums;
    Qdegs = qdegs;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function combine(uni10::Bond bd):");
  }
  return *this;
}

Bond combine(bondType tp, const std::vector<Bond>& bds){
  try{
    if((bds.size() == 0)){
      std::ostringstream err;
      err<<"There should be at least one bond in the input vector to be combined.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(bds.size() == 1){
      Bond bd = bds[0];
      bd.change(tp);
      return bd;
    }
    int bd_num = bds.size();
    Bond outBond1 = bds[bd_num - 1];
    Bond outBond2 = bds[bd_num - 2];
    int b = 0;
    outBond2.change(tp);
    outBond2.combine(outBond1);
    for(b = 0; b < bd_num - 2; b++){
      if(b % 2 == 0){
        outBond1 = bds[bd_num - 3 - b];
        outBond1.change(tp);
        outBond1.combine(outBond2);
      }
      else{
        outBond2 = bds[bd_num - 3 - b];
        outBond2.change(tp);
        outBond2.combine(outBond1);
      }
    }
    if(b % 2 == 0)
      return outBond2;
    else
      return outBond1;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function combine(bondType, std::vector<uni10::Bond>&):");
    return Bond();
  }

}
Bond combine(const std::vector<Bond>& bds){
  try{
    if((bds.size() == 0)){
      std::ostringstream err;
      err<<"There should be at least one bond in the input vector to be combined.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return combine(bds[0].m_type, bds);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function combine(std::vector<uni10::Bond>&):");
    return Bond();
  }
}
};	/* namespace uni10 */

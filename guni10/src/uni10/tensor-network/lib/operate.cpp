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
#include <uni10/tensor-network/UniTensor.h>
#include <uni10/numeric/uni10_lapack.h>
//using namespace uni10::datatype;
namespace uni10{
UniTensor& UniTensor::operator+= (const UniTensor& Tb){
	assert(status & Tb.status & HAVEELEM);
	assert(bonds == Tb.bonds);
	assert(blocks == Tb.blocks);
	vectorAdd(elem, Tb.elem, m_elemNum, ongpu, Tb.ongpu);
	return *this;
}
UniTensor operator+(const UniTensor& Ta, const UniTensor& Tb){
	assert(Ta.status & Tb.status & Ta.HAVEELEM);
	assert(Ta.bonds == Tb.bonds);
	assert(Ta.blocks == Tb.blocks);
	UniTensor Tc(Ta);
	vectorAdd(Tc.elem, Tb.elem, Tc.m_elemNum, Tc.ongpu, Tb.ongpu);
	return Tc;
}
UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb){
	assert(Ta.status & Tb.status & Ta.HAVEELEM);
	UniTensor cTa = Ta;
	UniTensor cTb = Tb;
	return contract(cTa, cTb, true);
}
UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast){
	assert(Ta.status & Tb.status & Ta.HAVEELEM);
	if(&Ta == &Tb){
		UniTensor Ttmp = Tb;
		return contract(Ta, Ttmp, fast);
	}
	if(Ta.status & Ta.HAVEBOND && Tb.status & Ta.HAVEBOND){
		int AbondNum = Ta.bonds.size();
		int BbondNum = Tb.bonds.size();
		std::vector<int> oldLabelA = Ta.labels;
		std::vector<int> oldLabelB = Tb.labels;
		int oldRnumA = Ta.RBondNum;
		int oldRnumB = Tb.RBondNum;
		std::vector<int> newLabelA;
		std::vector<int> interLabel;
		std::vector<int> newLabelB;
		std::vector<int> markB(BbondNum, 0);
		std::vector<int> newLabelC;
		bool match;
		for(int a = 0; a < AbondNum; a++){
			match = false;
			for(int b = 0; b < BbondNum; b++)
				if(Ta.labels[a] == Tb.labels[b]){
					markB[b] = 1;
					interLabel.push_back(Ta.labels[a]);
					newLabelB.push_back(Tb.labels[b]);
          bool BOND_DIMENSION_MATCHED = Ta.bonds[a].dim() == Tb.bonds[b].dim();
          assert(BOND_DIMENSION_MATCHED);
					match = true;
					break;
				}
			if(!match){
				newLabelA.push_back(Ta.labels[a]);
				newLabelC.push_back(Ta.labels[a]);
			}
		}
		for(int a = 0; a < interLabel.size(); a++)
			newLabelA.push_back(interLabel[a]);
		for(int b = 0; b < BbondNum; b++)
			if(markB[b] == 0){
				newLabelB.push_back(Tb.labels[b]);
				newLabelC.push_back(Tb.labels[b]);
			}
		int conBond = interLabel.size();
		Ta.permute(newLabelA, AbondNum - conBond);
		Tb.permute(newLabelB, conBond);
		std::vector<Bond> cBonds;
		for(int i = 0; i < AbondNum - conBond; i++)
			cBonds.push_back(Ta.bonds[i]);
		for(int i = conBond; i < BbondNum; i++)
			cBonds.push_back(Tb.bonds[i]);
		UniTensor Tc(cBonds);
		if(cBonds.size())
			Tc.setLabel(newLabelC);
		Block blockA, blockB, blockC;
		std::map<Qnum,Block>::iterator it;
		std::map<Qnum,Block>::iterator it2;
		for(it = Ta.blocks.begin() ; it != Ta.blocks.end(); it++){
			if((it2 = Tb.blocks.find(it->first)) != Tb.blocks.end()){
				blockA = it->second;
				blockB = it2->second;
				blockC = Tc.blocks[it->first];
				assert(blockA.Rnum == blockC.Rnum && blockB.Cnum == blockC.Cnum && blockA.Cnum == blockB.Rnum);
				matrixMul(blockA.elem, blockB.elem, blockA.Rnum, blockB.Cnum, blockA.Cnum, blockC.elem, Ta.ongpu, Tb.ongpu, Tc.ongpu);
			}
		}
		Tc.status |= Tc.HAVEELEM;

		if(conBond == 0){	//Outer product
			int idx = 0;
			for(int i = 0; i < oldRnumA; i++){
				newLabelC[idx] = oldLabelA[i];
				idx++;
			}
			for(int i = 0; i < oldRnumB; i++){
				newLabelC[idx] = oldLabelB[i];
				idx++;
			}
			for(int i = oldRnumA; i < AbondNum; i++){
				newLabelC[idx] = oldLabelA[i];
				idx++;
			}
			for(int i = oldRnumB; i < BbondNum; i++){
				newLabelC[idx] = oldLabelB[i];
				idx++;
			}
			Tc.permute(newLabelC, oldRnumA + oldRnumB);
		}

		if(!fast){
			Ta.permute(oldLabelA, oldRnumA);
			Tb.permute(oldLabelB, oldRnumB);
		}
		return Tc;
	}
	else if(Ta.status & Ta.HAVEBOND)
		return Ta * Tb[0];
	else if(Tb.status & Tb.HAVEBOND)
		return Ta[0] * Tb;
	else
		return UniTensor(Ta[0] * Tb[0]);
}

UniTensor& UniTensor::operator*=(const UniTensor& uT){
	return *this = *this * uT;
}

UniTensor& UniTensor::operator*= (double a){
	assert(status & HAVEELEM);
	vectorScal(a, elem, m_elemNum, ongpu);
	return *this;
}

UniTensor operator*(const UniTensor& Ta, double a){
	assert(Ta.status & Ta.HAVEELEM);
	UniTensor Tb(Ta);
	vectorScal(a, Tb.elem, Tb.m_elemNum, Tb.ongpu);
	return Tb;
}
UniTensor otimes(const UniTensor & Ta, const UniTensor& Tb){
	UniTensor T1 = Ta;
	UniTensor T2 = Tb;
	//int label1[T1.bondNum()];
	//int label2[T2.bondNum()];
  std::vector<int> label1(T1.bondNum());
  std::vector<int> label2(T2.bondNum());
	for(int i = 0; i < T1.bondNum(); i++){
		if(i < T1.inBondNum())
			label1[i] = i;
		else
			label1[i] = T2.inBondNum() + i;
	}
	for(int i = 0; i < T2.bondNum(); i++){
		if(i < T2.inBondNum())
			label2[i] = i + T1.inBondNum();
		else
			label2[i] = i + T1.bondNum();
	}

	T1.setLabel(label1);
	T2.setLabel(label2);
	return contract(T1, T2, true);
}
};	/* namespace uni10 */

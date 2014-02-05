#include "../../datatype/QnumF.h"
using namespace uni10::datatype;
#include "../../data-structure/Block.h"
#include "../../datatype/Bond.h"
#include "../../numeric/myLapack.h"
#include "../Matrix.h"
#include "../SyTensor.h"

void SyTensor_t::operator+= (const SyTensor_t& Tb){
	assert(bonds == Tb.bonds);
	assert(blocks == Tb.blocks);
	vecAdd(Tb.elem, elem, elemNum);
}
SyTensor_t operator+(const SyTensor_t& Ta, const SyTensor_t& Tb){
	assert(Ta.bonds == Tb.bonds);
	assert(Ta.blocks == Tb.blocks);
	SyTensor_t Tc(Ta);
	vecAdd(Tb.elem, Tc.elem, Ta.elemNum);
	return Tc;
}
SyTensor_t& SyTensor_t::operator*= (SyTensor_t& Tb){
	return *this = *this * Tb;
}

SyTensor_t operator* (SyTensor_t& Ta, SyTensor_t& Tb){
	assert(Ta.status & Tb.status & Ta.INIT);
	assert(Ta.status & Tb.status & Ta.HAVELABEL);
	assert(Ta.status & Tb.status & Ta.HAVEELEM);
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
	Ta.reshape(newLabelA, AbondNum - conBond);
	Tb.reshape(newLabelB, conBond);
	std::vector<Bond_t> cBonds;
	for(int i = 0; i < AbondNum - conBond; i++)
		cBonds.push_back(Ta.bonds[i]);
	for(int i = conBond; i < BbondNum; i++)
		cBonds.push_back(Tb.bonds[i]);
	SyTensor_t Tc(cBonds);
	Tc.addLabel(newLabelC);
	Block_t blockA, blockB, blockC;
	std::map<Qnum_t,Block_t>::iterator it; 
	std::map<Qnum_t,Block_t>::iterator it2; 
	for(it = Ta.blocks.begin() ; it != Ta.blocks.end(); it++){
		if((it2 = Tb.blocks.find(it->first)) != Tb.blocks.end()){
			blockA = it->second;
			blockB = it2->second;
			blockC = Tc.blocks[it->first]; 
			assert(blockA.Rnum == blockC.Rnum && blockB.Cnum == blockC.Cnum && blockA.Cnum == blockB.Rnum);
			myDgemm(blockA.elem, blockB.elem, blockA.Rnum, blockB.Cnum, blockA.Cnum, blockC.elem);
		}
	}
	Tc.status |= Tc.HAVEELEM;
	Ta.reshape(oldLabelA, oldRnumA);
	Tb.reshape(oldLabelB, oldRnumB);
	return Tc;
}

void SyTensor_t::operator*= (double a){
	vecScal(a, elem, elemNum);
}

SyTensor_t operator*(const SyTensor_t& Ta, double a){
	SyTensor_t Tb(Ta);
	vecScal(a, Tb.elem, Tb.elemNum);
	return Tb;
}

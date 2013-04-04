#include "SyTensor.h"

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
void SyTensor_t::operator*= (SyTensor_t& Tb){
	*this = *this * Tb;
}
SyTensor_t operator* (SyTensor_t& Ta, SyTensor_t& Tb){
	assert(Ta.status & Tb.status & INIT);
	assert(Ta.status & Tb.status & HAVELABEL);
	assert(Ta.status & Tb.status & HAVEELEM);
	int AbondNum = Ta.bonds.size();
	int BbondNum = Tb.bonds.size();
	vector<int> newLabelA;
	vector<int> interLabel;
	vector<int> newLabelB;
	vector<int> markB(BbondNum, 0);
	vector<int> newLabelC;
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
	assert(Ta.blocks.size() == Tb.blocks.size());
	vector<Bond_t> cBonds;
	for(int i = 0; i < AbondNum - conBond; i++)
		cBonds.push_back(Ta.bonds[i]);
	for(int i = conBond; i < BbondNum; i++)
		cBonds.push_back(Tb.bonds[i]);
	SyTensor_t Tc(cBonds);
	cout << Ta;
	cout << Tb;
	Tc.addLabel(newLabelC);
	Block_t blockA, blockB, blockC;
	map<Qnum_t,Block_t>::iterator it; 
	for(it = Ta.blocks.begin() ; it != Ta.blocks.end(); it++){
		blockA = it->second;
		blockB = Tb.blocks[it->first];
		blockC = Tc.blocks[it->first]; 
		cout<<it->first<<endl;
		cout<<blockA<<endl;
		cout<<blockB<<endl;
		cout<<blockC<<endl;
		assert(blockA.Rnum == blockC.Rnum && blockB.Cnum == blockC.Cnum && blockA.Cnum == blockB.Rnum);
		myDgemm(blockA.elem, blockB.elem, blockA.Rnum, blockB.Cnum, blockA.Cnum, blockC.elem);
	}
	Tc.status |= HAVEELEM;
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

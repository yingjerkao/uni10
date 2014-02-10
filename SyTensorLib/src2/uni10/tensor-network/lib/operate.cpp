#include <uni10/tensor-network/UniTensor.h>
#include <uni10/numeric/uni10_lapack.h>
//using namespace uni10::datatype;
namespace uni10{
UniTensor& UniTensor::operator+= (const UniTensor& Tb){
	assert(bonds == Tb.bonds);
	assert(blocks == Tb.blocks);
	vecAdd(Tb.elem, elem, m_elemNum);
	return *this;
}
UniTensor operator+(const UniTensor& Ta, const UniTensor& Tb){
	assert(Ta.bonds == Tb.bonds);
	assert(Ta.blocks == Tb.blocks);
	UniTensor Tc(Ta);
	vecAdd(Tb.elem, Tc.elem, Ta.m_elemNum);
	return Tc;
}
UniTensor& UniTensor::operator*= (UniTensor& Tb){
	return *this = *this * Tb;
}

UniTensor operator* (UniTensor& Ta, UniTensor& Tb){
	assert(Ta.status & Tb.status & Ta.INIT);
	//assert(Ta.status & Tb.status & Ta.HAVELABEL);
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
	Ta.permute(newLabelA, AbondNum - conBond);
	Tb.permute(newLabelB, conBond);
	std::vector<Bond> cBonds;
	for(int i = 0; i < AbondNum - conBond; i++)
		cBonds.push_back(Ta.bonds[i]);
	for(int i = conBond; i < BbondNum; i++)
		cBonds.push_back(Tb.bonds[i]);
	UniTensor Tc(cBonds);
	Tc.addLabel(newLabelC);
	Block blockA, blockB, blockC;
	std::map<Qnum,Block>::iterator it; 
	std::map<Qnum,Block>::iterator it2; 
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
	Ta.permute(oldLabelA, oldRnumA);
	Tb.permute(oldLabelB, oldRnumB);
	return Tc;
}

UniTensor& UniTensor::operator*= (double a){
	vecScal(a, elem, m_elemNum);
	return *this;
}

UniTensor operator*(const UniTensor& Ta, double a){
	UniTensor Tb(Ta);
	vecScal(a, Tb.elem, Tb.m_elemNum);
	return Tb;
}
};	/* namespace uni10 */	

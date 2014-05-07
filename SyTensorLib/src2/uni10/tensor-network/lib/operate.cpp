#include <uni10/tensor-network/UniTensor.h>
#include <uni10/numeric/uni10_lapack.h>
//using namespace uni10::datatype;
namespace uni10{
UniTensor& UniTensor::operator+= (const UniTensor& Tb){
	assert(status & Tb.status & HAVEELEM);
	assert(bonds == Tb.bonds);
	assert(blocks == Tb.blocks);
	vecAdd(Tb.elem, elem, m_elemNum);
	return *this;
}
UniTensor operator+(const UniTensor& Ta, const UniTensor& Tb){
	assert(Ta.status & Tb.status & Ta.HAVEELEM);
	assert(Ta.bonds == Tb.bonds);
	assert(Ta.blocks == Tb.blocks);
	UniTensor Tc(Ta);
	vecAdd(Tb.elem, Tc.elem, Ta.m_elemNum);
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
	vecScal(a, elem, m_elemNum);
	return *this;
}

UniTensor operator*(const UniTensor& Ta, double a){
	assert(Ta.status & Ta.HAVEELEM);
	UniTensor Tb(Ta);
	vecScal(a, Tb.elem, Tb.m_elemNum);
	return Tb;
}
UniTensor outer(const UniTensor & Ta, const UniTensor& Tb){
	UniTensor T1 = Ta;
	UniTensor T2 = Tb;
	int label1[T1.bondNum()];
	int label2[T2.bondNum()];
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

	T1.addLabel(label1);
	T2.addLabel(label2);
	return contract(T1, T2, true);
}
};	/* namespace uni10 */	

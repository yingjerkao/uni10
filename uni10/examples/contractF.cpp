#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"

bool elemCmp(SyTensor_t& ten1, SyTensor_t& ten2){
	int n1 = ten1.getElemNum();
	int n2 = ten2.getElemNum();
	double diff;
	if(n1 == n2){
		for(int i = 0; i < n1; i++){
			diff = fabs(ten1.elem[i] - ten2.elem[i]);
			if(diff > 1E-6)
				return false;
		}
	}
	else
		return false;
	return true;
}

int main(){
	/*
	Qnum_t q_20(-2, 0);
	Qnum_t q_21(-2, 1);
	Qnum_t q_10(-1, 0, 0);
	Qnum_t q_11(-1, 0, 1);
	*/
	Qnum_t q00(0, 0, 0);
	Qnum_t q01(0, 0, 1);
	/*
	Qnum_t q10(1, 0, 0);
	Qnum_t q11(1, 0, 1);
	Qnum_t q20(2, 0);
	Qnum_t q21(2, 1);
	*/
	vector<Bond_t> bonds;
	vector<Qnum_t> qnums;
	/*
	qnums.push_back(q_10);
	qnums.push_back(q_11);
	*/
	qnums.push_back(q00);
	qnums.push_back(q00);
	qnums.push_back(q01);
	qnums.push_back(q01);
	/*
	qnums.push_back(q10);
	qnums.push_back(q11);
	*/
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	
	/*
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4
					  };
	*/

	SyTensor_t Ta(bonds);
	SyTensor_t Tb(bonds);
	Ta.orthoRand();
	Tb.orthoRand();
	int labelA[] = {-1, 4, 3, -3};
	int labelB[] = {-4, 3, -5, 4};
	Ta.addLabel(labelA);
	Tb.addLabel(labelB);
	SyTensor_t Tc = Ta * Tb;
	vector<_Swap> swaps = Tb.exSwap(Ta);
	for(int s = 0; s < swaps.size(); s++)
		cout<<"swap: "<<swaps[s].b1 <<" <-> " << swaps[s].b2<<endl;
	Tb.addGate(swaps);
	SyTensor_t Tcr = Tb * Ta;
	int rsp_ord[] = {-1, -4, -3, -5};
	Tc.reshape(rsp_ord, 4);
	Tcr.reshape(rsp_ord, 4);
	//cout<<Tc;
	assert(elemCmp(Tc, Tcr));
}


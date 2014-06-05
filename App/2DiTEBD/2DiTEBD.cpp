#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;

const int CHI = 2;

void bondcat(int bidx, map<Qnum, Matrix>& L, UniTensor& T);
void updateU(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
double measure(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& MPS, Network& meas);
double expectation(UniTensor& L, UniTensor& R, UniTensor& Op, Network& MPS, Network& meas);

int main(){
	/*** Initialization ***/
	Qnum q0(0);
	Bond bdi_chi(BD_IN, CHI);
	Bond bdo_chi(BD_OUT, CHI);
	Bond bdi_d(BD_IN, 2);
	Bond bdo_d(BD_OUT, 2);

	vector<Bond> bond2;
	bond2.push_back(bdi_chi);
	bond2.push_back(bdo_chi);

	vector<Bond> bond5;	// For ALa, BLb
	bond5.push_back(bdi_d);
	bond5.push_back(bdo_chi);
	bond5.push_back(bdo_chi);
	bond5.push_back(bdo_chi);
	bond5.push_back(bdo_chi);

	vector<Bond> bond4;
	bond4.push_back(bdi_d);
	bond4.push_back(bdi_d);
	bond4.push_back(bdo_d);
	bond4.push_back(bdo_d);

	UniTensor A(bond5, "A");
	UniTensor B(bond5, "B");
	UniTensor H(bond4, "H0");
	map<Qnum, Matrix> LU;
	map<Qnum, Matrix> LR;
	map<Qnum, Matrix> LD;
	map<Qnum, Matrix> LL;

	Matrix I_chi(CHI, CHI, true);
	I_chi[0] = 0.3;
	I_chi[1] = 0.7;
	LU[q0] = I_chi;
	LR[q0] = I_chi;
	LD[q0] = I_chi;
	LL[q0] = I_chi;

	A.set_zero();
	B.set_zero();
	A[0] = 1;
	B[0] = 1;

	for(int i = 0; i < A.elemNum(); i++){
		A[i] = log(i % 10 + 1);
		B[i] = log(i % 10 + 1);
	}


	double g = 0.7;
	double H_elem[] = { -1,  g,  g,  0, 
				  	     g,  1,  0,  g,
					     g,  0,  1,  g,
	 				     0,  g,  g, -1 };
	H.addRawElem(H_elem);
	UniTensor expH(H.bond(), "expH");
	expH.putBlock(q0, takeExp(-0.1, H.getBlock(q0)));
	UniTensor AL = A;
	UniTensor BL = B;
	bondcat(1, LU, AL);
	bondcat(2, LR, AL);
	bondcat(1, LD, BL);
	bondcat(2, LL, BL);

	Network iTEBD("2DiTEBD.net");
	Network updateA("updateA.net");
	Network MPS("MPS.net");
	Network meas("measure.net");

	updateU(AL, BL, LU, LR, LD, LL, expH, iTEBD, updateA);
	/*
	update(BLb, ALa, Lb, La, H, iTEBD, updateA);
	cout<<"E = "<<measure(ALa, BLb, La, Lb, H, MPS, meas)<<endl;
	*/

}

void bondcat(int bidx, map<Qnum, Matrix>& L, UniTensor& T){
	vector<int> label = T.label();
	int inBondNum = T.inBondNum();
	vector<int> per_label(label.size());
	for(int i = 0; i < label.size(); i++)
		if(i == bidx)
			per_label[0] = label[i];
		else if(i < bidx)
			per_label[i + 1] = label[i];
		else
			per_label[i] = label[i];
	T.permute(per_label, 1);
	for(map<Qnum,Matrix>::const_iterator it=L.begin(); it!=L.end(); ++it)
		T.putBlock(it->first, it->second * T.getBlock(it->first));
	T.permute(label, inBondNum);
}

void bondrm(int bidx, map<Qnum, Matrix>& L, UniTensor& T){
	map<Qnum, Matrix> invL;
	for(map<Qnum,Matrix>::iterator it=L.begin(); it!=L.end(); ++it){
		Matrix mat(it->second.row(), it->second.col(), true);
		for(int i = 0; i < it->second.elemNum(); i++){
			if(it->second[i] == 0)
				mat[i] = 0;
			else
				mat[i] = 1 / it->second[i];
		}
		invL[it->first] = mat;
	}
	bondcat(bidx, invL, T);
}

void updateU(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	Qnum q0(0);
	UniTensor A = AL;
	UniTensor B = BL;
	bondcat(4, LR, B);
	iTEBD.putTensor("A", &A);
	iTEBD.putTensor("B", &B);
	iTEBD.putTensor("expH", &expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta = C;
	bondcat(2, LD, Theta);
	bondcat(3, LL, Theta);

	vector<Matrix> rets = Theta.getBlock(q0).svd();
	int dim = CHI < rets[1].row() ? CHI : rets[1].row();
	Matrix lambda(dim, dim, true);
	for(int i = 0; i < dim; i++)
		lambda[i] = rets[1][i];
	double norm = lambda.norm();
	lambda *= (1 / norm);
	LU[q0] = lambda;
	Bond bdi(BD_IN, dim);
	vector<Bond> bond5 = B.bond();
	vector<Bond> per_bond5(5);
	int order[] = {3, 0, 4, 1, 2};
	Bond newBond(BD_IN, dim);
	per_bond5[0] = newBond;
	for(int i = 1; i < 5; i++)
		per_bond5[i] = bond5[order[i]];
	B.assign(per_bond5);
	B.addLabel(order);
	B.elemSet(q0, rets[2].elem());
	int per_order[] = {0, 1, 2, 3, 4};
	B.permute(per_order, 1);

	updateA.putTensor("B", &B);
	updateA.putTensor("C", &C);
	AL = updateA.launch();
	AL *= (1 / norm);
	BL = B;
	bondrm(4, LR, BL);
	cout<<AL;
	cout<<BL;
}

double measure(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& MPS, Network& meas){
	Qnum q0(0);
	UniTensor A1(ALa);
	A1.permute(1);
	A1.putBlock(q0, Lb[q0] * A1.getBlock(q0));
	A1.permute(2);

	UniTensor B1(BLb);
	B1.permute(1);
	B1.putBlock(q0, La[q0] * B1.getBlock(q0));
	B1.permute(2);
	
	double val = expectation(A1, BLb, Op, MPS, meas);
	val += expectation(B1, ALa, Op, MPS, meas);
	return val / 2;
}

double expectation(UniTensor& L, UniTensor& R, UniTensor& Op, Network& MPS, Network& meas){
	Qnum q0(0);
	MPS.putTensor("L", &L);
	MPS.putTensor("R", &R);
	UniTensor psi = MPS.launch();
	double norm = psi.getBlock(q0).norm();
	norm *= norm;
	meas.putTensor("bra", &psi);
	meas.putTensor("ket", &psi);
	meas.putTensor("Op", &Op);
	return meas.launch()[0] / norm;
}

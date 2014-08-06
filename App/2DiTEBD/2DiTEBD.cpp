#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
const int N = 50;
const int CHI = 2;

void bondcat(int bidx, map<Qnum, Matrix>& L, UniTensor& T);
void updateU(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
void updateR(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
void updateD(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
void updateL(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
double measure(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& Ob, Network& measure_net, Network& norm_net);
double measure2(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& Ob, Network& iTEBD);
double bond_expectation(UniTensor& ALL, UniTensor& BLL, UniTensor& Ob, Network& measure_net, Network& norm_net);

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

	vector<Bond> bond5;	// For AL, BL
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
  I_chi.randomize();
	LU[q0] = I_chi;
	LR[q0] = I_chi;
	LD[q0] = I_chi;
	LL[q0] = I_chi;

	A.randomize();
	B.randomize();

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

	Network iTEBD_V("2DiTEBD_V.net");
	Network updateA_V("updateA_V.net");
	Network iTEBD_H("2DiTEBD_H.net");
	Network updateA_H("updateA_H.net");
	Network measure_net("measure.net");
	Network norm_net("norm.net");

  for(int step = 0; step < N; step++){
	  updateU(AL, BL, LU, LR, LD, LL, expH, iTEBD_V, updateA_V);
  	updateR(AL, BL, LU, LR, LD, LL, expH, iTEBD_H, updateA_H);
  	updateD(AL, BL, LU, LR, LD, LL, expH, iTEBD_V, updateA_V);
  	updateL(AL, BL, LU, LR, LD, LL, expH, iTEBD_H, updateA_H);
    cout<<"E = "<<setprecision(9)<<measure(AL, BL, LU, LR, LD, LL, H, measure_net, norm_net)<<endl;
  }
  cout<<"E = "<<setprecision(9)<<measure2(AL, BL, LU, LR, LD, LL, expH, iTEBD_V)<<endl;

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
  Matrix blk = B.getBlock(q0);
  blk.setElem(rets[2].getElem());
	B.putBlock(q0, blk);
	int back_order[] = {0, 1, 2, 3, 4};
	B.permute(back_order, 1);

	updateA.putTensor("B", &B);
	updateA.putTensor("C", &C);
	AL = updateA.launch();
	AL *= (1 / norm);
  AL.addLabel(back_order);
	BL = B;
	bondrm(4, LR, BL);

}

void updateR(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	Qnum q0(0);
	UniTensor A = AL;
	UniTensor B = BL;
	bondcat(3, LU, B);
	iTEBD.putTensor("A", &A);
	iTEBD.putTensor("B", &B);
	iTEBD.putTensor("expH", &expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta = C;
	bondcat(1, LD, Theta);
	bondcat(2, LL, Theta);

	vector<Matrix> rets = Theta.getBlock(q0).svd();
	int dim = CHI < rets[1].row() ? CHI : rets[1].row();
	Matrix lambda(dim, dim, true);
	for(int i = 0; i < dim; i++)
		lambda[i] = rets[1][i];
	double norm = lambda.norm();
	lambda *= (1 / norm);
	LR[q0] = lambda;
	Bond bdi(BD_IN, dim);
	vector<Bond> bond5 = B.bond();
	vector<Bond> per_bond5(5);
	int order[] = {4, 0, 1, 2, 3};
	Bond newBond(BD_IN, dim);
	per_bond5[0] = newBond;
	for(int i = 1; i < 5; i++)
		per_bond5[i] = bond5[order[i]];
	B.assign(per_bond5);
	B.addLabel(order);
  Matrix blk = B.getBlock(q0);
  blk.setElem(rets[2].getElem());
	B.putBlock(q0, blk);
	int back_order[] = {0, 1, 2, 3, 4};
	B.permute(back_order, 1);

	updateA.putTensor("B", &B);
	updateA.putTensor("C", &C);
	AL = updateA.launch();
	AL *= (1 / norm);
  AL.addLabel(back_order);
	BL = B;
	bondrm(3, LU, BL);
}

void updateD(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	updateU(BL, AL, LD, LL, LU, LR, expH, iTEBD, updateA);
}

void updateL(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	updateR(BL, AL, LD, LL, LU, LR, expH, iTEBD, updateA);
}

double measure2(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& expH, Network& iTEBD){
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

  UniTensor Theta2 = Theta;
  UniTensor val = Theta * Theta2;
  double delta = 0.1;
  return -log(val[0])/delta/2;

}

double measure(UniTensor& AL, UniTensor& BL, map<Qnum, Matrix>& LU, map<Qnum, Matrix>& LR, map<Qnum, Matrix>& LD, map<Qnum, Matrix>& LL, UniTensor& Ob, Network& measure_net, Network& norm_net){
  UniTensor ALL = AL;
  UniTensor BLL = BL;
	bondcat(3, LD, ALL);
	bondcat(4, LL, ALL);
	bondcat(3, LU, BLL);
	bondcat(4, LR, BLL);

  double val = 0;

  // measure up bond
  UniTensor BLL1 = BLL;
  bondrm(3, LU, BLL1);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);

  // measure right bond
  BLL1 = BLL;
  int rotateR_label[] = {0, 2, 3, 4, 1};
  ALL.permute(rotateR_label, 1);
  BLL1.permute(rotateR_label, 1);
  bondrm(3, LR, BLL1);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);

  // measure down bond
  BLL1 = BLL;
  int rotateD_label[] = {0, 3, 4, 1, 2};
  ALL.permute(rotateD_label, 1);
  BLL1.permute(rotateD_label, 1);
  bondrm(3, LD, BLL1);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);

  // measure down bond
  BLL1 = BLL;
  int rotateL_label[] = {0, 4, 1, 2, 3};
  ALL.permute(rotateL_label, 1);
  BLL1.permute(rotateL_label, 1);
  bondrm(3, LL, BLL1);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);
  return val / 4;
}

double bond_expectation(UniTensor& ALL, UniTensor& BLL, UniTensor& Ob, Network& measure_net, Network& norm_net){
	measure_net.putTensor("ALL", &ALL);
	measure_net.putTensor("BLL", &BLL);
	measure_net.putTensorT("ALLT", &ALL);
	measure_net.putTensorT("BLLT", &BLL);
	measure_net.putTensorT("Ob", &Ob);
	UniTensor val = measure_net.launch();

	norm_net.putTensor("ALL", &ALL);
	norm_net.putTensor("BLL", &BLL);
	norm_net.putTensorT("ALLT", &ALL);
	norm_net.putTensorT("BLLT", &BLL);
	UniTensor norm = norm_net.launch();
  return val[0] / norm[0];
}

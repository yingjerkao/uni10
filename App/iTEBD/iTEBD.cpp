#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

const int CHI = 20;

void update(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& U, Network& iTEBD, Network& updateA);
double measure(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& MPS, Network& meas);
double expectation(UniTensor& L, UniTensor& R, UniTensor& Op, Network& MPS, Network& meas);

int main(){
  // Define the parameters of the model / simulation
  const double J = 1.0;
  const double g = 0.5;
  int d = 2;
  double delta = 0.01;
  int N = 1000;

	/*** Initialization ***/
	Qnum q0(0);
	Bond bdi_chi(BD_IN, CHI);
	Bond bdo_chi(BD_OUT, CHI);
	Bond bdi_d(BD_IN, d);
	Bond bdo_d(BD_OUT, d);

	vector<Bond> bond2;
	bond2.push_back(bdi_chi);
	bond2.push_back(bdo_chi);

	vector<Bond> bond3;
	bond3.push_back(bdi_chi);
	bond3.push_back(bdi_chi);
	bond3.push_back(bdo_d);

	vector<Bond> bond4;
	bond4.push_back(bdi_d);
	bond4.push_back(bdi_d);
	bond4.push_back(bdo_d);
	bond4.push_back(bdo_d);

	UniTensor ALa(bond3, "ALa");
	UniTensor BLb(bond3, "BLb");
	UniTensor H(bond4, "H0");
	map<Qnum, Matrix> La;
	map<Qnum, Matrix> Lb;

	Matrix I_chi(CHI, CHI, true);
	for(int i = 0; i < I_chi.row(); i++)
		I_chi[i] = 1;
	La[q0] = I_chi;
	Lb[q0] = I_chi;
	srand(time(NULL));
	I_chi.randomize();
	cout<<I_chi;

	ALa.set_zero();
	BLb.set_zero();
	for(int i = 0; i < ALa.elemNum(); i++){
		ALa[i] = i * 0.1;
		BLb[i] = i * 0.1;
	}

	double H_elem[] = {
     J, -g/2, -g/2,    0,
  -g/2,   -J,    0, -g/2,
  -g/2,    0,   -J, -g/2,
	 0, -g/2, -g/2,    J};

	H.addRawElem(H_elem);
	UniTensor U(H.bond(), "U");
	U.putBlock(q0, takeExp(-delta, H.getBlock(q0)));

	Network iTEBD("iTEBD.net");
	Network updateA("updateA.net");
	Network MPS("MPS.net");
	Network meas("measure.net");

  for(int step = 0; step < N; step++){
	  update(ALa, BLb, La, Lb, U, iTEBD, updateA);
	  update(BLb, ALa, Lb, La, U, iTEBD, updateA);
  }
  cout<<"E = "<<setprecision(12)<<measure(ALa, BLb, La, Lb, H, MPS, meas)<<endl;

}

void update(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& expH, Network& iTEBD, Network& updateA){
	Qnum q0(0);
	iTEBD.putTensor("ALa", &ALa);
	iTEBD.putTensor("BLb", &BLb);
	iTEBD.putTensor("expH", &expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta(C.bond(), "Theta");
	Theta.putBlock(q0, Lb[q0] * C.getBlock(q0));
	Theta.permute(2);
	vector<Matrix> rets = Theta.getBlock(q0).svd();
	int dim = CHI < rets[1].row() ? CHI : rets[1].row();
	Matrix lambda(dim, dim, true);
	for(int i = 0; i < dim; i++)
		lambda[i] = rets[1][i];
	double norm = lambda.norm();
	lambda *= (1 / norm);
	La[q0] = lambda;
	Bond bdi(BD_IN, dim);
	vector<Bond> bond3 = ALa.bond();
	bond3[0] = bdi;
	BLb.assign(bond3);
	BLb.elemSet(q0, rets[2].elem());
	updateA.putTensor("BLb", &BLb);
	updateA.putTensor("C", &C);
	ALa = updateA.launch();
	ALa *= (1 / norm);
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

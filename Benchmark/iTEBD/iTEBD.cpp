#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>
#include <stdlib.h>
#include <time.h>

int CHI = 20;

void update(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& U, Network& iTEBD, Network& updateA);
double measure(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& MPS, Network& meas);
double expectation(UniTensor& L, UniTensor& R, UniTensor& Op, Network& MPS, Network& meas);
double measure2(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& Lb, UniTensor& expH, Network & iTEBD, double delta);

int main(int argc, char* argv[]){
	// Define the parameters of the model / simulation
	const double J = 1.0;
	const double g = 0.5;
	double delta = 0.01;
	int N = 10;
	if(argc > 1){
		sscanf(argv[1], "%d",&CHI);
	}
	if(argc > 2){
		sscanf(argv[2], "%d",&N);
	}
	if(CHI < 20 || CHI > 2000){
		cout<<"Fatal Errod\n";
		exit(0);
	}

	/*** Initialization ***/
	Qnum q0(0);
	UniTensor H("heisenberg.ham");
	vector<Bond> bondH = H.bond();
	UniTensor U(bondH, "U");
	U.putBlock(q0, takeExp(-delta, H.getBlock(q0)));

	Bond bdi_chi(BD_IN, CHI);
	vector<Bond> bond3;
	bond3.push_back(bdi_chi);
	bond3.push_back(bdi_chi);
	bond3.push_back(bondH[2]);

	UniTensor ALa(bond3, "ALa");
	UniTensor BLb(bond3, "BLb");
	map<Qnum, Matrix> La;
	map<Qnum, Matrix> Lb;

	Matrix I_chi(CHI, CHI, true);
	I_chi.randomize();
	La[q0] = I_chi;
	I_chi.randomize();
	Lb[q0] = I_chi;

	ALa.randomize();
	BLb.randomize();

	Network iTEBD("iTEBD.net");
	Network updateA("updateA.net");
	Network MPS("MPS.net");
	Network meas("measure.net");

	clock_t t;
	t = clock();
	for(int step = 0; step < N; step++){
		printf("step = %d\n", step);
		update(ALa, BLb, La, Lb, U, iTEBD, updateA);
		update(BLb, ALa, Lb, La, U, iTEBD, updateA);
	}
	t = clock() - t;
	printf ("It took %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	iTEBD.profile();
	updateA.profile();
	H.profile();
	//cout<<"E = "<<setprecision(12)<<measure2(ALa, BLb, Lb, U, iTEBD, delta)<<endl;
	//cout<<"E = "<<setprecision(12)<<measure(ALa, BLb, La, Lb, H, MPS, meas)<<endl;
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
	double norm = rets[1].resize(dim, dim).norm();
	rets[1] *= (1 / norm);
	La[q0] = rets[1];
	Bond bdi(BD_IN, dim);
	vector<Bond> bond3 = ALa.bond();
	bond3[0] = bdi;
	BLb.assign(bond3);
	Matrix blk = BLb.getBlock(q0);
	blk.setElem(rets[2].getElem(), rets[2].isOngpu());
	BLb.putBlock(q0, blk);
	updateA.putTensor("BLb", &BLb);
	updateA.putTensor("C", &C);
	ALa = updateA.launch();
	ALa *= (1 / norm);
}

double measure2(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& Lb, UniTensor& expH, Network & iTEBD, double delta){
	Qnum q0(0);
	iTEBD.putTensor("ALa", &ALa);
	iTEBD.putTensor("BLb", &BLb);
	iTEBD.putTensor("expH", &expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta(C.bond(), "Theta");
	Theta.putBlock(q0, Lb[q0] * C.getBlock(q0));
	Theta.permute(2);
  UniTensor Theta2 = Theta;
  UniTensor val = Theta * Theta2;
  return -log(val[0]) / delta / 2;
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
	meas.putTensorT("ket", &psi);
	meas.putTensor("Op", &Op);
	return meas.launch()[0] / norm;
}

#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

const int CHI = 5;

void update(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& U, Network& iTEBD, Network& updateA);
double measure(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& MPS, Network& meas);
double expectation(UniTensor& L, UniTensor& R, UniTensor& Op, Network& MPS, Network& meas);
double measure2(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& Lb, UniTensor& expH, Network & iTEBD, double delta);
vector<Qnum> qnumA();
vector<Qnum> qnumB();

int main(){
  // Define the parameters of the model / simulation
  const double J = 1.0;
  const double g = 0.5;
  int d = 2;
  double delta = 0.01;
  int N = 1000;

	/*** Initialization ***/
  UniTensor H("../Hamiltonian/heisenberg_U1.ham");
  vector<Bond> bondH = H.bond();
	UniTensor U(H.bond(), "U");
  map<Qnum, Matrix> blocks = H.getBlocks();
  for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); it++)
	  U.putBlock(it->first, takeExp(-delta, it->second));

	Bond bdi_A(BD_IN, qnumA());
	Bond bdi_B(BD_IN, qnumB());
	Bond bdo_A(BD_OUT, qnumA());
	Bond bdo_B(BD_OUT, qnumB());
	vector<Bond> bondA;
	bondA.push_back(bdi_B);
	bondA.push_back(bdo_A);
	bondA.push_back(bondH[2]);
	vector<Bond> bondB;
	bondB.push_back(bdi_A);
	bondB.push_back(bdo_B);
	bondB.push_back(bondH[2]);

	UniTensor ALa(bondA, "ALa");
	UniTensor BLb(bondB, "BLb");
	ALa.randomize();
	BLb.randomize();

	map<Qnum, Matrix> La;
	map<Qnum, Matrix> Lb;
	srand(time(NULL));
  map<Qnum, int> degs = bdi_A.degeneracy();
  for(map<Qnum, int>::iterator it = degs.begin(); it != degs.end(); it++){
	  Matrix diag(it->second, it->second, true);
    diag.randomize();
    La[it->first] = diag;
  }
  degs = bdi_B.degeneracy();
  for(map<Qnum, int>::iterator it = degs.begin(); it != degs.end(); it++){
	  Matrix diag(it->second, it->second, true);
    diag.randomize();
    Lb[it->first] = diag;
  }

	Network iTEBD("iTEBD_U1.net");
	Network updateA("updateA_U1.net");
	Network MPS("MPS.net");
	Network meas("measure.net");

  for(int step = 0; step < N; step++){
	  update(ALa, BLb, La, Lb, U, iTEBD, updateA);
	  update(BLb, ALa, Lb, La, U, iTEBD, updateA);
  }
	cout<<"E = "<<setprecision(12)<<measure2(ALa, BLb, Lb, U, iTEBD, delta)<<endl;
}

vector<Qnum> qnumA(){
  Qnum q0(0);
  Qnum q2(2);
  Qnum qm2(-2);
  vector<Qnum> qnums;
  int num_0 = 12;
  int num_2 = 6;
  for(int i = 0; i < num_0; i++)
    qnums.push_back(q0);
  for(int i = 0; i < num_2; i++)
    qnums.push_back(q2);
  for(int i = 0; i < num_2; i++)
    qnums.push_back(qm2);
  return qnums;
}

vector<Qnum> qnumB(){
  Qnum q1(1);
  Qnum qm1(-1);
  Qnum q3(3);
  Qnum qm3(-3);
  vector<Qnum> qnums;
  int num_1 = 9;
  int num_3 = 3;
  for(int i = 0; i < num_1; i++)
    qnums.push_back(q1);
  for(int i = 0; i < num_1; i++)
    qnums.push_back(qm1);
  for(int i = 0; i < num_3; i++)
    qnums.push_back(q3);
  for(int i = 0; i < num_3; i++)
    qnums.push_back(qm3);
  return qnums;
}

void update(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& expH, Network& iTEBD, Network& updateA){
	iTEBD.putTensor("ALa", &ALa);
	iTEBD.putTensor("BLb", &BLb);
	iTEBD.putTensor("expH", &expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta(C.bond(), "Theta");
  map<Qnum, Matrix> blocks = C.getBlocks();
  for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); it++)
	  Theta.putBlock(it->first, Lb[it->first] * it->second);
	Theta.permute(2);
  blocks = Theta.getBlocks();
  map<Qnum, double> norm;
  //cout<<Theta;
  for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); it++){
    map<Qnum, Matrix>::iterator itL = La.find(it->first);
    if(itL == La.end())
      continue;
    vector<Matrix> rets = it->second.svd();
    int dim = itL->second.col();
    Matrix lambda(dim, dim, true);
    for(int i = 0; i < dim; i++)
      lambda[i] = rets[1][i];
    norm[it->first] = lambda.norm();
    La[it->first] = lambda * (1.0 / norm[it->first]);
    BLb.elemSet(it->first, rets[2].elem());
  }
	updateA.putTensor("BLb", &BLb);
	updateA.putTensor("C", &C);
	ALa = updateA.launch();
  int per[] = {-1, -2, 3};
  ALa.permute(per, 2);
  blocks = ALa.getBlocks();
  for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); it++)
	  ALa.putBlock(it->first, (1.0 /  norm[it->first]) * it->second);
  for(map<Qnum, double>::iterator it = norm.begin(); it != norm.end(); it++){
    cout<<it->first<<": "<<it->second<<endl;
  }
  int per_back[] = {-1, 3, -2};
  ALa.permute(per_back, 1);
}

double measure2(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& Lb, UniTensor& expH, Network & iTEBD, double delta){
	iTEBD.putTensor("ALa", &ALa);
	iTEBD.putTensor("BLb", &BLb);
	iTEBD.putTensor("expH", &expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta(C.bond(), "Theta");
  map<Qnum, Matrix> blocks = C.getBlocks();
  for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); it++)
	  Theta.putBlock(it->first, Lb[it->first] * it->second);
	Theta.permute(2);
  cout<<Theta;
  UniTensor Theta2 = Theta;
  UniTensor val = Theta * Theta2;
  cout<<val;
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
	meas.putTensor("ket", &psi);
	meas.putTensor("Op", &Op);
	return meas.launch()[0] / norm;
}

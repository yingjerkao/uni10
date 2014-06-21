#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>


void update(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& iTEBD, Network& updateA);
double measure(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& MPS, Network& meas);
double expectation(UniTensor& L, UniTensor& R, UniTensor& Op, Network& MPS, Network& meas);
void bondcat(UniTensor& T, map<Qnum, Matrix>& L, int bidx);
void bondrm(UniTensor& T, map<Qnum, Matrix>& L, int bidx);
vector<Qnum> qnumA();
vector<Qnum> qnumB();

int main(){
  // Define the parameters of the model / simulation
  int d = 2;
  double delta = 0.01;
  int N = 2000;

  // Generate the two-site time evolution operator
  UniTensor H("../Hamiltonian/heisenberg_U1.ham");
  vector<Bond> bondH = H.bond();
	UniTensor U(H.bond(), "U");
  map<Qnum, Matrix> blocks = H.getBlocks();
  for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); ++it)
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

  vector<UniTensor> Gs;
  UniTensor G_A(bondA);
  G_A.randomize();
  Gs.push_back(G_A);
  UniTensor G_B(bondB);
  G_B.randomize();
  Gs.push_back(G_B);


	map<Qnum, Matrix> L_A;
	map<Qnum, Matrix> L_B;
	//srand(time(NULL));
  map<Qnum, int> degs = bdi_A.degeneracy();
  for(map<Qnum, int>::iterator it = degs.begin(); it != degs.end(); ++it){
	  Matrix diag(it->second, it->second, true);
    diag.randomize();
    L_A[it->first] = diag;
  }
  degs = bdi_B.degeneracy();
  for(map<Qnum, int>::iterator it = degs.begin(); it != degs.end(); ++it){
	  Matrix diag(it->second, it->second, true);
    diag.randomize();
    L_B[it->first] = diag;
  }
  vector< map<Qnum, Matrix> > Ls;
  Ls.push_back(L_A);
  Ls.push_back(L_B);


  // Perform the imaginary time evolution alternating on A and B bonds
  UniTensor theta;
  for(int step = 0; step < N; step++){
    // Construct theta
    cout<<"=========================== STEP "<<step<<" ===============================\n";
    int A = step % 2;
    int B = (step + 1) % 2;
    bondcat(Gs[A], Ls[A], 1);
    bondcat(Gs[A], Ls[B], 0);
    bondcat(Gs[B], Ls[B], 1);

    int ordA[] = {-1, 3, 1};
    int ordB[] = {3, -3, 2};
    int ordU[] = {1, 2, -2, -4};
    Gs[A].addLabel(ordA);
    Gs[B].addLabel(ordB);
    U.addLabel(ordU);
    theta = contract(Gs[A], Gs[B], true);
    theta *= U;
    int ordTheta[] = {-1, -2, -3, -4};
    theta.permute(ordTheta, 2);


    Gs[A].transpose();
	  blocks = theta.getBlocks();
  	double norm = 0;
    for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); ++it){
      map<Qnum, Matrix>::iterator itL = Ls[A].find(it->first);
      if(itL == Ls[A].end()){
        continue;
      }
		  int chi = itL->second.col();

      // SVD
      vector<Matrix> svd = theta.getBlock(it->first).svd();
      if(step == 1){
        cout<<svd[0];
        cout<<svd[1];
        cout<<svd[2];
      }
      svd[0].transpose();

      // Truncate
      Matrix sv(chi, chi, true);
      for(int i = 0; i < chi; i++){
        sv[i] = svd[1][i];
        norm += sv[i] * sv[i];
      }
      Ls[A][it->first] = sv;
      Gs[A].elemSet(it->first, svd[0].elem());
      Gs[B].elemSet(it->first, svd[2].elem());
    }

    norm = sqrt(norm);
    Gs[A].transpose();

    if(step == 1){
      cout<<Gs[A];
      for(map<Qnum, Matrix>::iterator it = Ls[A].begin(); it != Ls[A].end(); ++it){
        cout<<it->second;
      }
      cout<<Gs[B];
      exit(0);
    }

    Gs[A].permute(ordA, 1);
    cout<<"NORM: "<<norm<<endl;
	  for(map<Qnum, Matrix>::iterator it = Ls[A].begin(); it != Ls[A].end(); ++it){
		  it->second *= 1.0 / norm;
      //cout<<it->second;
      for(int i = 0; i < it->second.col(); i++)
        cout<<it->second[i]<<", ";
      cout<<endl;
    }
    bondrm(Gs[A], Ls[B], 0);
    bondrm(Gs[B], Ls[B], 1);

  }
  UniTensor theta2 = theta;
  UniTensor val = theta * theta2;
  cout<<"E = "<<setprecision(12)<<-log(val[0])/delta/2<<endl;
}

void bondcat(UniTensor& T, map<Qnum, Matrix>& L, int bidx){
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

void bondrm(UniTensor& T, map<Qnum, Matrix>& L, int bidx){
	map<Qnum, Matrix> invL;
	for(map<Qnum,Matrix>::iterator it=L.begin(); it!=L.end(); ++it){
		Matrix mat(it->second.row(), it->second.col(), true);
		for(int i = 0; i < it->second.elemNum(); i++){
      /*
			if(it->second[i] == 0)
				mat[i] = 0;
			else
      */
				mat[i] = 1 / it->second[i];
		}
		invL[it->first] = mat;
	}
	bondcat(T, invL, bidx);
}

vector<Qnum> qnumA(){
  Qnum q0(0);
  Qnum q2(2);
  Qnum qm2(-2);
  vector<Qnum> qnums;
  int num_0 = 8;
  int num_2 = 4;
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
  int num_1 = 6;
  int num_3 = 2;
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


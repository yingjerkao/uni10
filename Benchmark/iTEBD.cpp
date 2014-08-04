#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>


void update(UniTensor& ALa, UniTensor& BLb, map<Qnum, Matrix>& La, map<Qnum, Matrix>& Lb, UniTensor& Op, Network& iTEBD, Network& updateA);
void bondcat(UniTensor& T, map<Qnum, Matrix>& L, int bidx);
void bondrm(UniTensor& T, map<Qnum, Matrix>& L, int bidx);

int main(){
  // Define the parameters of the model / simulation
  Qnum q0(0);
  const double J = 1.0;
  const double g = 0.5;
  int chi = 20;
  int d = 2;
  double delta = 0.01;
  int N = 2000;
  vector<UniTensor> Gs;

	Bond bdi_chi(BD_IN, chi);
	Bond bdo_chi(BD_OUT, chi);
	Bond bdi_d(BD_IN, 2);
	Bond bdo_d(BD_OUT, 2);
  // bondG has dim. [d, chi, chi];
	vector<Bond> bondG;
	bondG.push_back(bdi_chi);
	bondG.push_back(bdi_chi);
	bondG.push_back(bdo_d);

  UniTensor G_t(bondG);
  G_t.randomize();
  Gs.push_back(G_t);
  G_t.randomize();
  Gs.push_back(G_t);

  map<Qnum, Matrix> lambda;
  vector< map<Qnum, Matrix> > Ls(2, lambda);
  Matrix L_m(chi, chi, true);
  L_m.randomize();
  Ls[0][q0] = L_m;
  L_m.randomize();
  Ls[1][q0] = L_m;


  // Generate the two-site time evolution operator
  /*
	double H_elem[] = {\
     J, -g/2, -g/2,    0,\
  -g/2,   -J,    0, -g/2,\
  -g/2,    0,   -J, -g/2,\
  0, -g/2, -g/2,    J};
  */
  double H_elem[] = {\
	  1.0/4,      0,      0,     0,
		 0, -1.0/4,  1.0/2,     0,
		 0,  1.0/2, -1.0/4,     0,
		 0,      0,      0, 1.0/4};
  vector<Bond> bondH;
  bondH.push_back(bdi_d);
  bondH.push_back(bdi_d);
  bondH.push_back(bdo_d);
  bondH.push_back(bdo_d);
  UniTensor H(bondH);
  H.addRawElem(H_elem);
  UniTensor U(bondH);
  U.putBlock(q0, takeExp(-delta, H.getBlock(q0)));
  UniTensor theta(9);
  cout<<theta;

  // Perform the imaginary time evolution alternating on A and B bonds
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

    // SVD
    vector<Matrix> svd = theta.getBlock(q0).svd();

    // Truncate
	svd[1].resize(chi, chi);
    double norm = svd[1].norm();
    cout<<"NORM: "<<norm<<endl;
    svd[1] *= (1.0 / norm);
    Ls[A][q0] = svd[1];

    Gs[A].putBlock(q0, svd[0].resize(svd[0].row(), chi));
    Gs[B].putBlock(q0, svd[2].resize(chi, svd[2].col()));
    Gs[A].permute(ordA, 2);

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

  for(map<Qnum,Matrix>::const_iterator it=L.begin(); it!=L.end(); ++it){
    T.putBlock(it->first, it->second * T.getBlock(it->first));
  }
  T.permute(label, inBondNum);
}

void bondrm(UniTensor& T, map<Qnum, Matrix>& L, int bidx){
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
	bondcat(T, invL, bidx);
}

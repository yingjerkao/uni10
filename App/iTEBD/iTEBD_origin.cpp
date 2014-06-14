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

int main(){
  // Define the parameters of the model / simulation
  Qnum q0(0);
  const double J = 1.0;
  const double g = 0.5;
  int chi = 5;
  int d = 2;
  double delta = 0.01;
  int N = 1000;
  vector<UniTensor> Gs;

	Bond bdi_chi(BD_IN, chi);
	Bond bdo_chi(BD_OUT, chi);
	Bond bdi_d(BD_IN, 2);
	Bond bdo_d(BD_OUT, 2);
  // bondG has dim. [d, chi, chi];
	vector<Bond> bondG;
	bondG.push_back(bdi_d);
	bondG.push_back(bdo_chi);
	bondG.push_back(bdo_chi);

  UniTensor G_t(bondG);
  G_t.randomize();
  for(int i = 0; i < G_t.elemNum(); i++)
    G_t[i] = 0.1;
  Gs.push_back(G_t);
  G_t.randomize();
  for(int i = 0; i < G_t.elemNum(); i++)
    G_t[i] = 0.2;
  Gs.push_back(G_t);

  map<Qnum, Matrix> lambda;
  vector< map<Qnum, Matrix> > Ls(2, lambda);
  Matrix L_m(chi, chi, true);
  for(int i = 0; i < L_m.elemNum(); i++)
    L_m[i] = 0.5;
  Ls[0][q0] = L_m;
  for(int i = 0; i < L_m.elemNum(); i++)
    L_m[i] = 0.7;
  Ls[1][q0] = L_m;

  // Generate the two-site time evolution operator
	double H_elem[] = {\
     J, -g/2, -g/2,    0,\
  -g/2,   -J,    0, -g/2,\
  -g/2,    0,   -J, -g/2,\
	 	 0, -g/2, -g/2,    J};
	vector<Bond> bondH;
	bondH.push_back(bdi_d);
	bondH.push_back(bdi_d);
	bondH.push_back(bdo_d);
	bondH.push_back(bdo_d);
  UniTensor H(bondH);
  H.addRawElem(H_elem);
  UniTensor U(bondH);
  U.putBlock(q0, takeExp(-delta, H.getBlock(q0)));

  // Perform the imaginary time evolution alternating on A and B bonds
  Network theta_net("theta.net");
  for(int step = 0; step < N; step++){
    // Construct theta
    int A = step % 2;
    int B = (step + 1) % 2;
    bondcat(Gs[A], Ls[B], 1);
    bondcat(Gs[A], Ls[A], 2);
    bondcat(Gs[B], Ls[B], 2);
    theta_net.putTensor("U", &U);
    theta_net.putTensor("A", &Gs[A]);
    theta_net.putTensor("B", &Gs[B]);
    UniTensor theta = theta_net.launch();

    // SVD
    vector<Matrix> svd = theta.getBlock(q0).svd();
    svd[0].transpose();

    // Truncate
    Matrix sv(chi, chi, true);
    for(int i = 0; i < chi; i++)
      sv[i] = svd[1][i];
    double norm = sv.norm();
    sv *= (1.0 / norm);
    Ls[A][q0] = sv;
    Gs[A].permute(2).transpose().elemSet(q0, svd[0].elem());
    Gs[A].transpose().permute(1);
    bondrm(Gs[A], Ls[B], 1);
    int per_label1[] = {1, 0, 2};
    Gs[B].permute(per_label1, 1).elemSet(q0, svd[2].elem());
    int per_label2[] = {0, 1, 2};
    Gs[B].permute(per_label2, 1);
    bondrm(Gs[A], Ls[B], 2);
  }
  int A = 0, B = 1;
  bondcat(Gs[A], Ls[B], 1);
  bondcat(Gs[A], Ls[A], 2);
  bondcat(Gs[B], Ls[B], 2);
  theta_net.putTensor("U", &U);
  theta_net.putTensor("A", &Gs[A]);
  theta_net.putTensor("B", &Gs[B]);
  UniTensor theta = theta_net.launch();
  double sum = 0;
  cout<<theta;
  for(int i = 0; i < theta.elemNum(); i++)
    sum += theta[i] * theta[i];
  cout<<"E = "<<log(sum)/delta/2<<endl;


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
			if(it->second[i] == 0)
				mat[i] = 0;
			else
				mat[i] = 1 / it->second[i];
		}
		invL[it->first] = mat;
	}
	bondcat(T, invL, bidx);
}

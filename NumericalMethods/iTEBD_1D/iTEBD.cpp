#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>
#include "iTEBD_tools.cpp"
#include "hamiltonian.cpp"


int main(){
  // Define the parameters of the model / simulation
  int chi = 30;
  double delta = 0.02;
  int N = 2000;
  UniTensor H = Heisenberg();

	Bond bdi_chi(BD_IN, chi);
	Bond bdo_chi(BD_OUT, chi);
	vector<Bond> bondG;
	bondG.push_back(bdi_chi);
	bondG.push_back(bdo_chi);
	bondG.push_back(H.bond(2));

  vector<UniTensor> Gs;
  UniTensor gamma(bondG);
  gamma.randomize();
  Gs.push_back(gamma);
  gamma.randomize();
  Gs.push_back(gamma);

  vector<Matrix> Ls ;
  Matrix lambda(chi, chi, true);
  lambda.randomize();
  Ls.push_back(lambda);
  lambda.randomize();
  Ls.push_back(lambda);

  UniTensor U(H.bond());
  U.putBlock(takeExp(-delta, H.getBlock()));
  UniTensor theta;

  int ordA[] = {-1, 3, 1};
  int ordB[] = {3, -3, 2};
  int ordU[] = {1, 2, -2, -4};

  // Perform the imaginary time evolution alternating on A and B bonds
  for(int step = 0; step < N; step++){
    // Construct theta
    int A = step % 2;
    int B = (step + 1) % 2;
    bondcat(Gs[A], Ls[A], 1);
    bondcat(Gs[A], Ls[B], 0);
    bondcat(Gs[B], Ls[B], 1);

    Gs[A].setLabel(ordA);
		Gs[B].setLabel(ordB);
		U.setLabel(ordU);
		theta = contract(Gs[A], Gs[B], true);
		theta *= U;
		int ordTheta[] = {-1, -2, -3, -4};
		theta.permute(ordTheta, 2);

		// SVD
		vector<Matrix> svd = theta.getBlock().svd();

		// Truncate
		double norm = svd[1].resize(chi, chi).norm();
		svd[1] *= (1.0 / norm);
		Ls[A] = svd[1];
		Gs[A].putBlock(svd[0].resize(svd[0].row(), chi));
		Gs[B].putBlock(svd[2].resize(chi, svd[2].col()));
		Gs[A].permute(ordA, 1);
    bondrm(Gs[A], Ls[B], 0);
    bondrm(Gs[B], Ls[B], 1);
  }
  UniTensor val = theta * theta;
  cout<<"E = "<<setprecision(12)<<-log(val[0])/delta/2<<endl;
}

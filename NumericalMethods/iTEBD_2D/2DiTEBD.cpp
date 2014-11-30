#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include "hamiltonian.cpp"
#include "2DiTEBD_tools.cpp"

int main(){
	/*** Initialization ***/
  const int chi = 2;
  const int N = 200;
  double tau = 0.01;
	UniTensor H = transverseIsing(0.5, 0.7);
  cout<< H;

	Bond bdi_chi(BD_IN, chi);
	Bond bdo_chi(BD_OUT, chi);

	vector<Bond> bond2;
	bond2.push_back(bdi_chi);
	bond2.push_back(bdo_chi);

	vector<Bond> bond5;	// For AL, BL
	bond5.push_back(H.bond(0));
	bond5.push_back(bdo_chi);
	bond5.push_back(bdo_chi);
	bond5.push_back(bdo_chi);
	bond5.push_back(bdo_chi);

	UniTensor AL(bond5, "A");
	UniTensor BL(bond5, "B");
	AL.randomize();
	BL.randomize();

	Matrix I_chi(chi, chi, true);
  I_chi.randomize();
	Matrix LU = I_chi;
	Matrix LR = I_chi;
	Matrix LD = I_chi;
	Matrix LL = I_chi;

	bondcat(AL, LU, 1);
	bondcat(AL, LR, 2);
	bondcat(BL, LD, 1);
	bondcat(BL, LL, 2);

	UniTensor U(H.bond(), "U");
	U.putBlock(takeExp(-tau, H.getBlock()));

	Network iTEBD_V("2DiTEBD_V.net");
	Network updateA_V("updateA_V.net");
	Network iTEBD_H("2DiTEBD_H.net");
	Network updateA_H("updateA_H.net");
	Network measure_net("measure.net");
	Network norm_net("norm.net");

  for(int step = 0; step < N; step++){
	  updateU(AL, BL, LU, LR, LD, LL, U, iTEBD_V, updateA_V);
  	updateR(AL, BL, LU, LR, LD, LL, U, iTEBD_H, updateA_H);
  	updateD(AL, BL, LU, LR, LD, LL, U, iTEBD_V, updateA_V);
  	updateL(AL, BL, LU, LR, LD, LL, U, iTEBD_H, updateA_H);
    cout<<"E = "<<setprecision(9)<<measure(AL, BL, LU, LR, LD, LL, H, measure_net, norm_net)<<endl;
  }
  cout<<"E = "<<setprecision(9)<<measure2(AL, BL, LU, LR, LD, LL, U, iTEBD_V)<<endl;
  int ll[] = {-1, -2, -3, -4, -5};
  AL.setLabel(ll);
  AL.permute( 2);
  cout<<AL;
}

#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include "DMRG_tools.cpp"
#include "hamiltonian.cpp"

int main(){
	/*** Initialization ***/
  const int chi = 20;
  const int N = 100;
	UniTensor H0 = Heisenberg();

	vector<Bond> bond2;
	bond2.push_back(H0.bond(0));
	bond2.push_back(H0.bond(2));

	UniTensor Id(bond2, "Id");
	Id.identity();

  vector<UniTensor> HLs;
  vector<UniTensor> HRs;
	HLs.push_back(H0);
	HRs.push_back(H0);
  vector<UniTensor> As;
  vector<UniTensor> Bs;
  As.push_back(Id);
  Bs.push_back(Id);

	Network HLn("HL.net");
	Network HRn("HR.net");

	/*** END initilization ***/
  Matrix psi;
  Matrix pre_psi;
  double diff = 0;
  vector<Matrix> Ls(2);
  bool stable = false;

  double Ep = 0, E0;
	for(int l = 1; l < N; l++){
		UniTensor SB = combineH(H0, HLs[l - 1], HRs[l - 1]);
    int iter;
    UniTensor GS = findGS(SB, E0, psi, iter);

    UniTensor A, B;
    int D = updateMPS(GS, chi, A, B, Ls[1]);
    As.push_back(A);
    Bs.push_back(B);

    if(pre_psi.elemNum() == psi.elemNum()){
      diff = 1 - abs((psi * pre_psi.transpose())[0]);
      psi = trialState(A, B, Ls);
    }
    else if(D == chi)
      stable = true;

    cout<<"N = "<< 2 * (l+1) <<", D = " << D << setprecision(10) << ", E = " << E0  << ", e = " << E0 / (2 * (l + 1)) <<", iter = "<<iter<<", dE = "<<(E0 - Ep)/2<<", diff = "<<diff<<endl;
    UniTensor newHL, newHR;
    updateH(HLs[l - 1], HRs[l - 1], A, B, H0, H0, HLn, HRn, newHL, newHR);
    HLs.push_back(newHL);
    HRs.push_back(newHR);

    pre_psi = psi;
    Ls[0] = Ls[1];
    Ep = E0;
	}
}

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
  const int N = 20;
	//UniTensor H0 = theModel(1, 0, 0, 1.5, 0.1, 0);
	UniTensor H0 = Heisenberg(1);
  //UniTensor H0 = transverseIsing(0.5, 0);

	vector<Bond> bond2;
	bond2.push_back(H0.bond(0));
	bond2.push_back(H0.bond(2));

	UniTensor Id(bond2, "Id");
	Id.identity();

  vector<UniTensor> As;
  vector<UniTensor> Bs;
  vector<Matrix> Ls;
  vector<UniTensor> HLs;
  vector<UniTensor> HRs;
	HLs.push_back(H0);
	HRs.push_back(H0);
  As.push_back(Id);
  Bs.push_back(Id);
  Ls.push_back(Id.getBlock());

	Network HLn("HL.net");
	Network HRn("HR.net");

  Matrix psi;
	/*** END initilization ***/

  double Ep = 0, E0;
	for(int l = 1; l < N; l++){
		UniTensor SB = combineH(H0, HLs[l - 1], HRs[l - 1]);

    int iter;
    UniTensor GS = findGS(SB, E0, psi, iter);

    UniTensor A, B;
    Matrix L;
    int D = updateMPS(GS, chi, A, B, L);
    As.push_back(A);
    Bs.push_back(B);
    Ls.push_back(L);

    cout<<"N = "<< 2 * (l+1) <<", D = " << D << setprecision(10) << ", E = " << E0  << ", e = " << E0 / (2 * (l + 1)) <<", iter = "<<iter<<", dE = "<<(E0 - Ep)/2<<endl;
    UniTensor newHL, newHR;
    updateH(HLs[l - 1], HRs[l - 1], A, B, H0, H0, HLn, HRn, newHL, newHR);
    HLs.push_back(newHL);
    HRs.push_back(newHR);
    Ep = E0;
	}
  sweep(N, chi, N-1, 1, H0, HLs, HRs, As, Bs, Ls, HLn, HRn);

  Network normL("normL.net");
  Network normR("normR.net");
  Network expOb("expOb.net");
  cout<<H0;
  cout<<"expH = "<<mpsExp2s(As, Bs, Ls[N - 1], N, H0, normL, normR, expOb)<<endl;
  cout<<"mpsNorm = "<<mpsNorm(As, Bs, Ls[N - 1], normL, normR)<<endl;

  double up_elem[] = {0, 1, 0};
  double dn_elem[] = {0, 1.0, 0};
  vector<Bond> bond1(1, H0.bond(0));
  UniTensor up_st(bond1);
  UniTensor dn_st(bond1);
  up_st.setElem(up_elem);
  dn_st.setElem(dn_elem);
  cout<<up_st<<dn_st;

  vector<UniTensor> prdState;
  for(int i = 0; i < 2*N; i++){
    if(i % 2 == 0)
      prdState.push_back(up_st);
    else
      prdState.push_back(dn_st);
  }
  Network prdL("expPrdL.net");
  Network prdR("expPrdR.net");
  cout<<"expPrd = "<<mpsPrdS(As, Bs, Ls[N - 1], prdState, prdL, prdR)<<endl;
}

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

void sweep(int N, int chi, int range, int times, vector<UniTensor>& H0s, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn);
int main(){
	/*** Initialization ***/
  const int chi = 20;
  const int N = 20;
  vector<UniTensor> H0s;
  Matrix rand(1, 2*N-1);
  rand.randomize();
  rand.randomize();
  for(int i = 0; i < 2*N-1; i++){
    cout<<2 * rand[i]<<", ";
    H0s.push_back(Heisenberg(0.5, 2 * rand[i]));
  }
  cout<<endl;
  UniTensor H0 = Heisenberg(0.5, 1);

	vector<Bond> bond2;
	bond2.push_back(H0.bond(0));
	bond2.push_back(H0.bond(2));

	UniTensor Id(bond2, "Id");
	Id.identity();

  vector<UniTensor> HLs;
  vector<UniTensor> HRs;
	HLs.push_back(H0s[hidx(N, Left, 0)]);
	HRs.push_back(H0s[hidx(N, Right, 0)]);
  vector<UniTensor> As;
  vector<UniTensor> Bs;
  As.push_back(Id);
  Bs.push_back(Id);

	Network HLn("../HL.net");
	Network HRn("../HR.net");

  Matrix psi;
	/*** END initilization ***/

	for(int l = 1; l < N; l++){
		UniTensor SB = combineH(H0, HLs[l - 1], HRs[l - 1]);

    double E0;
    int iter;
    UniTensor GS = findGS(SB, E0, psi, iter);

    UniTensor A, B;
    int D = updateMPS(GS, chi, A, B);
    As.push_back(A);
    Bs.push_back(B);

    cout<<"N = "<< 2 * (l+1) <<", D = " << D << setprecision(10) << ", E = " << E0  << ", e = " << E0 / (2 * (l + 1)) <<", iter = "<<iter<<endl;
    UniTensor newHL, newHR;
    updateH(HLs[l - 1], HRs[l - 1], A, B, H0s[hidx(N, Left, l)], H0s[hidx(N, Right, l)], HLn, HRn, newHL, newHR);
    HLs.push_back(newHL);
    HRs.push_back(newHR);
	}
  cout<<"size = "<<HLs.size()<<endl;
  sweep(N, chi, N-1, 3, H0s, HLs, HRs, HLn, HRn);
}



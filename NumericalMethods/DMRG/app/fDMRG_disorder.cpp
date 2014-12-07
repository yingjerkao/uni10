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
  const int chi = 20;
  const int N = 20;
  Matrix rand(1, 2*N-1);

  double dJ = 0.2;
  int sNum = 200;
  srand((int)dJ * chi * N * 1773);
  for(int s = 0; s < sNum; s++){
    /*** Initialization ***/
    rand.randomize();
    rand *= dJ;
    vector<UniTensor> H0s;
    for(int i = 0; i < 2*N-1; i++){
      double J = 1 + rand[i] - (dJ/2);
      H0s.push_back(Heisenberg(0.5, J));
    }
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

    double E0;
    for(int l = 1; l < N; l++){
      UniTensor SB = combineH(H0, HLs[l - 1], HRs[l - 1]);

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
    cout<<setprecision(10)<<E0<<", ";
    cout<<sweep(N, chi, N-1, 2, H0s, HLs, HRs, HLn, HRn);
    cout<<endl;
  }
}



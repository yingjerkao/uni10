#include <iostream>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include "hamiltonian.h"

int main(){
/*
  cout<<Heisenberg_U1();
  cout<<Heisenberg(1);
  Matrix mx = matSx(1);
  cout<<mx;
  cout<<mx.eigh()[0];
  cout<<mx.eigh()[1];
  cout<<theModel(0.5, 0, 0, 1.5, 0.1, 0);
  cout<<theModel(1, 0, 0, 1.5, 0.1, 0);
*/
  UniTensor H = Heisenberg();
  vector<Bond> bond2;
  bond2.push_back(H.bond(0));
  bond2.push_back(H.bond(2));
  UniTensor Id(bond2);
  Id.identity();
  UniTensor H01 = otimes(otimes(H, Id), Id);
  UniTensor H30 = H01;
  int per_lab[] = {1,2,3,0,5,6,7,4};
  H30.permute(per_lab, 4);
  UniTensor H12 = otimes(Id, otimes(H, Id));
  UniTensor H23 = otimes(Id, otimes(Id, H));
//J-term
  UniTensor HJ = H01 + H30 + H12 + H23;
//Q-term
  UniTensor Id2(H.bond());
  Id2.identity();
  UniTensor Sij = (H + (-1.0/4.0) * Id2);
  UniTensor H0132 = otimes(Sij, Sij);
  UniTensor H0312 = H0132;
  H0312.permute(per_lab, 4);
  UniTensor HQ = H0132 + H0312;
  cout<<JQmodel(1, 7);
}

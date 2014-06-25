#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){

  Bond bd_border(BD_IN, 10);
  Bond bd_internal(BD_IN, 9);
  Bond bd_physical(BD_IN, 8);
  vector<Bond> C_bds;
  vector<Bond> T_bds;
  vector<Bond> A_bds;
  vector<Bond> O_bds;
  C_bds.push_back(bd_border); C_bds.push_back(bd_border);
  T_bds.push_back(bd_border); T_bds.push_back(bd_border); T_bds.push_back(bd_internal); T_bds.push_back(bd_internal);
  A_bds.push_back(bd_internal); A_bds.push_back(bd_internal); A_bds.push_back(bd_internal); A_bds.push_back(bd_internal); A_bds.push_back(bd_physical);
  O_bds.push_back(bd_physical); O_bds.push_back(bd_physical); O_bds.push_back(bd_physical); O_bds.push_back(bd_physical);
  UniTensor C(C_bds, "C");
  UniTensor T(T_bds, "T");
  UniTensor A(A_bds, "A");
  UniTensor O(O_bds, "O");

  //Network CTM("CTM_fixed.net");
  Network CTM("CTM_suggest.net");
  CTM.putTensor("C1", &C);
  CTM.putTensor("C2", &C);
  CTM.putTensor("C3", &C);
  CTM.putTensor("C4", &C);
  CTM.putTensor("T1a", &T);
  CTM.putTensor("T1b", &T);
  CTM.putTensor("T2a", &T);
  CTM.putTensor("T2b", &T);
  CTM.putTensor("T3a", &T);
  CTM.putTensor("T3b", &T);
  CTM.putTensor("T4a", &T);
  CTM.putTensor("T4b", &T);
  CTM.putTensor("A1", &A);
  CTM.putTensor("A1T", &A);
  CTM.putTensor("A2", &A);
  CTM.putTensor("A2T", &A);
  CTM.putTensor("B1", &A);
  CTM.putTensor("B1T", &A);
  CTM.putTensor("B2", &A);
  CTM.putTensor("B2T", &A);
  CTM.putTensor("O", &O);

  //cout<<C<<T<<A<<O;
  cout<<CTM;
  cout<<"Max: "<<CTM.max_tensor_elemNum()<<endl;
  cout<<"Sum: "<<CTM.sum_of_memory_usage()<<endl;
  cout<<"Usage: "<<CTM.memory_requirement()<<endl;
  CTM.profile();
  C.profile();
}

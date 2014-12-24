#include<Accelerate/Accelerate.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<iostream>
#include<iomanip>
#include<cassert>
#include<fstream>
#include<cmath>
#include<map>
using namespace std;
#include"uni10.hpp"
using namespace uni10;
#define d 2     //dimensions of physical bonds
#define D 4     //dimensions of virtual bonds
#define GR 3    //Geographical Representation
#define Htype 2 //2 difference Hamiltonian Representations( H-up; H-down )
#define LAMBDA_MIN 1.0E-12
#include"KagomeLib.cpp"

int main(){
  UniTensor H = HeisenbergOneHalf(2);
  Bond bdi_d(BD_IN, 2);
  Bond bdo_d(BD_OUT, 2);
  vector<Bond> bond2;
  bond2.push_back(bdi_d);
  bond2.push_back(bdo_d);

  UniTensor Id(bond2);
  Id.identity();
//  cout << Id << endl;

  UniTensor buf = otimes(H, Id);

  int i[] = {0,1,2,3,4,5};
  int a[] = {1,2,0,4,5,3};
  int b[] = {2,0,1,5,3,4};

  UniTensor U1 = buf.permute(a, 3);
  UniTensor U2 = buf.permute(b, 3);
  UniTensor U3 = buf.permute(i, 3);

  UniTensor H3 = U1 + U2 + U3;

  cout << H3 << endl;
}

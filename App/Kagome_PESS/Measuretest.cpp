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
#define GR 4    //Geographical Representation
#define Htype 1 //2 difference Hamiltonian Representations( H-up; H-down )
#define LAMBDA_MIN 1.0E-12
#include"KagomeLib.cpp"

int main(){

  UniTensor Tu = initTensor(d, D, GR);
  UniTensor Td = initTensor(d, D, GR);

  int n_site = GR;
  UniTensor H = HeisenbergOneHalf(n_site);
  UniTensor expH = Exp_dT(H, 0.01);

  Network measure("Obs.net");
   cout << measure << endl;
  measure.putTensor("Tu", &Tu);
  measure.putTensorT("TuT", &Tu);
  measure.putTensor("H", &H);
  UniTensor M = measure.launch();

//   cout << Tu.transpose() << endl;
//    cout << M << endl;
}

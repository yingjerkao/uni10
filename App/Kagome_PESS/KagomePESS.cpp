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
#define D 4    //dimensions of virtual bonds
#define GR 3    //Geographical Representation
#define Htype 2 //2 difference Hamiltonian Representations( H-up; H-down )
#define LAMBDA_MIN 1.0E-12
#include"KagomeLib.cpp"
#include"Measurement.cpp"
//#include"writeNetFile.cpp"

int main(){

	cout.setf(ios::showpoint);
	cout.precision(30);
	cout.setf(ios::fixed);

  Qnum q0(0);
  Bond bdi_d(BD_IN, d);
  Bond bdo_d(BD_OUT, d);
  Bond bdi_D(BD_IN, D);
  Bond bdo_D(BD_OUT, D);
/*
 *             0                                     1         2
 *             |                                      \       /
 *            / \                                      \-- --/
 *           /   \                                      \   /
 *          /-- --\                                      \ /
 *         /       \   -- Hamiltonian of                  |         -- Hamiltonian of
 *        1         2        Upware Triangle              0                 Downware Triangle
 */


  int n_site = GR;
  UniTensor Hu = HeisenbergOneHalf(3);
  UniTensor H = HeisenbergOneHalf(2);
  UniTensor expHu(Hu.bond());
  expHu.putBlock(takeExp(-0.001, Hu.getBlock()));

  //InitTensor U
  vector<Bond> bondU3;
  bondU3.push_back(bdi_d);
  bondU3.push_back(bdi_D);
  bondU3.push_back(bdo_D);

  UniTensor U0(bondU3);
  U0.orthoRand();
//  double normU = U0.getBlock(q0).norm();
//  Matrix U_elem(inBondDim(U0), outBondDim(U0));
//  U_elem = U0.getBlock(q0) * (1/normU);
//  U0.putBlock(q0, U_elem);

  UniTensor U1 = U0;
  UniTensor U2 = U0;
  UniTensor U3 = U0;
  UniTensor U4 = U0;

  vector<Bond> bondS3;
  bondS3.push_back(bdi_D);
  bondS3.push_back(bdo_D);
  bondS3.push_back(bdo_D);

  UniTensor S0(bondS3);
  S0.orthoRand();
//  double normS = S0.getBlock(q0).norm();
//  Matrix S_elem(inBondDim(S0), outBondDim(S0));
//  S_elem = S0.getBlock(q0) * (1/normS);
//  S0.putBlock(q0, S_elem);
//  cout << S0 << endl;
  map<string, UniTensor> As;
  As["U0"] = U0;
  As["U1"] = U1;
  As["U2"] = U2;
  As["S0"] = S0;


  //Permute U2 for Initializing And
/*
  vector<int> initlabelU2 = U2.label();
  int labelU2[] = {0,2,1};
  U2.permute(labelU2, 2);

  map<string, UniTensor> Adn;
  Adn["U0"] = U2;
  Adn["U1"] = U3;
  Adn["U2"] = U4;
  Adn["S0"] = S0;
*/

  //InitLambda La, Lb, Lc
  map<string, UniTensor> lAlpha = initLambda(D, GR);
  map<string, UniTensor> lBeta  = lAlpha;
  //Delcare subbonds type of U inBond.
  /*
  vector<Bond> bondU;
  bondU.push_back(bdi_d);
  bondU.push_back(bdi_D);

  */
/* Upware Triangle
     D[2]  d[1]                  d[-1] d[-3] d[-5]
      \   \                       |     |     |
              -- D[6]           -- -- -- -- -- --
          Td                    |     expH      |
              -- d[5]           -- -- -- -- -- --
       /  /                       |     |     |
     d[3]  D[4]                  d[1]  d[3]  d[5]
*/
/* Downware Trangle
     D[2]  d[1]                  d[-1] d[-5] d[-3]
      \   \                       |     |     |
              -- D[6]           -- -- -- -- -- --
          Td                    |     expH      |
              -- d[5]           -- -- -- -- -- --
       /  /                       |     |     |
     d[3]  D[4]                  d[1]  d[5]  d[3]
*/

    //Update PESS
  double Efu, Efd;

  Network simplexA("A.net");
  Network simplexB("B.net");

	for(int i = 0; i < 10000; i++){
      simpleUpdatePESS(true , As, lAlpha, lBeta, expHu, "KagomePESSnet/",simplexA);
      simpleUpdatePESS(false , As, lAlpha, lBeta, expHu, "KagomePESSnet/",simplexB);
      cout << "========================" << endl;
      cout << measureObs(lAlpha, As, "T", Hu , "H", "KagomePESSnet/") << endl;
  }
}

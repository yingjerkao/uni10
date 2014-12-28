#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include "PESS_tools.cpp"
#include "hamiltonian.cpp"

int main(){
  const int chi = 10;
  const int modeNum = 4;
  const double delta = 0.05;
  const int N = 10000;

  UniTensor squareH = JQmodel(-1, -1.4);//periodicHamiltonian(modeNum, Heisenberg());
  cout<<squareH;
  UniTensor expH(squareH.bond());
  expH.putBlock(takeExp(-delta, squareH.const_getBlock()));

  vector<Bond> bonds;
  bonds.push_back(squareH.bond(0));
  bonds.push_back(Bond(BD_IN, chi));
  bonds.push_back(Bond(BD_OUT, chi));
  vector<UniTensor>Us(modeNum, UniTensor(bonds));
  for(int i = 0; i < Us.size(); i++)
    Us[i].randomize();

  vector<Bond> bondC;
  bondC.push_back(Bond(BD_IN, chi));
  bondC.push_back(Bond(BD_OUT, chi));
  bondC.push_back(Bond(BD_OUT, chi));
  bondC.push_back(Bond(BD_OUT, chi));
  vector<UniTensor>Cs(2, UniTensor(bondC));
  for(int i = 0; i < Cs.size(); i++)
    Cs[i].randomize();

  vector<Matrix> Ls(2 * modeNum, Matrix(chi, chi, true));
  for(int i = 0; i < Ls.size(); i++)
    Ls[i].randomize();

  Network simplexUp("simplex4Out.net");
  Network simplexDn("simplex4In.net");
  Network stateUp("state4Out.net");
  Network stateDn("state4In.net");
  Network measure3("measure4Ob4.net");
  Network measure2("measure4Ob2.net");
  UniTensor H0 = Heisenberg();
  for(int n = 0; n < N; n++){
    simpleUpdate(true, Us, Cs, Ls, expH, simplexUp);
    cout<<measureObs(true, squareH, Us, Cs, Ls, stateUp, measure3) * 2.0 / 4<<", ";
    cout<<measureObs(true, H0, Us, Cs, Ls, stateUp, measure2) * 2<<", ";
    simpleUpdate(false, Us, Cs, Ls, expH, simplexDn);
    cout<<measureObs(false, squareH, Us, Cs, Ls, stateDn, measure3) * 2.0 / 4<<", ";
    cout<<measureObs(false, H0, Us, Cs, Ls, stateDn, measure2) * 2<<endl;
  }
  return 0;
}

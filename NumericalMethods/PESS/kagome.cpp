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

/*
 *             0                                     1         2
 *             |                                      \       /
 *            / \                                      \-- --/
 *           /   \                                      \   /
 *          /-- --\                                      \ /
 *         /       \   -- Hamiltonian of                  |         -- Hamiltonian of
 *        1         2        Upware Triangle              0                 Downware Triangle
 */

int main(){
  const int chi = 13;
  const int modeNum = 3;
  const double delta = 0.05;
  const int N = 10000;

  UniTensor triH = periodicHamiltonian(modeNum, Heisenberg());

  cout<<triH;
  UniTensor expH(triH.bond());
  expH.putBlock(takeExp(-delta, triH.const_getBlock()));

  vector<Bond> bonds;
  bonds.push_back(triH.bond(0));
  bonds.push_back(Bond(BD_IN, chi));
  bonds.push_back(Bond(BD_OUT, chi));
  vector<UniTensor>Us(3, UniTensor(bonds));
  for(int i = 0; i < Us.size(); i++)
    Us[i].randomize();

  bonds[0] = Bond(BD_IN, chi);
  bonds[1] = Bond(BD_OUT, chi);
  vector<UniTensor>Cs(2, UniTensor(bonds));
  for(int i = 0; i < Cs.size(); i++)
    Cs[i].randomize();

  vector<Matrix> Ls(6, Matrix(chi, chi, true));
  for(int i = 0; i < Ls.size(); i++)
    Ls[i].randomize();

  Network simplexUp("simplex3Out.net");
  Network simplexDn("simplex3In.net");
  Network stateUp("state3Out.net");
  Network stateDn("state3In.net");
  Network measure3("measure3Ob3.net");
  Network measure2("measure3Ob2.net");
  UniTensor H0 = Heisenberg();
  for(int n = 0; n < N; n++){
    simpleUpdate(true, Us, Cs, Ls, expH, simplexUp);
    cout<<measureObs(true, triH, Us, Cs, Ls, stateUp, measure3) * 2 / 3<<", ";
    cout<<measureObs(true, H0, Us, Cs, Ls, stateUp, measure2) * 2<<", ";
    simpleUpdate(false, Us, Cs, Ls, expH, simplexDn);
    cout<<measureObs(false, triH, Us, Cs, Ls, stateDn, measure3) * 2 / 3<<", ";
    cout<<measureObs(false, H0, Us, Cs, Ls, stateDn, measure2) * 2<<endl;
  }








  return 0;
}

#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){
  CMatrix mat(4, 5);
  mat.randomize();
  vector<CMatrix> rets = mat.svd();
  //cout<<rets[0];
  //cout<<rets[1];
  //cout<<rets[2];
  /*
  cout<<mat;
  cout<<rets[0]* rets[1] * rets[2];
  mat.orthoRand();
  //cout<<mat;
  CMatrix matT = mat;
  matT.cTranspose();
  //cout<<mat * matT;

  */

  CMatrix H(10, 10);
  H.randomize();
  CMatrix HH = H;
  HH.cTranspose();
  CMatrix S = H * HH + HH * H;
  vector<CMatrix>eigs = S.eigh();
  cout<<S;
  cout<<eigs[0];
  CMatrix U = eigs[1];
  CMatrix UT = U;
  UT.cTranspose();
  cout<<U;
  CMatrix psi(1, 10, eigs[1].getElem());
  psi.cTranspose();
  cout<<psi;

  double E0;
  Matrix rev(S.col(), 1);
  rev.randomize();
  CMatrix ev = rev;
  S.lanczosEigh(E0, ev);
  cout<<E0<<", "<<ev.norm()<<endl;
  cout<<ev;
  complex<double> r = psi[0] / ev[0];
  cout<<"ratio = "<<r<<endl;
  cout<<ev * r;

}

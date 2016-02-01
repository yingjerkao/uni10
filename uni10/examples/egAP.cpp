#include <iostream>
#include <assert.h>
#include <map>
#include <time.h>
#include "uni10.hpp"
#include "uni10/numeric/arpack/uni10_arpack.h"
using namespace uni10;

int main(){
  Matrix mat(RTYPE, 4, 5);
  mat.randomize();
  std::vector<Matrix> rets = mat.svd();

  Matrix H(RTYPE, 10, 10);
  H.randomize();
  Matrix HH = H;
  HH.transpose();
  Matrix S = H * HH + HH * H;
  std::vector<Matrix> eigs = S.eigh();
  std::cout << "Eigenvalues" << std::endl;
  std::cout << eigs[0] << std::endl;
  Matrix U = eigs[1];
  std::cout << "Eigenvectors" << std::endl;
  U.transpose();
  std::cout << U << std::endl;

  double E0;
  Matrix rev(RTYPE, S.col(), 1);
  rev.randomize();
  Matrix ev = rev;
  lanczosEigh(S, E0, ev);
  std::cout << "Eigenvalue from arpack" << std::endl;
  std::cout << E0 << std::endl;
  std::cout << "GS Eigenvector from arpack" << std::endl;
  std::cout << ev << std::endl;

  std::cout << "Complex Example" << std::endl;
  Matrix CH(CTYPE, 10, 10);
  CH.randomize();
  Matrix CHH = CH;
  CHH.cTranspose();
  Matrix CS = CH * CHH + CHH * CH;
  std::vector<Matrix> Ceigs = CS.eigh();
  std::cout << "Eigenvalues" << std::endl;
  std::cout << Ceigs[0] << std::endl;
  Matrix CU = Ceigs[1];
  std::cout << "Eigenvectors" << std::endl;
  CU.cTranspose();
  std::cout << CU << std::endl;

  Matrix crev(CTYPE, CS.col(), 1);
  crev.randomize();
  // std::cout << CS << std::endl;
  lanczosEigh(CS, E0, crev);
  std::cout << "Eigenvalue from arpack" << std::endl;
  std::cout << E0 << std::endl;
  std::cout << "GS Eigenvector from arpack" << std::endl;
  std::cout << crev << std::endl;
}

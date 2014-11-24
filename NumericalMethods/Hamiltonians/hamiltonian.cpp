#include <cmath>
#include <assert.h>
#include <exception>
#include "hamiltonian.h"
using namespace std;
using namespace uni10;
const int CURRENT_SUPPORTED_SPIN_DIM = 2;

Matrix matSp(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  double* mat_elem = (double*)malloc(dim * dim * sizeof(double));
  memset(mat_elem, 0, dim * dim * sizeof(double));
  if(dim == 2){
    mat_elem[1] = 1;
  }
  Matrix sp(dim, dim, mat_elem);
  free(mat_elem);
  return sp;
}

Matrix matSm(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  double* mat_elem = (double*)malloc(dim * dim * sizeof(double));
  memset(mat_elem, 0, dim * dim * sizeof(double));
  if(dim == 2){
    mat_elem[2] = 1;
  }
  Matrix sm(dim, dim, mat_elem);
  free(mat_elem);
  return sm;
}

Matrix matSx(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  double* mat_elem = (double*)malloc(dim * dim * sizeof(double));
  memset(mat_elem, 0, dim * dim * sizeof(double));
  if(dim == 2){
    mat_elem[1] = 1;
    mat_elem[2] = 1;
  }
  Matrix sp(dim, dim, mat_elem);
  free(mat_elem);
  return sp;
}

Matrix matSz(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  double* mat_elem = (double*)malloc(dim * dim * sizeof(double));
  memset(mat_elem, 0, dim * dim * sizeof(double));
  if(dim == 2){
    mat_elem[0] = 0.5;
    mat_elem[3] = -0.5;
  }
  Matrix sz(dim, dim, mat_elem);
  free(mat_elem);
  return sz;
}

UniTensor Heisenberg(float spin){
  Matrix sp = matSp(spin);
  Matrix sm = matSm(spin);
  Matrix sz = matSz(spin);
  Matrix ham = otimes(sz, sz);
  ham += 0.5 * (otimes(sp, sm) + otimes(sm, sp));
  Bond bdi = spin_bond(spin, BD_IN);
  Bond bdo = spin_bond(spin, BD_OUT);
  vector<Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  UniTensor H(bonds, "Heisenberg");
  H.putBlock(ham);
  return H;
}

UniTensor Heisenberg_U1(float spin){
  Matrix sp = matSp(spin);
  Matrix sm = matSm(spin);
  Matrix sz = matSz(spin);
  Matrix ham = otimes(sz, sz);
  ham += 0.5 * (otimes(sp, sm) + otimes(sm, sp));
  const bool U1 = true;
  Bond bdi = spin_bond(spin, BD_IN, U1);
  Bond bdo = spin_bond(spin, BD_OUT, U1);
  vector<Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  UniTensor H(bonds, "Heisenberg");
  H.setRawElem(ham.getElem());
  return H;
}

UniTensor transverseIsing(float spin, float h){
  Matrix sx = matSx(spin);
  Matrix sz = matSz(spin);
  Matrix I(sx.row(), sx.col(), true);
  I.identity();
  Matrix ham = 2 * otimes(sz, sz);
  ham += otimes((h/2) * sx, I);
  ham += otimes(I, (h/2) * sx);
  Bond bdi = spin_bond(spin, BD_IN);
  Bond bdo = spin_bond(spin, BD_OUT);
  vector<Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  UniTensor H(bonds, "transverseIsing");
  H.putBlock(ham);
  return H;
}

void spin_check(float spin){
  if(!(spin > 0 && floor(2 * spin) == 2 * spin)){
      std::ostringstream err;
      err<<"The given spin is not half integer.";
      throw std::runtime_error(err.str());
  }
}

Bond spin_bond(float spin, bondType btype, bool U1){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(U1){
    vector<Qnum> qnums(dim);
    int halfint = true;
    if(spin == floor(spin))
      halfint = false;
    for(int i = 0; i < dim; i++){
      int s = spin - i;
      if(halfint){
        s = spin + 0.5 - i;
        if(s <= 0)
          s--;
      }
      qnums[i] = Qnum(s);
    }
    return Bond(btype, qnums);
  }
  else{
    return Bond(btype, dim);
  }
}


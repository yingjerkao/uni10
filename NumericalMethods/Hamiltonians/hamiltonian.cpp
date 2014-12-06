#include <cmath>
#include <assert.h>
#include <exception>
#include "hamiltonian.h"
using namespace std;
const int CURRENT_SUPPORTED_SPIN_DIM = 2;

uni10::Matrix matSp(float spin){
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
  uni10::Matrix sp(dim, dim, mat_elem);
  free(mat_elem);
  return sp;
}

uni10::Matrix matSm(float spin){
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
  uni10::Matrix sm(dim, dim, mat_elem);
  free(mat_elem);
  return sm;
}

uni10::Matrix matSx(float spin){
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
    mat_elem[1] = 0.5;
    mat_elem[2] = 0.5;
  }
  uni10::Matrix sp(dim, dim, mat_elem);
  free(mat_elem);
  return sp;
}

uni10::Matrix matSz(float spin){
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
  uni10::Matrix sz(dim, dim, mat_elem);
  free(mat_elem);
  return sz;
}

uni10::UniTensor Heisenberg(float spin, double J){
  uni10::Matrix sp = matSp(spin);
  uni10::Matrix sm = matSm(spin);
  uni10::Matrix sz = matSz(spin);
  uni10::Matrix ham = uni10::otimes(sz, sz);
  ham += 0.5 * (uni10::otimes(sp, sm) + uni10::otimes(sm, sp));
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor H(bonds, "Heisenberg");
  H.putBlock(ham);
  return J * H;
}

uni10::UniTensor Heisenberg_U1(float spin, double J){
  const bool U1 = true;
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN, U1);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT, U1);
  vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor H(bonds, "Heisenberg");
  H.setRawElem(Heisenberg(spin, J).getBlock().getElem());
  return H;
}

uni10::UniTensor transverseIsing(float spin, float h){
  uni10::Matrix sx = matSx(spin);
  uni10::Matrix sz = matSz(spin);
  uni10::Matrix I(sx.row(), sx.col(), true);
  I.identity();
  uni10::Matrix ham = uni10::otimes(2*sz, 2*sz); // otimes(sigma_z, sizga_z);
  ham += uni10::otimes((h/2) * 2*sx, I);
  ham += uni10::otimes(I, (h/2) * 2*sx);
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor H(bonds, "transverseIsing");
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

uni10::Bond spin_bond(float spin, uni10::bondType btype, bool U1){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(U1){
    vector<uni10::Qnum> qnums(dim);
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
      qnums[i] = uni10::Qnum(s);
    }
    return uni10::Bond(btype, qnums);
  }
  else{
    return uni10::Bond(btype, dim);
  }
}


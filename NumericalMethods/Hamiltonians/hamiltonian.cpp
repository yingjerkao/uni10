#include <cmath>
#include <assert.h>
#include <exception>
#include "hamiltonian.h"
using namespace std;
const int CURRENT_SUPPORTED_SPIN_DIM = 3;

uni10::Matrix matSp(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    double mat_elem[] = {\
      0, 1,\
      0, 0};
    return uni10::Matrix(dim, dim, mat_elem);
  }
  else if(dim == 3){
    double mat_elem[] = {\
      0.0, 1.0, 0.0,\
      0.0, 0.0, 1.0,\
      0.0, 0.0, 0.0};
    return sqrt(2) * uni10::Matrix(dim, dim, mat_elem);
  }

  return uni10::Matrix();
}

uni10::Matrix matSm(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    double mat_elem[] = {\
      0, 0,\
      1, 0};
    return uni10::Matrix(dim, dim, mat_elem);
  }
  else if(dim == 3){
    double mat_elem[] = {\
      0.0, 0.0, 0.0,\
      1.0, 0.0, 0.0,\
      0.0, 1.0, 0.0};
    return sqrt(2) * uni10::Matrix(dim, dim, mat_elem);
  }
  return uni10::Matrix();
}

uni10::Matrix matSx(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    double mat_elem[] = {\
      0,   0.5,\
      0.5, 0  };
    return uni10::Matrix(dim, dim, mat_elem);
  }
  else if(dim == 3){
    double mat_elem[] = {\
      0.0, 1.0, 0.0,\
      1.0, 0.0, 1.0,\
      0.0, 1.0, 0.0};
    return (1.0 / sqrt(2)) * uni10::Matrix(dim, dim, mat_elem);
  }
  return uni10::Matrix();
}

uni10::Matrix matSz(float spin){
  spin_check(spin);
  int dim = spin * 2 + 1;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    double mat_elem[] = {\
      0.5,  0,\
      0,   -0.5  };
    return uni10::Matrix(dim, dim, mat_elem);
  }
  else if(dim == 3){
    double mat_elem[] = {\
      1.0, 0.0,  0.0,\
      0.0, 0.0,  0.0,\
      0.0, 0.0, -1.0};
    return uni10::Matrix(dim, dim, mat_elem);
  }
  return uni10::Matrix();
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

// arXiv:1409.8616v1, eq(S1), page 6
uni10::UniTensor theModel(float spin, int i, double delta, double Dz, double hz, double dx){
  uni10::Matrix sx = matSx(spin);
  uni10::Matrix sz = matSz(spin);
  uni10::Matrix sp = matSp(spin);
  uni10::Matrix sm = matSm(spin);
  uni10::Matrix I(sz.row(), sz.col(), true);
  I.identity();
  uni10::Matrix ham = (1.0 + delta*pow(-1, i)) * (otimes(sz, sz) + (1.0/2.0) * (otimes(sp, sm) + otimes(sm, sp)))\
             + (Dz/2.0) * (otimes(sz*sz, I) + otimes(I, sz*sz))\
             + (hz/2.0) * ((pow(-1, i) * otimes(sz, I)) + (pow(-1, i+1) * otimes(I, sz)))\
             + dx * (otimes(sx, sz) + (-1) * otimes(sz, sx));

  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor H(bonds, "The model");
  H.putBlock(ham);
  return H;
}

uni10::UniTensor JQmodel(double J, double Q){
  float spin = 0.5;
  uni10::UniTensor H = Heisenberg(spin);
  vector<uni10::Bond> bond2;
  bond2.push_back(H.bond(0));
  bond2.push_back(H.bond(2));
  uni10::UniTensor Id(bond2);
  Id.identity();
  uni10::UniTensor H01 = uni10::otimes(uni10::otimes(H, Id), Id);
  uni10::UniTensor H30 = H01;
  int per_lab[] = {1,2,3,0,5,6,7,4};
  H30.permute(per_lab, 4);
  uni10::UniTensor H12 = uni10::otimes(Id, uni10::otimes(H, Id));
  uni10::UniTensor H23 = uni10::otimes(Id, uni10::otimes(Id, H));
//J-term
  uni10::UniTensor HJ = H01 + H30 + H12 + H23;
//Q-term
  uni10::UniTensor Id2(H.bond());
  Id2.identity();
  uni10::UniTensor Sij = (H + (-1.0/4.0) * Id2);
  uni10::UniTensor H0132 = otimes(Sij, Sij);
  uni10::UniTensor H0312 = H0132;
  H0312.permute(per_lab, 4);
  uni10::UniTensor HQ = H0132 + H0312;
  return ((-J) * HJ) + ((-Q) * HQ);
}

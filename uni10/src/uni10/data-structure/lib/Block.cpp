/****************************************************************************
*  @file Block.cpp
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Implementation file of Block class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#include <uni10/tensor-network/Matrix.h>
#ifndef UNI10_PURE_REAL
#include <uni10/tensor-network/CMatrix.h>
#endif


#ifndef UNI10_DTYPE
#define UNI10_DTYPE double
#endif
#ifndef UNI10_MATRIX
#define UNI10_MATRIX Matrix
#endif
#ifndef UNI10_BLOCK
#define UNI10_BLOCK Block
#endif

namespace uni10{
UNI10_BLOCK::UNI10_BLOCK(): Rnum(0), Cnum(0), diag(false), ongpu(false), m_elem(NULL){}
UNI10_BLOCK::UNI10_BLOCK(size_t _Rnum, size_t _Cnum, bool _diag): Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL){}
UNI10_BLOCK::UNI10_BLOCK(const UNI10_BLOCK& _b): Rnum(_b.Rnum), Cnum(_b.Cnum), diag(_b.diag), ongpu(_b.ongpu), m_elem(_b.m_elem){}

size_t UNI10_BLOCK::row()const{return Rnum;}
size_t UNI10_BLOCK::col()const{return Cnum;}
bool UNI10_BLOCK::isDiag()const{return diag;}
bool UNI10_BLOCK::isOngpu()const{return ongpu;}
size_t UNI10_BLOCK::elemNum()const{
  if(diag)
    return (Rnum < Cnum ? Rnum : Cnum);
  else
    return Rnum * Cnum;
}

UNI10_DTYPE* UNI10_BLOCK::getElem()const{return m_elem;}

UNI10_MATRIX UNI10_BLOCK::getDiag()const{
  try{
    if(diag)
      return *this;
    else{
      UNI10_MATRIX D(Rnum, Cnum, true, ongpu);
      ::uni10::getDiag(m_elem, D.getElem(), Rnum, Cnum, D.elemNum(), ongpu, D.isOngpu());
      return D;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::getDiag():");
    return Matrix();
  }
}

void UNI10_BLOCK::save(const std::string& fname)const{
  try{
    FILE *fp = fopen(fname.c_str(), "w");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in writing to file '"<<fname<<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UNI10_DTYPE* elem = m_elem;
    if(ongpu){
      elem = (UNI10_DTYPE*)malloc(elemNum() * sizeof(UNI10_DTYPE));
      elemCopy(elem, m_elem, elemNum() * sizeof(UNI10_DTYPE), false, ongpu);
    }
    fwrite(elem, sizeof(UNI10_DTYPE), elemNum(), fp);
    fclose(fp);
    if(ongpu)
      free(elem);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::save(std::string&):");
  }
}

UNI10_DTYPE UNI10_BLOCK::operator[](size_t idx)const{
  try{
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return getElemAt(idx, m_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::operator[](size_t):");
    return 0;
  }
}

UNI10_DTYPE UNI10_BLOCK::at(size_t r, size_t c)const{
  try{
    if(!((r < Rnum) && (c < Cnum))){
      std::ostringstream err;
      err<<"The input indices are out of range.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag){
      if(!(r == c && r < elemNum())){
        std::ostringstream err;
        err<<"The matrix is diagonal, there is no off-diagonal element.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      return getElemAt(r, m_elem, ongpu);
    }
    else
      return getElemAt(r * Cnum + c, m_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::at(size_t, size_t):");
    return 0;
  }
}

std::vector<UNI10_MATRIX> UNI10_BLOCK::eigh()const{
  std::vector<UNI10_MATRIX> outs;
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a diagonal matrix. Need not to do so.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    //GPU_NOT_READY
    outs.push_back(UNI10_MATRIX(Rnum, Cnum, true, ongpu));
    outs.push_back(UNI10_MATRIX(Rnum, Cnum, false, ongpu));
    eigSyDecompose(m_elem, Rnum, outs[0].m_elem, outs[1].m_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::eigh():");
  }
	return outs;
}

#ifndef UNI10_PURE_REAL
std::vector<CMatrix> UNI10_BLOCK::eig()const{
  std::vector<CMatrix> outs;
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a diagonal matrix. Need not to do so.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    //GPU_NOT_READY
    outs.push_back(CMatrix(Rnum, Cnum, true, ongpu));
    outs.push_back(CMatrix(Rnum, Cnum, false, ongpu));
    eigDecompose(m_elem, Rnum, outs[0].m_elem, outs[1].m_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::eigh():");
  }
	return outs;
}
#endif

std::vector<UNI10_MATRIX> UNI10_BLOCK::svd()const{
	std::vector<UNI10_MATRIX> outs;
  try{
	if(diag){
    std::ostringstream err;
    err<<"Cannot perform singular value decomposition on a diagonal matrix. Need not to do so.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	size_t min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
  //GPU_NOT_READY
	outs.push_back(UNI10_MATRIX(Rnum, min, false, ongpu));
  outs.push_back(UNI10_MATRIX(min, min, true, ongpu));
	outs.push_back(UNI10_MATRIX(min, Cnum, false, ongpu));
	matrixSVD(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem, outs[2].m_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::svd():");
  }
	return outs;
}
UNI10_MATRIX UNI10_BLOCK::inverse()const{
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform inversion on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UNI10_MATRIX invM(*this);
    assert(ongpu == invM.isOngpu());
    matrixInv(invM.m_elem, Rnum, invM.diag, invM.ongpu);
    return invM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::inverse():");
    return UNI10_MATRIX();
  }
}

size_t UNI10_BLOCK::lanczosEigh(double& E0, UNI10_MATRIX& psi, size_t max_iter, double err_tol)const{
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform Lanczos algorithm to find the lowest eigen value and eigen vector on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(Rnum == psi.elemNum())){
      std::ostringstream err;
      err<<"Error in Lanczos initial vector psi. The vector dimension does not match with the number of the columns.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(ongpu && !psi.ongpu){
      if(!psi.toGPU()){
        std::ostringstream err;
        err<<"Error when allocating GPU global memory.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    }
    size_t iter = max_iter;
    if(!lanczosEV(m_elem, psi.m_elem, Rnum, iter, err_tol, E0, psi.m_elem, ongpu)){
      std::ostringstream err;
      err<<"Lanczos algorithm fails in converging.";;
      throw std::runtime_error(exception_msg(err.str()));
    }
    return iter;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::lanczosEigh(double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
    return 0;
  }
}

double UNI10_BLOCK::norm()const{
  try{
	  return vectorNorm(m_elem, elemNum(), 1, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::norm():");
    return 0;
  }
}

UNI10_DTYPE UNI10_BLOCK::sum()const{
  try{
	  return vectorSum(m_elem, elemNum(), 1, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::sum():");
    return 0;
  }
}

UNI10_DTYPE UNI10_BLOCK::trace()const{
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform trace on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag)
      return vectorSum(m_elem, elemNum(), 1, ongpu);
    else
      return vectorSum(m_elem, Cnum, Cnum + 1, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::trace():");
    return 0;
  }
}
UNI10_MATRIX operator* (const UNI10_BLOCK& Ma, const UNI10_BLOCK& Mb){
  try{
    if(!(Ma.Cnum == Mb.Rnum)){
      std::ostringstream err;
      err<<"The dimensions of the two matrices do not match for matrix multiplication.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((!Ma.diag) && (!Mb.diag)){
      UNI10_MATRIX Mc(Ma.Rnum, Mb.Cnum);
      matrixMul(Ma.m_elem, Mb.m_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.m_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
      return Mc;
    }
    else if(Ma.diag && (!Mb.diag)){
      UNI10_MATRIX Mc(Mb);
      Mc.resize(Ma.Rnum, Mb.Cnum);
      diagRowMul(Mc.m_elem, Ma.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
      return Mc;
    }
    else if((!Ma.diag) && Mb.diag){
      UNI10_MATRIX Mc(Ma);
      Mc.resize(Ma.Rnum, Mb.Cnum);
      diagColMul(Mc.m_elem, Mb.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
      return Mc;
    }
    else{
      UNI10_MATRIX Mc(Ma.Rnum, Mb.Cnum, true);
      Mc.set_zero();
      size_t min = std::min(Ma.elemNum(), Mb.elemNum());
      elemCopy(Mc.m_elem, Ma.m_elem, min * sizeof(UNI10_DTYPE), Mc.ongpu, Ma.ongpu);
      vectorMul(Mc.m_elem, Mb.m_elem, min, Mc.ongpu, Mb.ongpu);
      return Mc;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::Matrix&, uni10::Matrix&):");
    return UNI10_BLOCK();
  }
}
UNI10_MATRIX operator*(const UNI10_BLOCK& Ma, double a){
  try{
    UNI10_MATRIX Mb(Ma);
    vectorScal(a, Mb.m_elem, Mb.elemNum(), Mb.isOngpu());
    return Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::Matrix&, double):");
    return UNI10_BLOCK();
  }
}
UNI10_MATRIX operator*(double a, const UNI10_BLOCK& Ma){return Ma * a;}

#ifndef UNI10_PURE_REAL
CMatrix operator*(const UNI10_BLOCK& Ma, const std::complex<double>& a){
  try{
    CMatrix Mb(Ma);
    vectorScal(a, Mb.getElem(), Mb.elemNum(), Mb.isOngpu());
    return Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::Matrix&, double):");
    return UNI10_BLOCK();
  }
}
CMatrix operator*(const std::complex<double>& a, const UNI10_BLOCK& Ma){return Ma * a;}
#endif

UNI10_MATRIX operator+(const UNI10_BLOCK& Ma, const UNI10_BLOCK& Mb){
  try{
    UNI10_MATRIX Mc(Ma);
    vectorAdd(Mc.m_elem, Mb.m_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
    return Mc;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator+(uni10::Matrix&, uni10::Matrix&):");
    return UNI10_BLOCK();
  }
}

bool operator== (const UNI10_BLOCK& m1, const UNI10_BLOCK& m2){
  try{
    double diff;
    if(m1.elemNum() == m2.elemNum()){
      for(size_t i = 0; i < m1.elemNum(); i++){
        diff = std::abs(m1.m_elem[i] - m2.m_elem[i]);
        if(diff > 1E-12)
          return false;
      }
    }
    else
      return false;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator==(uni10::Matrix&, uni10::Matrix&):");
  }
  return true;
}

UNI10_BLOCK::~UNI10_BLOCK(){}
std::ostream& operator<< (std::ostream& os, const UNI10_BLOCK& b){
  try{
    os << b.Rnum << " x " << b.Cnum << " = " << b.elemNum();
    if(b.diag)
      os << ", Diagonal";
    if(b.ongpu)
      os<< ", onGPU";
    os <<std::endl << std::endl;
    UNI10_DTYPE* elem;
    if(b.ongpu){
      elem = (UNI10_DTYPE*)malloc(b.elemNum() * sizeof(UNI10_DTYPE));
      elemCopy(elem, b.m_elem, b.elemNum() * sizeof(UNI10_DTYPE), false, b.ongpu);
    }
    else
      elem = b.m_elem;
    for(size_t i = 0; i < b.Rnum; i++){
      for(size_t j = 0; j < b.Cnum; j++)
        if(b.diag){
#ifdef UNI10_COMPLEX
          if(i == j)
            os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i];
          else
            os << std::setw(17) << std::fixed << std::setprecision(3) << 0.0;
        }
        else
          os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i * b.Cnum + j];
#else
          if(i == j)
            os << std::setw(7) << std::fixed << std::setprecision(3) << elem[i];
          else
            os << std::setw(7) << std::fixed << std::setprecision(3) << 0.0;
        }
        else
          os << std::setw(7) << std::fixed << std::setprecision(3) << elem[i * b.Cnum + j];
#endif
      os << std::endl << std::endl;
    }
    if(b.ongpu)
      free(elem);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator<<(std::ostream&, uni10::Matrix&):");
  }
  return os;
}

};	/* namespace uni10 */
#ifdef UNI10_BLOCK
#undef UNI10_BLOCK
#endif
#ifdef UNI10_MATRIX
#undef UNI10_MATRIX
#endif
#ifdef UNI10_DTYPE
#undef UNI10_DTYPE
#endif

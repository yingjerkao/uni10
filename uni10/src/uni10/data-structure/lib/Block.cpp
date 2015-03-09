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
#include <uni10/tensor-network/CMatrix.h>

typedef double Real;


namespace uni10{
Block::Block(): Rnum(0), Cnum(0), diag(false), ongpu(false), m_elem(NULL){}
Block::Block(size_t _Rnum, size_t _Cnum, bool _diag): Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL){}
Block::Block(const Block& _b): Rnum(_b.Rnum), Cnum(_b.Cnum), diag(_b.diag), ongpu(_b.ongpu), m_elem(_b.m_elem){}

size_t Block::row()const{return Rnum;}
size_t Block::col()const{return Cnum;}
bool Block::isDiag()const{return diag;}
bool Block::isOngpu()const{return ongpu;}
size_t Block::elemNum()const{
  if(diag)
    return (Rnum < Cnum ? Rnum : Cnum);
  else
    return Rnum * Cnum;
}

Real* Block::getElem()const{return m_elem;}

Matrix Block::getDiag()const{
  try{
    if(diag)
      return *this;
    else{
      Matrix D(Rnum, Cnum, true, ongpu);
      ::uni10::getDiag(m_elem, D.getElem(), Rnum, Cnum, D.elemNum(), ongpu, D.isOngpu());
      return D;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::getDiag():");
    return Matrix();
  }
}

void Block::save(const std::string& fname)const{
  try{
    FILE *fp = fopen(fname.c_str(), "w");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in writing to file '"<<fname<<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    Real* elem = m_elem;
    if(ongpu){
      elem = (Real*)malloc(elemNum() * sizeof(Real));
      elemCopy(elem, m_elem, elemNum() * sizeof(Real), false, ongpu);
    }
    fwrite(elem, sizeof(Real), elemNum(), fp);
    fclose(fp);
    if(ongpu)
      free(elem);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::save(std::string&):");
  }
}

Real Block::operator[](size_t idx)const{
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

Real Block::at(size_t r, size_t c)const{
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

std::vector<Matrix> Block::eigh()const{
  std::vector<Matrix> outs;
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
    outs.push_back(Matrix(Rnum, Cnum, true, ongpu));
    outs.push_back(Matrix(Rnum, Cnum, false, ongpu));
    Matrix Eig(Rnum, Cnum, true, ongpu);
    eigSyDecompose(m_elem, Rnum, Eig.m_elem, outs[1].m_elem, ongpu);
    outs[0] = Eig;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::eigh():");
  }
	return outs;
}

std::vector<CMatrix> Block::eig()const{
  std::vector<CMatrix> outs;
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a diagonal matrix. Nothing to do.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    //GPU_NOT_READY
    outs.push_back(CMatrix(Rnum, Cnum, true, ongpu));
    outs.push_back(CMatrix(Rnum, Cnum, false, ongpu));
    eigDecompose(m_elem, Rnum, outs[0].m_elem, outs[1].m_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::eig():");
  }
	return outs;
}


std::vector<Matrix> Block::svd()const{
	std::vector<Matrix> outs;
  try{
	if(diag){
    std::ostringstream err;
    err<<"Cannot perform singular value decomposition on a diagonal matrix. Nothing to do.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	size_t min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
  //GPU_NOT_READY
	outs.push_back(Matrix(Rnum, min, false, ongpu));
  outs.push_back(Matrix(min, min, true, ongpu));
	outs.push_back(Matrix(min, Cnum, false, ongpu));
	matrixSVD(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem, outs[2].m_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::svd():");
  }
	return outs;
}
Matrix Block::inverse()const{
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform inversion on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    Matrix invM(*this);
    assert(ongpu == invM.isOngpu());
    matrixInv(invM.m_elem, Rnum, invM.diag, invM.ongpu);
    return invM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::inverse():");
    return Matrix();
  }
}

size_t Block::lanczosEigh(double& E0, Matrix& psi, size_t max_iter, double err_tol)const{
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

double Block::norm()const{
  try{
	  return vectorNorm(m_elem, elemNum(), 1, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::norm():");
    return 0;
  }
}

Real Block::sum()const{
  try{
	  return vectorSum(m_elem, elemNum(), 1, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::sum():");
    return 0;
  }
}

Real Block::trace()const{
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
Matrix operator* (const Block& Ma, const Block& Mb){
  try{
    if(!(Ma.Cnum == Mb.Rnum)){
      std::ostringstream err;
      err<<"The dimensions of the two matrices do not match for matrix multiplication.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((!Ma.diag) && (!Mb.diag)){
      Matrix Mc(Ma.Rnum, Mb.Cnum);
      matrixMul(Ma.m_elem, Mb.m_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.m_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
      return Mc;
    }
    else if(Ma.diag && (!Mb.diag)){
      Matrix Mc(Mb);
      Mc.resize(Ma.Rnum, Mb.Cnum);
      diagRowMul(Mc.m_elem, Ma.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
      return Mc;
    }
    else if((!Ma.diag) && Mb.diag){
      Matrix Mc(Ma);
      Mc.resize(Ma.Rnum, Mb.Cnum);
      diagColMul(Mc.m_elem, Mb.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
      return Mc;
    }
    else{
      Matrix Mc(Ma.Rnum, Mb.Cnum, true);
      Mc.set_zero();
      size_t min = std::min(Ma.elemNum(), Mb.elemNum());
      elemCopy(Mc.m_elem, Ma.m_elem, min * sizeof(Real), Mc.ongpu, Ma.ongpu);
      vectorMul(Mc.m_elem, Mb.m_elem, min, Mc.ongpu, Mb.ongpu);
      return Mc;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::Matrix&, uni10::Matrix&):");
    return Block();
  }
}
Matrix operator*(const Block& Ma, double a){
  try{
    Matrix Mb(Ma);
    vectorScal(a, Mb.m_elem, Mb.elemNum(), Mb.isOngpu());
    return Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::Matrix&, double):");
    return Block();
  }
}
Matrix operator*(double a, const Block& Ma){return Ma * a;}

#ifndef UNI10_PURE_REAL
CMatrix operator*(const Block& Ma, const std::complex<double>& a){
  try{
    CMatrix Mb(Ma);
    vectorScal(a, Mb.getElem(), Mb.elemNum(), Mb.isOngpu());
    return Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::Matrix&, double):");
    return Block();
  }
}
CMatrix operator*(const std::complex<double>& a, const Block& Ma){return Ma * a;}
#endif

Matrix operator+(const Block& Ma, const Block& Mb){
  try{
    Matrix Mc(Ma);
    vectorAdd(Mc.m_elem, Mb.m_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
    return Mc;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator+(uni10::Matrix&, uni10::Matrix&):");
    return Block();
  }
}

bool operator== (const Block& m1, const Block& m2){
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

Block::~Block(){}
std::ostream& operator<< (std::ostream& os, const Block& b){
  try{
    os << b.Rnum << " x " << b.Cnum << " = " << b.elemNum();
    if(b.diag)
      os << ", Diagonal";
    if(b.ongpu)
      os<< ", onGPU";
    os <<std::endl << std::endl;
    Real* elem;
    if(b.ongpu){
      elem = (Real*)malloc(b.elemNum() * sizeof(Real));
      elemCopy(elem, b.m_elem, b.elemNum() * sizeof(Real), false, b.ongpu);
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
#ifdef Block
#undef Block
#endif
#ifdef Matrix
#undef Matrix
#endif
#ifdef Real
#undef Real
#endif

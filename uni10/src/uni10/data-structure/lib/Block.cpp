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
typedef std::complex<double> Complex;

namespace uni10{
  
  std::ostream& operator<< (std::ostream& os, const matrixType& tp){
    try{
      if(tp == REAL)
        os << "The matrix type is REAL.";
      if(tp == COMPLEX)
        os << "The matrix type is COMPLEX.";
      if(tp == EMPTY)
        os << "This matrix is EMPTY.";
      os << std::endl << std::endl;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator<<(std::ostream&, uni10::matrixType&):");
    }
    return os;
  }

  Block::Block():m_type(EMPTY), Rnum(0), Cnum(0), diag(false), ongpu(false), cm_elem(NULL){}

  Block::Block(size_t _Rnum, size_t _Cnum, bool _diag): m_type(EMPTY), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL), cm_elem(NULL){}

  Block::Block(matrixType _tp, size_t _Rnum, size_t _Cnum, bool _diag): m_type(_tp), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL), cm_elem(NULL){}

  Block::Block(const Block& _b): m_type(_b.m_type), Rnum(_b.Rnum), Cnum(_b.Cnum), diag(_b.diag), ongpu(_b.ongpu), m_elem(_b.m_elem), cm_elem(_b.cm_elem){}
  
  Block::~Block(){}

  Real* Block::getElem()const{return m_elem;}
  
  Real* Block::getRealElem()const{return m_elem;}

  Complex* Block::getComplexElem()const{return cm_elem;}

  void Block::RtoC(){
    try{
      if(m_type == EMPTY){
        std::ostringstream err;
        err<<"This matrix is EMPTY. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(m_type == REAL){
        m_type = COMPLEX;
        cm_elem = (Complex*)elemAlloc(elemNum() * sizeof(Complex), ongpu);
        elemCast(cm_elem, m_elem, elemNum(), ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
        m_elem = NULL;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::RtoC():");
    }
  }

  matrixType Block::getType()const{
    return m_type;
  }

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
  
  std::vector<Matrix> Block::qr()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QR decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      matrixType tp = ( m_type == REAL ) ? REAL : COMPLEX;
      outs.push_back(Matrix(tp, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(tp, Cnum, Cnum, false, ongpu));
      if(!diag){
        if(m_type == REAL)
          matrixQR(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem);
        if(m_type == COMPLEX)
          matrixQR(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        if(m_type == REAL){
          for(int i = 0; i < min; i++)
            tmpR[i*min+i] = m_elem[i];
          matrixQR(tmpR, min, min, outs[0].m_elem, outs[1].m_elem);
        }
        if(m_type == COMPLEX){ 
          for(int i = 0; i < min; i++)
            tmpC[i*min+i] = cm_elem[i];
          matrixQR(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem);
        }
        free(tmpC);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::qr():");
    }
    return outs;
  }

  std::vector<Matrix> Block::rq()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      matrixType tp = ( m_type == REAL ) ? REAL : COMPLEX;
      outs.push_back(Matrix(tp, Rnum, Rnum, false, ongpu)); //r
      outs.push_back(Matrix(tp, Rnum, Cnum, false, ongpu)); //q
      if(!diag){
        if(m_type == REAL)
          matrixRQ(m_elem, Rnum, Cnum, outs[1].m_elem, outs[0].m_elem);
        if(m_type == COMPLEX){
          matrixRQ(cm_elem, Rnum, Cnum, outs[1].cm_elem, outs[0].cm_elem);
        }
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        if(m_type == REAL){
          for(int i = 0; i < min; i++)
            tmpR[i*min+i] = m_elem[i];
          matrixRQ(tmpR, min, min, outs[1].m_elem, outs[0].m_elem);
        }
        if(m_type == COMPLEX){ 
          for(int i = 0; i < min; i++)
            tmpC[i*min+i] = cm_elem[i];
          matrixRQ(tmpC, min, min, outs[1].cm_elem, outs[0].cm_elem);
        }
        free(tmpC);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::rq():");
    }
    return outs;
  }

  std::vector<Matrix> Block::ql()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QL decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      matrixType tp = ( m_type == REAL ) ? REAL : COMPLEX;
      outs.push_back(Matrix(tp, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(tp, Cnum, Cnum, false, ongpu));
      if(!diag){
        if(m_type == REAL)
          matrixQL(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem);
        if(m_type == COMPLEX)
          matrixQL(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        if(m_type == REAL){
          for(int i = 0; i < min; i++)
            tmpR[i*min+i] = m_elem[i];
          matrixQL(tmpR, min, min, outs[0].m_elem, outs[1].m_elem);
        }
        if(m_type == COMPLEX){ 
          for(int i = 0; i < min; i++)
            tmpC[i*min+i] = cm_elem[i];
          matrixQL(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem);
        }
        free(tmpC);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::ql():");
    }
    return outs;
  }

  std::vector<Matrix> Block::lq()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      } 
      matrixType tp = ( m_type == REAL ) ? REAL : COMPLEX;
      outs.push_back(Matrix(tp, Rnum, Rnum, false, ongpu));
      outs.push_back(Matrix(tp, Rnum, Cnum, false, ongpu));
      if(!diag){
        if(m_type == REAL)
          matrixLQ(m_elem, Rnum, Cnum, outs[1].m_elem, outs[0].m_elem);
        if(m_type == COMPLEX)
          matrixLQ(cm_elem, Rnum, Cnum, outs[1].cm_elem, outs[0].cm_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        if(m_type == REAL){
          for(int i = 0; i < min; i++)
            tmpR[i*min+i] = m_elem[i];
          matrixLQ(tmpR, min, min, outs[1].m_elem, outs[0].m_elem);
        }
        if(m_type == COMPLEX){ 
          for(int i = 0; i < min; i++)
            tmpC[i*min+i] = cm_elem[i];
          matrixLQ(tmpC, min, min, outs[1].cm_elem, outs[0].cm_elem);
        }
        free(tmpC);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::lq():");
    }
    return outs;
  }

  std::vector<Matrix> Block::svd()const{
    std::vector<Matrix> outs;
    try{
    /*  
      if(diag){
        std::ostringstream err;
        err<<"Cannot perform singular value decomposition on a diagonal matrix. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    */
      size_t min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
      //GPU_NOT_READY
      matrixType tp = ( m_type == REAL ) ? REAL : COMPLEX;
      outs.push_back(Matrix(tp, Rnum, min, false, ongpu));
      outs.push_back(Matrix(tp, min, min, true, ongpu));
      outs.push_back(Matrix(tp, min, Cnum, false, ongpu));
      if(!diag){
        if(m_type == REAL)
          matrixSVD(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem, outs[2].m_elem, ongpu);
        if(m_type == COMPLEX)
          matrixSVD(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem, outs[2].cm_elem, ongpu);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        if(m_type == REAL){
          for(int i = 0; i < min; i++)
            tmpR[i*min+i] = m_elem[i];
          matrixSVD(tmpR, min, min, outs[0].m_elem, outs[1].m_elem, outs[2].m_elem, ongpu);
        }
        if(m_type == COMPLEX){ 
          for(int i = 0; i < min; i++)
            tmpC[i*min+i] = cm_elem[i];
          matrixSVD(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem, outs[2].cm_elem, ongpu);
        }
        free(tmpC);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::svd():");
    }
    return outs;
  }
  
  Matrix Block::getDiag()const{
    try{
      if(diag)
        return *this;
      else{
        matrixType tp = ( m_type == REAL ) ? REAL : COMPLEX;
        Matrix D(tp, Rnum, Cnum, true, ongpu);
        if(m_type == REAL)
          ::uni10::getDiag(m_elem, D.getRealElem(), Rnum, Cnum, D.elemNum(), ongpu, D.isOngpu());
        if(m_type == COMPLEX)
          ::uni10::getDiag(cm_elem, D.getComplexElem(), Rnum, Cnum, D.elemNum(), ongpu, D.isOngpu());
        return D;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::getDiag():");
      return Matrix();
    }
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
      if(m_type == REAL)
        matrixInv(invM.m_elem, Rnum, invM.diag, invM.ongpu);
      if(m_type == COMPLEX)
        matrixInv(invM.cm_elem, Rnum, invM.diag, invM.ongpu);
      return invM;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::inverse():");
      return Matrix();
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
      matrixType tp = ( m_type == REAL ) ? REAL : COMPLEX;
      outs.push_back(Matrix(tp, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(tp, Rnum, Cnum, false, ongpu));
      Matrix Eig(REAL, Rnum, Cnum, true, ongpu);
      if(m_type == REAL)
        eigSyDecompose(m_elem, Rnum, Eig.m_elem, outs[1].m_elem, ongpu);
      if(m_type == COMPLEX)
        eigSyDecompose(cm_elem, Rnum, Eig.m_elem, outs[1].cm_elem, ongpu);
      outs[0] = Eig;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eigh():");
    }
    return outs;
  }
  
  void Block::save(const std::string& fname)const{
    try{
      FILE *fp = fopen(fname.c_str(), "w");
      if(!(fp != NULL)){
        std::ostringstream err;
        err<<"Error in writing to file '"<<fname<<"'.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      fwrite(&m_type, sizeof(m_type), 1, fp);
      fwrite(&Rnum, sizeof(Rnum), 1, fp);
      fwrite(&Cnum, sizeof(Cnum), 1, fp);
      fwrite(&diag, sizeof(diag), 1, fp);
      fwrite(&ongpu, sizeof(ongpu), 1, fp);
      if(m_type == REAL){
        Real* elem = m_elem;
        if(ongpu){
          elem = (Real*)malloc(elemNum() * sizeof(Real));
          elemCopy(elem, m_elem, elemNum() * sizeof(Real), false, ongpu);
        }
        fwrite(elem, sizeof(Real), elemNum(), fp);
        if(ongpu)
          free(elem);
      }
      if(m_type == COMPLEX){
        Complex* elem = cm_elem;
        if(ongpu){
          elem = (Complex*)malloc(elemNum() * sizeof(Complex));
          elemCopy(elem, cm_elem, elemNum() * sizeof(Complex), false, ongpu);
        }
        fwrite(elem, sizeof(Complex), elemNum(), fp);
        if(ongpu)
          free(elem);
      }
      fclose(fp);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::save(std::string&):");
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
      if(m_type == REAL){
        if(!lanczosEV(m_elem, psi.m_elem, Rnum, iter, err_tol, E0, psi.m_elem, ongpu)){
          std::ostringstream err;
          err<<"Lanczos algorithm fails in converging.";;
          throw std::runtime_error(exception_msg(err.str()));
        }
      }
      if(m_type == COMPLEX){
        if(!lanczosEV(cm_elem, psi.cm_elem, Rnum, iter, err_tol, E0, psi.cm_elem, ongpu)){
          std::ostringstream err;
          err<<"Lanczos algorithm fails in converging.";;
          throw std::runtime_error(exception_msg(err.str()));
        }
      }
      return iter;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::lanczosEigh(double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
      return 0;
    }
  }
  
  bool operator== (const Block& m1, const Block& m2){
    try{
      if(m1.m_type != m2.m_type)
        return false;
      else{
        double diff;
        if(m1.elemNum() == m2.elemNum()){
          for(size_t i = 0; i < m1.elemNum(); i++){
            diff = (m1.m_type == REAL) ? std::abs(m1.m_elem[i] - m2.m_elem[i]) : std::abs(m1.cm_elem[i] -m2.cm_elem[i]);
            if(diff > 1E-12)
              return false;
          }
        }
        else
          return false;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator==(uni10::Matrix&, uni10::Matrix&):");
    }
    return true;
  }

  Matrix operator* (const Block& Ma, const Block& Mb){
    try{
      if(!(Ma.Cnum == Mb.Rnum)){
        std::ostringstream err;
        err<<"The dimensions of the two matrices do not match for matrix multiplication.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      // R*R
      if(Ma.m_type == REAL && Mb.m_type == REAL){
        if((!Ma.diag) && (!Mb.diag)){
          Matrix Mc(REAL, Ma.Rnum, Mb.Cnum);
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
          Matrix Mc(REAL, Ma.Rnum, Mb.Cnum, true);
          Mc.set_zero();
          size_t min = std::min(Ma.elemNum(), Mb.elemNum());
          elemCopy(Mc.m_elem, Ma.m_elem, min * sizeof(Real), Mc.ongpu, Ma.ongpu);
          vectorMul(Mc.m_elem, Mb.m_elem, min, Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      if(Ma.m_type == COMPLEX && Mb.m_type == COMPLEX){
        if((!Ma.diag) && (!Mb.diag)){
          Matrix Mc(COMPLEX, Ma.Rnum, Mb.Cnum);
          matrixMul(Ma.cm_elem, Mb.cm_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.cm_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
          return Mc;
        }
        else if(Ma.diag && (!Mb.diag)){
          Matrix Mc(Mb);
          Mc.resize(Ma.Rnum, Mb.Cnum);
          diagRowMul(Mc.cm_elem, Ma.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else if((!Ma.diag) && Mb.diag){
          Matrix Mc(Ma);
          Mc.resize(Ma.Rnum, Mb.Cnum);
          diagColMul(Mc.cm_elem, Mb.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else{
          Matrix Mc(COMPLEX, Ma.Rnum, Mb.Cnum, true);
          Mc.set_zero();
          size_t min = std::min(Ma.elemNum(), Mb.elemNum());
          elemCopy(Mc.cm_elem, Ma.cm_elem, min * sizeof(Complex), Mc.ongpu, Ma.ongpu);
          vectorMul(Mc.cm_elem, Mb.cm_elem, min, Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      
      if(Ma.m_type == REAL && Mb.m_type == COMPLEX){
        Matrix _Ma(Ma);
        _Ma.RtoC();
        if((!_Ma.diag) && (!Mb.diag)){
          Matrix Mc(COMPLEX, _Ma.Rnum, Mb.Cnum);
          matrixMul(_Ma.cm_elem, Mb.cm_elem, _Ma.Rnum, Mb.Cnum, _Ma.Cnum, Mc.cm_elem, _Ma.ongpu, Mb.ongpu, Mc.ongpu);
          return Mc;
        }
        else if(_Ma.diag && (!Mb.diag)){
          Matrix Mc(Mb);
          Mc.resize(_Ma.Rnum, Mb.Cnum);
          diagRowMul(Mc.cm_elem, _Ma.cm_elem, Mc.Rnum, Mc.Cnum, _Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else if((!_Ma.diag) && Mb.diag){
          Matrix Mc(_Ma);
          Mc.resize(_Ma.Rnum, Mb.Cnum);
          diagColMul(Mc.cm_elem, Mb.cm_elem, Mc.Rnum, Mc.Cnum, _Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else{
          Matrix Mc(COMPLEX, _Ma.Rnum, Mb.Cnum, true);
          Mc.set_zero();
          size_t min = std::min(_Ma.elemNum(), Mb.elemNum());
          elemCopy(Mc.cm_elem, _Ma.cm_elem, min * sizeof(Complex), Mc.ongpu, _Ma.ongpu);
          vectorMul(Mc.cm_elem, Mb.cm_elem, min, Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      if(Ma.m_type == COMPLEX && Ma.m_type == REAL){
        Matrix _Mb(Mb);
        _Mb.RtoC();
        if((!Ma.diag) && (!_Mb.diag)){
          Matrix Mc(COMPLEX, Ma.Rnum, _Mb.Cnum);
          matrixMul(Ma.cm_elem, _Mb.cm_elem, Ma.Rnum, _Mb.Cnum, Ma.Cnum, Mc.cm_elem, Ma.ongpu, _Mb.ongpu, Mc.ongpu);
          return Mc;
        }
        else if(Ma.diag && (!_Mb.diag)){
          Matrix Mc(_Mb);
          Mc.resize(Ma.Rnum, _Mb.Cnum);
          diagRowMul(Mc.cm_elem, Ma.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else if((!Ma.diag) && _Mb.diag){
          Matrix Mc(Ma);
          Mc.resize(Ma.Rnum, _Mb.Cnum);
          diagColMul(Mc.cm_elem, _Mb.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else{
          Matrix Mc(COMPLEX, Ma.Rnum, _Mb.Cnum, true);
          Mc.set_zero();
          size_t min = std::min(Ma.elemNum(), _Mb.elemNum());
          elemCopy(Mc.cm_elem, Ma.cm_elem, min * sizeof(Complex), Mc.ongpu, Ma.ongpu);
          vectorMul(Mc.cm_elem, _Mb.cm_elem, min, Mc.ongpu, _Mb.ongpu);
          return Mc;
        }
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
      if(Ma.m_type == REAL)
        vectorScal(a, Mb.m_elem, Mb.elemNum(), Mb.isOngpu());
      if(Ma.m_type == COMPLEX)
        vectorScal(a, Mb.cm_elem, Mb.elemNum(), Mb.isOngpu());
      return Mb;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::Matrix&, double):");
      return Block();
    }
  }
  Matrix operator*(double a, const Block& Ma){return Ma * a;}

  Matrix operator*(const Block& Ma, const std::complex<double>& a){
    try{
      Matrix Mb(Ma);
      if(a.imag() == 0){
        double _a = a.real();
        if(Ma.m_type == REAL) 
          vectorScal(_a, Mb.m_elem, Mb.elemNum(), Mb.isOngpu());
        if(Ma.m_type == COMPLEX)
          vectorScal(_a, Mb.cm_elem, Mb.elemNum(), Mb.isOngpu());
      }
      else{
        Mb.RtoC();
        vectorScal(a, Mb.cm_elem, Mb.elemNum(), Mb.isOngpu());
      }
      return Mb;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::Matrix&, Complex):");
      return Block();
    }
  }
  Matrix operator*(const std::complex<double>& a, const Block& Ma){return Ma * a;}
  
  Matrix operator+(const Block& Ma, const Block& Mb){
    /*if ((Ma.diag && !Mb.diag) || (!Ma.diag && Mb.diag) ){
      std::ostringstream err;
      err<<"Fatal error(code = T1). Addition of a diagonal and a regular matrices is not implemented.";
      throw std::runtime_error(exception_msg(err.str())); ;
      };*/

    try{
      if(Ma.m_type == REAL && Mb.m_type == REAL){
        if (Ma.diag && !Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(REAL, Ma.Rnum,Ma.Cnum);
          setDiag(Mc.m_elem,Ma.m_elem,Mc.Rnum,Mc.Cnum,Ma.Rnum,Mc.ongpu,Ma.ongpu);
          vectorAdd(Mc.m_elem, Mb.m_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        } else if (Mb.diag && !Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;

          Matrix Mc(REAL, Mb.Rnum, Mb.Cnum);
          setDiag(Mc.m_elem,Mb.m_elem,Mc.Rnum,Mc.Cnum,Mb.Cnum,Mc.ongpu,Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.m_elem, Ma.m_elem, Mc.elemNum(), Mc.ongpu, Ma.ongpu);

          return Mc;
        } else {
          Matrix Mc(Ma);
          vectorAdd(Mc.m_elem, Mb.m_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      if(Ma.m_type == COMPLEX && Mb.m_type == COMPLEX){
        if (Ma.diag && !Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(COMPLEX, Ma.Rnum,Ma.Cnum);
          setDiag(Mc.cm_elem,Ma.cm_elem,Mc.Rnum,Mc.Cnum,Ma.Rnum,Mc.ongpu,Ma.ongpu);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        } else if (Mb.diag && !Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;

          Matrix Mc(COMPLEX, Mb.Rnum,Mb.Cnum);
          setDiag(Mc.cm_elem,Mb.cm_elem,Mc.Rnum,Mc.Cnum,Mb.Cnum,Mc.ongpu,Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.cm_elem, Ma.cm_elem, Mc.elemNum(), Mc.ongpu, Ma.ongpu);

          return Mc;
        } else {
          Matrix Mc(Ma);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      if(Ma.m_type == REAL && Mb.m_type == COMPLEX){
        Matrix _Ma(Ma);
        _Ma.RtoC();
        if (_Ma.diag && !Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(COMPLEX, _Ma.Rnum,_Ma.Cnum);
          setDiag(Mc.cm_elem,_Ma.cm_elem,Mc.Rnum,Mc.Cnum,_Ma.Rnum,Mc.ongpu,_Ma.ongpu);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        } else if (Mb.diag && !_Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;

          Matrix Mc(COMPLEX, Mb.Rnum,Mb.Cnum);
          setDiag(Mc.cm_elem,Mb.cm_elem,Mc.Rnum,Mc.Cnum,Mb.Cnum,Mc.ongpu,Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.cm_elem, _Ma.cm_elem, Mc.elemNum(), Mc.ongpu, _Ma.ongpu);

          return Mc;
        } else {
          Matrix Mc(_Ma);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      if(Ma.m_type == COMPLEX && Mb.m_type == REAL){
        Matrix _Mb(Mb);
        _Mb.RtoC();
        if (Ma.diag && !_Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(COMPLEX, Ma.Rnum,Ma.Cnum);
          setDiag(Mc.cm_elem,Ma.cm_elem,Mc.Rnum,Mc.Cnum,Ma.Rnum,Mc.ongpu,Ma.ongpu);
          vectorAdd(Mc.cm_elem, _Mb.cm_elem, Mc.elemNum(), Mc.ongpu, _Mb.ongpu);
          return Mc;
        } else if (_Mb.diag && !Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;

          Matrix Mc(COMPLEX, _Mb.Rnum,_Mb.Cnum);
          setDiag(Mc.cm_elem,_Mb.cm_elem,Mc.Rnum,Mc.Cnum,_Mb.Cnum,Mc.ongpu,_Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.cm_elem, Ma.cm_elem, Mc.elemNum(), Mc.ongpu, Ma.ongpu);

          return Mc;
        } else {
          Matrix Mc(Ma);
          vectorAdd(Mc.cm_elem, _Mb.cm_elem, Mc.elemNum(), Mc.ongpu, _Mb.ongpu);
          return Mc;
        }
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator+(uni10::Matrix&, uni10::Matrix&):");
      return Block();
    }
  }
  
  std::ostream& operator<< (std::ostream& os, const Block& b){
    try{
      os << b.Rnum << " x " << b.Cnum << " = " << b.elemNum();
      if(b.m_type == REAL)
        os << ", REAL";
      if(b.m_type == COMPLEX)
        os << ", COMPLEX";
      if(b.diag)
        os << ", Diagonal";
      if(b.ongpu)
        os<< ", onGPU";
      os <<std::endl << std::endl;
      if(b.m_type == REAL){
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
              //#ifdef UNI10_COMPLEX
              if(i == j)
                os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i];
              else
                os << std::setw(17) << std::fixed << std::setprecision(3) << 0.0;
            }
            else
              os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i * b.Cnum + j];
          os << std::endl << std::endl;
        }
        if(b.ongpu)
          free(elem);
      }
      if(b.m_type == COMPLEX){
        Complex* elem;
        if(b.ongpu){
          elem = (Complex*)malloc(b.elemNum() * sizeof(Complex));
          elemCopy(elem, b.cm_elem, b.elemNum() * sizeof(Complex), false, b.ongpu);
        }
        else
          elem = b.cm_elem;
        for(size_t i = 0; i < b.Rnum; i++){
          for(size_t j = 0; j < b.Cnum; j++)
            if(b.diag){
              if(i == j)
                os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i];
              else
                os << std::setw(17) << std::fixed << std::setprecision(3) << 0.0;
            }
            else
              os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i * b.Cnum + j];
          os << std::endl << std::endl;
        }
        if(b.ongpu)
          free(elem);
      }
      if(b.m_type == EMPTY)
        os << "This matrix is EMPTY" <<std::endl << std::endl;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator<<(std::ostream&, uni10::Matrix&):");
    }
    return os;
  }
  
  
  Real Block::norm()const{
    try{
      if(m_type == EMPTY){
        std::ostringstream err;
        err<<"This matirx is empty" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(m_type == REAL)
        return vectorNorm(m_elem, elemNum(), 1, ongpu);
      if(m_type == COMPLEX)
        return vectorNorm(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::norm():");
      return 0;
    }
  }
  
  Complex Block::trace()const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform trace on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(m_type == REAL){
        if(diag)
          return Complex(vectorSum(m_elem, elemNum(), 1, ongpu), 0);
        else
          return Complex(vectorSum(m_elem, Cnum, Cnum + 1, ongpu), 0);
      }
      if(m_type == COMPLEX){
        if(diag)
          return vectorSum(cm_elem, elemNum(), 1, ongpu);
        else
          return vectorSum(cm_elem, Cnum, Cnum + 1, ongpu);
      }
    }catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::trace():");
      return 0;
    }
  }

  Complex Block::sum()const{
    try{
      if(m_type == REAL)
        return Complex(vectorSum(m_elem, elemNum(), 1, ongpu), 0);
      if(m_type == COMPLEX)
        return vectorSum(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::sum():");
      return 0;
    }
  }
  
  Complex Block::operator[](size_t idx)const{
    try{
      if(!(idx < elemNum())){
        std::ostringstream err;
        err<<"Index exceeds the number of elements("<<elemNum()<<").";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(m_type == REAL)
        return Complex(getElemAt(idx, m_elem, ongpu), 0);
      if(m_type == COMPLEX)
        return getElemAt(idx, cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::operator[](size_t):");
      return 0;
    }
  }

  Complex Block::at(size_t r, size_t c)const{
    try{
      if(!((r < Rnum) && (c < Cnum))){
        std::ostringstream err;
        err<<"The input indices are out of range.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(m_type == REAL){
        if(diag){
          if(!(r == c && r < elemNum())){
            std::ostringstream err;
            err<<"The matrix is diagonal, there is no off-diagonal element.";
            throw std::runtime_error(exception_msg(err.str()));
          }
          return Complex(getElemAt(r, m_elem, ongpu), 0);
        }
        else
          return Complex(getElemAt(r * Cnum + c, m_elem, ongpu), 0);
      }
      if(m_type == COMPLEX){
        if(diag){
          if(!(r == c && r < elemNum())){
            std::ostringstream err;
            err<<"The matrix is diagonal, there is no off-diagonal element.";
            throw std::runtime_error(exception_msg(err.str()));
          }
          return Complex(getElemAt(r, m_elem, ongpu), 0);
        }
        else
          return Complex(getElemAt(r * Cnum + c, m_elem, ongpu), 0);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::at(size_t, size_t):");
      return 0;
    }
  }

/************************************************************************/
/************************/

  std::vector<Matrix> Block::eig()const{
    std::vector<Matrix> outs;
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
      outs.push_back(Matrix(COMPLEX, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(COMPLEX, Rnum, Cnum, false, ongpu));
      if(m_type == REAL)
        eigDecompose(m_elem, Rnum, outs[0].cm_elem, outs[1].cm_elem, ongpu);
      if(m_type == COMPLEX)
        eigDecompose(cm_elem, Rnum, outs[0].cm_elem, outs[1].cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eig():");
    }
    return outs;
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

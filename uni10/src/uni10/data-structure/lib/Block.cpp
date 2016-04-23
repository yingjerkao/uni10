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
#include <uni10/numeric/lapack/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#include <uni10/tensor-network/Matrix.h>



namespace uni10{

  std::ostream& operator<< (std::ostream& os, const Block& b){
    try{
      os << b.Rnum << " x " << b.Cnum << " = " << b.elemNum();
      if(b.typeID() == 0){
        os <<std::endl << std::endl;
        os << "This matrix is EMPTY" <<std::endl << std::endl;
      }
      if(b.typeID() == 1){
        os << ", REAL";
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
      if(b.typeID() == 2){
        os << ", COMPLEX";
        if(b.diag)
          os << ", Diagonal";
        if(b.ongpu)
          os<< ", onGPU";
        os <<std::endl << std::endl;
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
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator<<(std::ostream&, uni10::Matrix&):");
    }
    return os;
  }

  bool operator== (const Block& m1, const Block& m2){
    try{
      if(m1.typeID() != m2.typeID())
        return false;
      else{
        double diff;
        if(m1.elemNum() == m2.elemNum()){
          for(size_t i = 0; i < m1.elemNum(); i++){
            diff = (m1.typeID() == 1) ? std::abs(m1.m_elem[i] - m2.m_elem[i]) : std::abs(m1.cm_elem[i] -m2.cm_elem[i]);
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
      if(Ma.typeID() == 1 && Mb.typeID() == 1)
        return RDotR(Ma, Mb);
      else if(Ma.typeID() == 2 && Mb.typeID() == 2)
        return CDotC(Ma, Mb);
      else if(Ma.typeID() == 1 && Mb.typeID() == 2)
        return RDotC(Ma, Mb);
      else if(Ma.typeID() == 2 && Mb.typeID() == 1)
        return CDotR(Ma, Mb);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::Matrix&, uni10::Matrix&):");
      return Block();
    }
    return Matrix();
  }

  Matrix operator*(const Block& Ma, double a){
    try{
      Matrix Mb(Ma);
      if(Ma.typeID() != 0)
        (Ma.typeID() == 1) ? vectorScal(a, Mb.m_elem, Mb.elemNum(), Mb.isOngpu()): vectorScal(a, Mb.cm_elem, Mb.elemNum(), Mb.isOngpu());
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
        return a.real() * Mb;
      }else{
        RtoC(Mb);
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
    try{
      Matrix Mc = Ma;
      Mc += Mb; 
      return Mc;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator+(uni10::Matrix&, uni10::Matrix&):");
      return Block();
    }
    return Block();
  }

  // Default Real
  Block::Block():r_flag(RTYPE), c_flag(CNULL), m_elem(NULL), cm_elem(NULL), Rnum(0), Cnum(0), diag(false), ongpu(false){}

  Block::Block(size_t _Rnum, size_t _Cnum, bool _diag): r_flag(RTYPE), c_flag(CNULL), m_elem(NULL), cm_elem(NULL), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false){}

  Block::Block(const Block& _b): r_flag(_b.r_flag), c_flag(_b.c_flag), m_elem(_b.m_elem), cm_elem(_b.cm_elem), Rnum(_b.Rnum), Cnum(_b.Cnum), diag(_b.diag), ongpu(_b.ongpu){}

  Block::~Block(){}

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

  int Block::typeID()const{
    return r_flag + c_flag;
  }

  void Block::savePrototype(const std::string& fname)const{
    FILE *fp = fopen(fname.c_str(), "w");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in writing to file '"<<fname<<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 0){
      fwrite(&r_flag, sizeof(r_flag), 1, fp);
      fwrite(&c_flag, sizeof(c_flag), 1, fp);
      fwrite(&Rnum, sizeof(Rnum), 1, fp);
      fwrite(&Cnum, sizeof(Cnum), 1, fp);
      fwrite(&diag, sizeof(diag), 1, fp);
      fwrite(&ongpu, sizeof(ongpu), 1, fp);
    }
    fclose(fp);
  }

  void Block::save(const std::string& fname)const{
    try{
      if(typeID() == 1)
        save(RTYPE, fname);
      else if(typeID() == 2)
        save(CTYPE, fname);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::save(std::string&):");
    }
  }

  std::vector<Matrix> Block::qr()const{
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QR decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return qr(RTYPE);
      else if(typeID() == 2)
        return qr(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::qr():");
    }
    return std::vector<Matrix>();
  }

  std::vector<Matrix> Block::rq()const{
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return rq(RTYPE);
      else if(typeID() == 2)
        return rq(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::rq():");
    }
    return std::vector<Matrix>();
  }

  std::vector<Matrix> Block::ql()const{
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QL decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return ql(RTYPE);
      else if(typeID() == 2)
        return ql(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::ql():");
    }
    return std::vector<Matrix>();
  }

  std::vector<Matrix> Block::lq()const{
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return lq(RTYPE);
      else if(typeID() == 2)
        return lq(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::lq():");
    }
    return std::vector<Matrix>();
  }

  std::vector<Matrix> Block::svd()const{
    try{
      if(typeID() == 1)
        return svd(RTYPE);
      else if(typeID() == 2)
        return svd(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::svd():");
    }
    return std::vector<Matrix>();
  }

  Real Block::norm()const{
    try{
      if(typeID() == 1)
        return vectorNorm(m_elem, elemNum(), 1, ongpu);
      else if(typeID() == 2)
        return vectorNorm(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::norm():");
    }
    return 0;
  }

  Matrix Block::inverse()const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform inversion on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return inverse(RTYPE); 
      else if(typeID() == 2)
        return inverse(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::inverse():");
    }
    return Matrix();
  }

  Matrix Block::getDiag()const{
    try{
      int _typeID = this->typeID();
      if(_typeID == 1)
        return getDiag(RTYPE);
      else if(_typeID == 2)
        return getDiag(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::getDiag():");
      return Matrix();
    }
    return Matrix();
  }

  Real Block::trace()const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform trace on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return trace(RTYPE);
      else if(typeID() == 2){
        std::ostringstream err;
        err<<"This matrix is Complex. Please use trace(uni10::cflag) instead.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    }catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::trace():");
    }
    return 0;
  }

  Real Block::sum()const{
    try{
      if(typeID() == 1)
        return vectorSum(m_elem, elemNum(), 1, ongpu);
      else if(typeID() == 2){
        std::ostringstream err;
        err<<"This matrix is Complex. Please use sum(uni10::cflag) instead.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::sum():");
    }
    return 0.;
  }

  std::vector<Matrix> Block::eig()const{
    try{
      if(typeID() == 1)
        return eig(RTYPE);
      else if(typeID() == 2)
        return eig(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eig():");
    }
    return std::vector<Matrix>();
  }

  std::vector<Matrix> Block::eigh()const{
    try{
      if(typeID() == 1)
        return eigh(RTYPE);
      else if(typeID() == 2)
        return eigh(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eigh():");
    }
    return std::vector<Matrix>();
  }

  Real Block::at(size_t r, size_t c)const{
    try{
      if(typeID() == 1)
        return at(RTYPE, r, c);
      else if(typeID() == 2){
        std::ostringstream err;
        err<<"This matrix is Complex. Please use at(uni10::cflag, size_t, size_t) instead.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::at(size_t, size_t):");
    }
    return 0;
  }

  bool Block::CelemIsNULL()const{return cm_elem == NULL;}

  bool Block::RelemIsNULL()const{return m_elem == NULL;}

  size_t lanczosEigh(rflag tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol){
    try{
      throwTypeError(tp);
      if(!(ori_mat.Rnum == ori_mat.Cnum)){
        std::ostringstream err;
        err<<"Cannot perform Lanczos algorithm to find the lowest eigen value and eigen vector on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(!(ori_mat.Rnum == psi.elemNum())){
        std::ostringstream err;
        err<<"Error in Lanczos initial vector psi. The vector dimension does not match with the number of the columns.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(ori_mat.ongpu && !psi.ongpu){
        if(!psi.toGPU()){
          std::ostringstream err;
          err<<"Error when allocating GPU global memory.";
          throw std::runtime_error(exception_msg(err.str()));
        }
      }
      size_t iter = max_iter;
      if(!lanczosEV(ori_mat.m_elem, psi.m_elem, ori_mat.Rnum, iter, err_tol, E0, psi.m_elem, ori_mat.ongpu)){
        std::ostringstream err;
        err<<"Lanczos algorithm fails in converging.";;
        throw std::runtime_error(exception_msg(err.str()));
      }
      return iter;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::lanczosEigh(uni10::rflag, double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
    }
    return 0;
  }

  size_t lanczosEigh(cflag tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol){
    try{
      throwTypeError(tp);
      if(!(ori_mat.Rnum == ori_mat.Cnum)){
        std::ostringstream err;
        err<<"Cannot perform Lanczos algorithm to find the lowest eigen value and eigen vector on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(!(ori_mat.Rnum == psi.elemNum())){
        std::ostringstream err;
        err<<"Error in Lanczos initial vector psi. The vector dimension does not match with the number of the columns.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(ori_mat.ongpu && !psi.ongpu){
        if(!psi.toGPU()){
          std::ostringstream err;
          err<<"Error when allocating GPU global memory.";
          throw std::runtime_error(exception_msg(err.str()));
        }
      }
      size_t iter = max_iter;
      if(!lanczosEV(ori_mat.cm_elem, psi.cm_elem, ori_mat.Rnum, iter, err_tol, E0, psi.cm_elem, ori_mat.ongpu)){
        std::ostringstream err;
        err<<"Lanczos algorithm fails in converging.";;
        throw std::runtime_error(exception_msg(err.str()));
      }
      return iter;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::lanczosEigh(uni10::cflag, double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
      return 0;
    }
  }

  size_t lanczosEigh(Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol){
    try{
      if(ori_mat.typeID() == 1)
        return lanczosEigh(RTYPE, ori_mat, E0, psi, max_iter, err_tol);
      else if(ori_mat.typeID() == 2)
        return lanczosEigh(CTYPE, ori_mat, E0, psi, max_iter, err_tol);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::lanczosEigh(double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
    }
    return 0;
  }

};	/* namespace uni10 */
/*
#ifdef Block
#undef Block
#endif
#ifdef Matrix
#undef Matrix
#endif
#ifdef Real
#undef Real
#endif
*/

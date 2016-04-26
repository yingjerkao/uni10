/****************************************************************************
*  @file BlockComplex.cpp
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

  Block::Block(cflag tp, size_t _Rnum, size_t _Cnum, bool _diag): r_flag(RNULL), c_flag(tp), m_elem(NULL), cm_elem(NULL), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false){}

  Complex Block::operator()(size_t idx)const{
    try{
      if(!(idx < elemNum())){
        std::ostringstream err;
        err<<"Index exceeds the number of elements("<<elemNum()<<").";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 0 || typeID() == 1){
        std::ostringstream err;
        err<<"This matrix is EMPTY or REAL. If it's REAL, please use operator[] instead of operator().";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 2)
        return getElemAt(idx, cm_elem, ongpu), 0;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::operator()(size_t):");
    }
    return Complex();
  };

  Complex* Block::getElem(cflag tp)const{
    throwTypeError(tp);
    try{
      if(typeID() == 1){
        std::ostringstream err;
        err<<"This matrix is REAL. Please use getElem() or getElem(uni10::rflag ) instead.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      return cm_elem;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::getElem(uni10::cflag ):");
    }
    return cm_elem;
  }

  void Block::save(cflag tp, const std::string& fname)const{
    try{
      throwTypeError(tp);
      FILE *fp = fopen(fname.c_str(), "w");
      if(!(fp != NULL)){
        std::ostringstream err;
        err<<"Error in writing to file '"<<fname<<"'.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      fwrite(&r_flag, sizeof(r_flag), 1, fp);
      fwrite(&c_flag, sizeof(c_flag), 1, fp);
      fwrite(&Rnum, sizeof(Rnum), 1, fp);
      fwrite(&Cnum, sizeof(Cnum), 1, fp);
      fwrite(&diag, sizeof(diag), 1, fp);
      fwrite(&ongpu, sizeof(ongpu), 1, fp);
      Complex* elem = cm_elem;
      if(ongpu){
        elem = (Complex*)malloc(elemNum() * sizeof(Complex));
        elemCopy(elem, cm_elem, elemNum() * sizeof(Complex), false, ongpu);
      }
      fwrite(elem, sizeof(Complex), elemNum(), fp);
      if(ongpu)
        free(elem);
      fclose(fp);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::save(uni10::cflag, std::string&):");
    }
  }

  std::vector<Matrix> Block::qr(cflag tp)const{
    std::vector<Matrix> outs;
    try{
      throwTypeError(tp);
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QR decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(CTYPE, Cnum, Cnum, false, ongpu));
      if(!diag)
          matrixQR(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem, ongpu);
      else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(size_t i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixQR(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem, ongpu);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::qr(uni10::cflag ):");
    }
    return outs;
  }

  std::vector<Matrix> Block::rq(cflag tp)const{
    std::vector<Matrix> outs;
    try{
      throwTypeError(tp);
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(CTYPE, Rnum, Rnum, false, ongpu)); //r
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu)); //q
      if(!diag){
        matrixRQ(cm_elem, Rnum, Cnum, outs[1].cm_elem, outs[0].cm_elem, ongpu);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(size_t i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixRQ(tmpC, min, min, outs[1].cm_elem, outs[0].cm_elem, ongpu);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::rq(uni10::cflag ):");
    }
    return outs;
  }

  std::vector<Matrix> Block::lq(cflag tp)const{
    std::vector<Matrix> outs;
    try{
      throwTypeError(tp);
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(CTYPE, Rnum, Rnum, false, ongpu));
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      if(!diag){
        matrixLQ(cm_elem, Rnum, Cnum, outs[1].cm_elem, outs[0].cm_elem, ongpu);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(size_t i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixLQ(tmpC, min, min, outs[1].cm_elem, outs[0].cm_elem, ongpu);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::lq(uni10::cflag ):");
    }
    return outs;
  }

  std::vector<Matrix> Block::ql(cflag tp)const{
    std::vector<Matrix> outs;
    try{
      throwTypeError(tp);
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QL decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(CTYPE, Cnum, Cnum, false, ongpu));
      if(!diag){
        matrixQL(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem, ongpu);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(size_t i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixQL(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem, ongpu);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::ql(uni10::cflag ):");
    }
    return outs;
  }

  std::vector<Matrix> Block::svd(cflag tp)const{
    std::vector<Matrix> outs;
    try{
      throwTypeError(tp);
    /*
      if(diag){
        std::ostringstream err;
        err<<"Cannot perform singular value decomposition on a diagonal matrix. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    */
      size_t min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
      //GPU_NOT_READY
      outs.push_back(Matrix(CTYPE, Rnum, min, false, ongpu));
      outs.push_back(Matrix(CTYPE, min, min, true, ongpu));
      outs.push_back(Matrix(CTYPE, min, Cnum, false, ongpu));
      if(!diag){
        matrixSVD(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem, outs[2].cm_elem, ongpu);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(size_t i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixSVD(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem, outs[2].cm_elem, ongpu);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::svd(uni10::cflag ):");
    }
    return outs;
  }

  Real Block::norm(cflag tp)const{
    try{
      throwTypeError(tp);
      if(typeID() == 0){
        std::ostringstream err;
        err<<"This matirx is empty" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
        throw std::runtime_error(exception_msg(err.str()));
      }
      return vectorNorm(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::norm(uni10::cflag ):");
    }
    return 0;
  }

  Matrix Block::inverse(cflag tp)const{
    try{
      throwTypeError(tp);
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform inversion on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      Matrix invM(*this);
      assert(ongpu == invM.isOngpu());
      matrixInv(invM.cm_elem, Rnum, invM.diag, invM.ongpu);
      return invM;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::inverse(uni10::cflag ):");
      return Matrix();
    }
  }

  Matrix Block::getDiag(cflag tp)const{
    try{
      throwTypeError(tp);
      if(diag)
        return *this;
      else{
        Matrix D(CTYPE, Rnum, Cnum, true, ongpu);
        ::uni10::getDiag(cm_elem, D.getElem(CTYPE), Rnum, Cnum, D.elemNum(), ongpu, D.isOngpu());
        return D;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::getDiag(uni10::cflag ):");
      return Matrix();
    }
  }

  Complex Block::trace(cflag tp)const{
    try{
      throwTypeError(tp);
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform trace on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(diag)
        return vectorSum(cm_elem, elemNum(), 1, ongpu);
      else
        return vectorSum(cm_elem, Cnum, Cnum + 1, ongpu);
    }catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::trace(uni10::cflag ):");
      return 0;
    }
  }

  Complex Block::sum(cflag tp)const{
    try{
      throwTypeError(tp);
      return vectorSum(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::sum(uni10::cflag ):");
    }
    return 0;
  }

  std::vector<Matrix> Block::eig(cflag tp)const{
    std::vector<Matrix> outs;
    try{
      throwTypeError(tp);
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
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      eigDecompose(cm_elem, Rnum, outs[0].cm_elem, outs[1].cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eig(uni10::cflag ):");
    }
    return outs;
  }

  std::vector<Matrix> Block::eigh(cflag tp)const{
    std::vector<Matrix> outs;
    try{
      throwTypeError(tp);
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
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      Matrix Eig(RTYPE, Rnum, Cnum, true, ongpu);
      eigSyDecompose(cm_elem, Rnum, Eig.m_elem, outs[1].cm_elem, ongpu);
      outs[0] = Eig;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eigh(uni10::cflag ):");
    }
    return outs;
  }

  Complex Block::at(cflag tp, size_t r, size_t c)const{
    try{
      throwTypeError(tp);
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
        return getElemAt(r, cm_elem, ongpu);
      }
      else
        return getElemAt(r * Cnum + c, cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::at(uni10::cflag, size_t, size_t):");
    }
    return Complex();
  }

};	/* namespace uni10 */


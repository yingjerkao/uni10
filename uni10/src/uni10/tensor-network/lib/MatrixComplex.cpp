/****************************************************************************
*  @file CMakeLists.txt
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
*  @brief Implementation file for Matrix class
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/tools/uni10_tools.h>
#include <uni10/numeric/lapack/uni10_lapack.h>
#include <uni10/tensor-network/Matrix.h>


namespace uni10{

  void Matrix::init(cflag tp, bool _ongpu){

    if(elemNum()){
      if(_ongpu)	// Try to allocate GPU memory
        cm_elem = (Complex*)elemAlloc(elemNum() * sizeof(Complex), ongpu);
      else{
        cm_elem = (Complex*)elemAllocForce(elemNum() * sizeof(Complex), false);
        ongpu = false;
      }
    }
    m_elem = NULL;

  }

  void Matrix::init(const Complex* _elem, bool src_ongpu){

    init(CTYPE, ongpu);
    elemCopy(cm_elem, _elem, elemNum() * sizeof(Complex), ongpu, src_ongpu);

  }

Matrix::Matrix(cflag tp, size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu):Block(tp, _Rnum, _Cnum, _diag){
  try{
    throwTypeError(tp);
    init(tp, _ongpu);
    if(elemNum())
      elemBzero(cm_elem, elemNum() * sizeof(Complex), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(uni10::cflag, size_t, size_t, bool ,bool = false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const Complex* _elem, bool _diag, bool _ongpu, bool src_ongpu): Block(CTYPE, _Rnum, _Cnum, _diag){
  try{
    ongpu = _ongpu;
    init(_elem, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, double*, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const std::vector<Complex>& _elem, bool _diag, bool _ongpu, bool src_ongpu): Block(CTYPE, _Rnum, _Cnum, _diag){
  try{
    ongpu = _ongpu;
    if(_diag == false && _Rnum*_Cnum != _elem.size() ){
      std::ostringstream err;
      err<<"Number of the input elements is: " << _elem.size() <<", and it doesn't match to the size of matrix: "\
        << _Rnum*_Cnum << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if( _diag == true && std::min(_Rnum, _Cnum) != _elem.size()){
      std::ostringstream err;
      err<<"Number of the input elements is: " << _elem.size() <<", and it doesn't match to the min(Rnum, Cnum) of matrix: "\
        << std::min(_Rnum, _Cnum) << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    init(&_elem[0], src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, std::vector<Complex>&, bool=false):");
  }
}

Matrix::Matrix(cflag tp, const std::string& fname){
  try{
    throwTypeError(tp);
    FILE *fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    fread(&r_flag, sizeof(r_flag), 1, fp);
    fread(&c_flag, sizeof(c_flag), 1, fp);
    fread(&Rnum, sizeof(Rnum), 1, fp);
    fread(&Cnum, sizeof(Cnum), 1, fp);
    fread(&diag, sizeof(diag), 1, fp);
    fread(&ongpu, sizeof(ongpu), 1, fp);
    init(tp, ongpu);
    if(elemNum())
      elemBzero(cm_elem, elemNum() * sizeof(Complex), ongpu);
    Complex* elem = cm_elem;
    if(ongpu)
      elem = (Complex*)malloc(elemNum() * sizeof(Complex));
    fread(elem, sizeof(Complex), elemNum(), fp);
    if(ongpu){
      elemCopy(cm_elem, elem, elemNum() * sizeof(Complex), ongpu, false);
      free(elem);
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(uni10::cflag, std::string& )");
  }
}

void Matrix::setElem(const std::vector<Complex>& elem, bool src_ongpu){
  try{
    if(diag == false && Rnum*Cnum != elem.size() ){
      std::ostringstream err;
      err<<"Number of the input elements is: " << elem.size() <<", and it doesn't match to the size of matrix: "\
        << Rnum*Cnum << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag == true && std::min(Rnum, Cnum) != elem.size()){
      std::ostringstream err;
      err<<"Number of the input elements is: " << elem.size() <<", and it doesn't match to the min(Rnum, Cnum) of matrix: "\
        << std::min(Rnum, Cnum) << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    setElem(&elem[0], src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(std::vector<Complex>&, bool=false):");
  }
}

void Matrix::setElem(const Complex* elem, bool src_ongpu){
  try{
    if(typeID() == 1){
      r_flag = RNULL;
      c_flag = CTYPE;
      init(CTYPE, ongpu);
    }
    elemCopy(cm_elem, elem, elemNum() * sizeof(Complex), ongpu, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(const Complex*, bool=false):");
  }
}

void Matrix::identity(cflag tp){
  try{
    throwTypeError(tp);
    diag = true;
    if(cm_elem != NULL)
      elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);

    cm_elem = (Complex*)elemAlloc(elemNum() * sizeof(Complex), ongpu);
    Complex* elemI = (Complex*)malloc(elemNum() * sizeof(Complex));
    for(size_t i = 0; i < elemNum(); i++)
      elemI[i] = 1;
    this->setElem(elemI, false);
    free(elemI);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::identity(uni10::cflag ):");
  }
}

void Matrix::set_zero(cflag tp){
  try{
    throwTypeError(tp); 
    if(elemNum())
      elemBzero(cm_elem, elemNum() * sizeof(Complex), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::set_zero(uni10::cflag ):");
  }
}

void Matrix::randomize(cflag tp){
  try{
    throwTypeError(tp);
    if(!ongpu)
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    elemRand(cm_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::randomize(uni10::cflag ):");
  }
}

void Matrix::orthoRand(cflag tp){
  try{
    throwTypeError(tp);
    if(!ongpu)
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    if(!diag)
      orthoRandomize(cm_elem, Rnum, Cnum, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::orthoRand(uni10::cflag ):");
  }
}

Matrix& Matrix::normalize(cflag tp){
  try{
    throwTypeError(tp);
    Real norm = vectorNorm(cm_elem, elemNum(), 1, ongpu);
    vectorScal((1./norm), cm_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::normalize(uni10::cflag ):");
  }
  return *this;
}

Matrix& Matrix::transpose(cflag tp){
  try{
    throwTypeError(tp);
    if(!ongpu)
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    if(!(diag || Rnum == 1 || Cnum == 1))
      setTranspose(cm_elem, Rnum, Cnum, ongpu);
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose(uni10::cflag ):");
  }
  return *this;
}

Matrix& Matrix::cTranspose(cflag tp){
  try{
    throwTypeError(tp);
    if(!ongpu)
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    if(diag || Rnum == 1 || Cnum == 1)
      this->conj();
    else
      setCTranspose(cm_elem, Rnum, Cnum, ongpu);
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose(uni10::cflag ):");
  }
  return *this;
}

Matrix& Matrix::conj(cflag tp){
  throwTypeError(tp);
  setConjugate(cm_elem, elemNum(), ongpu);
  return *this;
}

Matrix& Matrix::resize(cflag tp, size_t row, size_t col){
  try{
    throwTypeError(tp);
    if(diag){
      size_t _elemNum = row < col ? row : col;
      if(_elemNum > elemNum()){
        bool des_ongpu;
        Complex* elem = (Complex*)elemAlloc(_elemNum * sizeof(Complex), des_ongpu);
        elemBzero(elem, _elemNum * sizeof(Complex), des_ongpu);
        elemCopy(elem, cm_elem, elemNum() * sizeof(Complex), des_ongpu, ongpu);
        if(cm_elem != NULL)
          elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
        cm_elem = elem;
        ongpu = des_ongpu;
      }
      else
        shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Complex), ongpu);
      Rnum = row;
      Cnum = col;
    }
    else{
      if(col == Cnum){
        size_t _elemNum = row * col;
        if(row > Rnum){
          bool des_ongpu;
          Complex* elem = (Complex*)elemAlloc(_elemNum * sizeof(Complex), des_ongpu);
          elemBzero(elem, _elemNum * sizeof(Complex), des_ongpu);
          elemCopy(elem, cm_elem, elemNum() * sizeof(Complex), des_ongpu, ongpu);
          if(cm_elem != NULL)
            elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
          cm_elem = elem;
          ongpu = des_ongpu;
        }
        else
          shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Complex), ongpu);
        Rnum = row;
      }
      else{
        size_t data_row = row < Rnum ? row : Rnum;
        size_t data_col = col < Cnum ? col : Cnum;
        bool des_ongpu;
        Complex* elem = (Complex*)elemAlloc(row * col * sizeof(Complex), des_ongpu);
        elemBzero(elem, row * col * sizeof(Complex), des_ongpu);
        for(size_t r = 0; r < data_row; r++)
          elemCopy(&(elem[r * col]), &(cm_elem[r * Cnum]), data_col * sizeof(Complex), des_ongpu, ongpu);
        if(cm_elem != NULL)
          elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
        cm_elem = elem;
        ongpu = des_ongpu;
        Rnum = row;
        Cnum = col;
      }
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::resize(uni10::cflag, size_t, size_t):");
  }
  return *this;
}

void Matrix::assign(cflag tp, size_t _Rnum, size_t _Cnum){
  try{
    throwTypeError(tp);
    Matrix M(CTYPE, _Rnum, _Cnum, diag, ongpu);
    *this = M;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::assign(uni10::cflag, size_t ,size_t ):");
  }
}

bool Matrix::toGPU(cflag tp){
  throwTypeError(tp);
  if(!ongpu)
    cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
  return ongpu;
}

Complex& Matrix::at(cflag tp, size_t idx){
  try{
    throwTypeError(tp);
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of the matrix elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    cm_elem = (Complex*)mvCPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    return cm_elem[idx];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::at(uni10::cflag , size_t ):");
    return cm_elem[0];
  }
}

Complex* Matrix::getHostElem(cflag tp){
  try{
    throwTypeError(tp);
    if(ongpu){
      cm_elem = (Complex*)mvCPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::getHostElem(uni10::cflag ):");
  }
  return cm_elem;
}

};	/* namespace uni10 */

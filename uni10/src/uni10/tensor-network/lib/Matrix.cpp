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
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/tensor-network/CMatrix.h>

typedef double Real;
typedef std::complex<double> Complex;

namespace uni10{
Matrix::Matrix(): Block(){}

void Matrix::matrixElemFree(){
  if(m_elem != NULL)
    elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
  if(cm_elem != NULL)
    elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
}

void Matrix::init(const Real* _elem, bool src_ongpu){
  m_type = REAL;
  init(true, REAL);
  elemCopy(m_elem, _elem, elemNum() * sizeof(Real), ongpu, src_ongpu);
}

void Matrix::init(const Complex* _elem, bool src_ongpu){
  m_type = COMPLEX;
  init(true, COMPLEX);
  elemCopy(cm_elem, _elem, elemNum() * sizeof(Complex), ongpu, src_ongpu);
}

void Matrix::init(bool _ongpu, matrixType tp){
  if(tp == REAL){
    if(elemNum()){
      if(_ongpu)	// Try to allocate GPU memory
        m_elem = (Real*)elemAlloc(elemNum() * sizeof(Real), ongpu);
      else{
        m_elem = (Real*)elemAllocForce(elemNum() * sizeof(Real), false);
        ongpu = false;
      }
    }
  }
  if(tp == COMPLEX){
    if(elemNum()){
      if(_ongpu)	// Try to allocate GPU memory
        cm_elem = (Complex*)elemAlloc(elemNum() * sizeof(Complex), ongpu);
      else{
        cm_elem = (Complex*)elemAllocForce(elemNum() * sizeof(Complex), false);
        ongpu = false;
      }
    }
  }
}

void Matrix::init(const Real* _m_elem, const Complex* _cm_elem, bool src_ongpu){
  assert(!(_m_elem != NULL && _cm_elem != NULL));
  if(_m_elem != NULL){
    init(_m_elem, src_ongpu);
  }
  if(_cm_elem != NULL){
    init(_cm_elem, src_ongpu);
  }
}

Matrix::Matrix(const Matrix& _m): Block(_m.Rnum, _m.Cnum, _m.diag){
  try{
    init(_m.m_elem, _m.cm_elem, _m.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(uni10::Matrix&):");
  }
}

Matrix::Matrix(const Block& _b): Block(_b){
  try{
    init(_b.m_elem, _b.cm_elem, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(uni10::Block&):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const Real* _elem, bool _diag, bool src_ongpu): Block(_Rnum, _Cnum, _diag){
  try{
    init(_elem, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, double*, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const Complex* _elem, bool _diag, bool src_ongpu): Block(_Rnum, _Cnum, _diag){
  try{
    init(_elem, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, double*, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const std::vector<Real>& _elem, bool _diag, bool src_ongpu): Block(_Rnum, _Cnum, _diag){
  try{
    if(_diag == false && _Rnum*_Cnum != _elem.size()){
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
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, std::vector<double>&, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const std::vector<Complex>& _elem, bool _diag, bool src_ongpu): Block(_Rnum, _Cnum, _diag){
  try{
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

Matrix::Matrix(matrixType tp, size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu):Block(tp, _Rnum, _Cnum, _diag){
  try{
    if(tp == REAL){
      init(_ongpu, tp);
      if(elemNum())
        elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
    }
    if(tp == COMPLEX){
      init(_ongpu, tp);
      if(elemNum())
        elemBzero(cm_elem, elemNum() * sizeof(Complex), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, bool=false):");
  }
}

Matrix::Matrix(const std::string& fname){
  try{
    FILE *fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    fread(&m_type, sizeof(m_type), 1, fp);
    fread(&Rnum, sizeof(Rnum), 1, fp);
    fread(&Cnum, sizeof(Cnum), 1, fp);
    fread(&diag, sizeof(diag), 1, fp);
    fread(&ongpu, sizeof(ongpu), 1, fp);
    if(m_type == REAL){
      init(ongpu, REAL);
      if(elemNum())
        elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
    }
    if(m_type == COMPLEX){
      init(ongpu, COMPLEX);
      if(elemNum())
        elemBzero(cm_elem, elemNum() * sizeof(Complex), ongpu);
    }
    if(m_type == REAL){
      Real* elem = m_elem;
      if(ongpu)
        elem = (Real*)malloc(elemNum() * sizeof(Real));
      fread(elem, sizeof(Real), elemNum(), fp);
      if(ongpu){
        elemCopy(m_elem, elem, elemNum() * sizeof(Real), ongpu, false);
        free(elem);
      }
    }
    if(m_type == COMPLEX){
      Complex* elem = cm_elem;
      if(ongpu)
        elem = (Complex*)malloc(elemNum() * sizeof(Complex));
      fread(elem, sizeof(Complex), elemNum(), fp);
      if(ongpu){
        elemCopy(cm_elem, elem, elemNum() * sizeof(Complex), ongpu, false);
        free(elem);
      }
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(std::string& )");
  }
}

Matrix& Matrix::operator=(const Matrix& _m){
  try{
    Rnum = _m.Rnum;
    Cnum = _m.Cnum;
    diag = _m.diag;
    matrixElemFree();
    init(_m.m_elem, _m.cm_elem, _m.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator=(uni10::Matrix&):");
  }
  return *this;
}

Matrix& Matrix::operator=(const Block& _b){
  try{
    Rnum = _b.Rnum;
    Cnum = _b.Cnum;
    diag = _b.diag;
    matrixElemFree();
    init(_b.m_elem, _b.cm_elem, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator=(uni10::Block&):");
  }
  return *this;
}

Matrix::~Matrix(){
  try{
    matrixElemFree();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In destructor Matrix::~Matrix():");
  }
}

void Matrix::setElem(const std::vector<Real>& elem, bool _ongpu){
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
    setElem(&elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(std::vector<double>&, bool=false):");
  }
}

void Matrix::setElem(const Real* elem, bool _ongpu){
  try{
    if(m_type == COMPLEX){
      std::ostringstream err;
      err<<"Can't set REAL elements in a COMPLEX matrix" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(m_type == EMPTY)
      init(elem, _ongpu);
    m_type = REAL;
    elemCopy(m_elem, elem, elemNum() * sizeof(Real), ongpu, _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(double*, bool=false):");
  }
}

void Matrix::setElem(const std::vector<Complex>& elem, bool _ongpu){
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
    setElem(&elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(std::vector<Complex>&, bool=false):");
  }
}

void Matrix::setElem(const Complex* elem, bool _ongpu){
  assert(m_type == EMPTY || m_type == COMPLEX);
  try{
    if(m_type == REAL){
      std::ostringstream err;
      err<<"Can't set Complex elements in a REAL matrix" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(m_type == EMPTY)
      init(elem, _ongpu);
    m_type = COMPLEX;
    elemCopy(cm_elem, elem, elemNum() * sizeof(Complex), ongpu, _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(Complex*, bool=false):");
  }
}

void Matrix::identity(){
  try{
    if(m_type == EMPTY){
      std::ostringstream err;
      err<<"Haven't defined the type of matrix." << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    diag = true;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
    if(cm_elem != NULL)
      elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
    
    if(m_type == REAL){
      m_elem = (Real*)elemAlloc(elemNum() * sizeof(Real), ongpu);
      Real* elemI = (Real*)malloc(elemNum() * sizeof(Real));
      for(int i = 0; i < elemNum(); i++)
        elemI[i] = 1;
      this->setElem(elemI, false);
    }
    if(m_type == COMPLEX){
      cm_elem = (Complex*)elemAlloc(elemNum() * sizeof(Complex), ongpu);
      Complex* elemI = (Complex*)malloc(elemNum() * sizeof(Complex));
      for(int i = 0; i < elemNum(); i++)
        elemI[i] = 1;
      this->setElem(elemI, false);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::identity():");
  }
}

void Matrix::set_zero(){
  try{
    if(m_type == EMPTY){
      std::ostringstream err;
      err<<"Haven't defined the type of matrix." << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(m_type == REAL){
      if(elemNum())
        elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
    }
    if(m_type == COMPLEX){
      if(elemNum())
        elemBzero(cm_elem, elemNum() * sizeof(Complex), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::set_zero():");
  }
}

void Matrix::randomize(){
  try{
    if(m_type == EMPTY){
      std::ostringstream err;
      err<<"Haven't defined the type of matrix." << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(m_type == REAL){
      if(!ongpu)
        m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
      elemRand(m_elem, elemNum(), ongpu);
    }
    if(m_type == COMPLEX){
      if(!ongpu)
        cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
      elemRand(cm_elem, elemNum(), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::randomize():");
  }
}

void Matrix::orthoRand(){
  try{
    if(m_type == EMPTY){
      std::ostringstream err;
      err<<"Haven't defined the type of matrix." << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(m_type == REAL){
      if(!ongpu)
        m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
      if(!diag)
        orthoRandomize(m_elem, Rnum, Cnum, ongpu);
    }
    if(m_type == COMPLEX){
      if(!ongpu)
        cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
      if(!diag)
        orthoRandomize(cm_elem, Rnum, Cnum, ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::orthoRand():");
  }
}

Matrix& Matrix::transpose(){
  assert(m_type == REAL || m_type == COMPLEX);
  try{
    if(m_type == REAL){
      if(!ongpu)
        m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
      if(!(diag || Rnum == 1 || Cnum == 1))
        setTranspose(m_elem, Rnum, Cnum, ongpu);
    }
    if(m_type == COMPLEX){
      if(!ongpu)
        cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
      if(!(diag || Rnum == 1 || Cnum == 1))
        setTranspose(cm_elem, Rnum, Cnum, ongpu);
    }
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose():");
  }
  return *this;
}

Matrix& Matrix::cTranspose(){
  assert(m_type == REAL || m_type == COMPLEX);
  try{
    if(m_type == REAL){
      if(!ongpu)
        m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
      if(diag || Rnum == 1 || Cnum == 1)
        this->conj();
      else
        setCTranspose(m_elem, Rnum, Cnum, ongpu);
    }
    if(m_type == COMPLEX){
      if(!ongpu)
        cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
      if(diag || Rnum == 1 || Cnum == 1)
        this->conj();
      else
        setCTranspose(cm_elem, Rnum, Cnum, ongpu);
    }
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose():");
  }
  return *this;
}

Matrix& Matrix::resize(size_t row, size_t col){
  try{
    if(diag){
      size_t _elemNum = row < col ? row : col;
      if(_elemNum > elemNum()){
        bool des_ongpu;
        if(m_type == REAL){
          Real* elem = (Real*)elemAlloc(_elemNum * sizeof(Real), des_ongpu);
          elemBzero(elem, _elemNum * sizeof(Real), des_ongpu);
          elemCopy(elem, m_elem, elemNum() * sizeof(Real), des_ongpu, ongpu);
          if(m_elem != NULL)
            elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
          m_elem = elem;
        }
        if(m_type == COMPLEX){
          Complex* elem = (Complex*)elemAlloc(_elemNum * sizeof(Complex), des_ongpu);
          elemBzero(elem, _elemNum * sizeof(Complex), des_ongpu);
          elemCopy(elem, cm_elem, elemNum() * sizeof(Complex), des_ongpu, ongpu);
          if(cm_elem != NULL)
            elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
          cm_elem = elem;
        }
        ongpu = des_ongpu;
      }
      else{
        if(m_type == REAL)
          shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Real), ongpu);
        if(m_type == COMPLEX)
          shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Complex), ongpu);
      }
      Rnum = row;
      Cnum = col;
    }
    else{
      if(col == Cnum){
        size_t _elemNum = row * col;
        if(row > Rnum){
          bool des_ongpu;
          if(m_type == REAL){
            Real* elem = (Real*)elemAlloc(_elemNum * sizeof(Real), des_ongpu);
            elemBzero(elem, _elemNum * sizeof(Real), des_ongpu);
            elemCopy(elem, m_elem, elemNum() * sizeof(Real), des_ongpu, ongpu);
            if(m_elem != NULL)
              elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
            m_elem = elem;
          }
          if(m_type == COMPLEX){
            Complex* elem = (Complex*)elemAlloc(_elemNum * sizeof(Complex), des_ongpu);
            elemBzero(elem, _elemNum * sizeof(Complex), des_ongpu);
            elemCopy(elem, cm_elem, elemNum() * sizeof(Complex), des_ongpu, ongpu);
            if(cm_elem != NULL)
              elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
            cm_elem = elem;
          }
          ongpu = des_ongpu;
        }
        else{
          if(m_type == REAL)
            shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Real), ongpu);
          if(m_type == REAL)
            shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Complex), ongpu);
        }
        Rnum = row;
      }
      else{
        size_t data_row = row < Rnum ? row : Rnum;
        size_t data_col = col < Cnum ? col : Cnum;
        bool des_ongpu;
        if(m_type == REAL){
          Real* elem = (Real*)elemAlloc(row * col * sizeof(Real), des_ongpu);
          elemBzero(elem, row * col * sizeof(Real), des_ongpu);
          for(size_t r = 0; r < data_row; r++)
            elemCopy(&(elem[r * col]), &(m_elem[r * Cnum]), data_col * sizeof(Real), des_ongpu, ongpu);
          if(m_elem != NULL)
            elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
          m_elem = elem;
        }
        if(m_type == COMPLEX){
          Complex* elem = (Complex*)elemAlloc(row * col * sizeof(Complex), des_ongpu);
          elemBzero(elem, row * col * sizeof(Complex), des_ongpu);
          for(size_t r = 0; r < data_row; r++)
            elemCopy(&(elem[r * col]), &(cm_elem[r * Cnum]), data_col * sizeof(Complex), des_ongpu, ongpu);
          if(cm_elem != NULL)
            elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
          cm_elem = elem;
        }
        ongpu = des_ongpu;
        Rnum = row;
        Cnum = col;
      }
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::resize(size_t, size_t):");
  }
  return *this;
}

bool Matrix::toGPU(){
  if(m_type == REAL){
    if(!ongpu)
      m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
  }
  if(m_type == COMPLEX){
    if(!ongpu)
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
  }
  return ongpu;
}

double Matrix::max(bool on_gpu){
  try{
    if(m_type != REAL){
      std::ostringstream err;
      if(m_type == COMPLEX)
        err<< "Can't comparision. The type of matirx is COMPLEX" << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      if(m_type == EMPTY)
        err<< "Can't comparision. The type of matirx is EMPTY" << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return elemMax(m_elem,  elemNum(),  on_gpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::max(bool ):");
  }
}

void Matrix::load(const std::string& fname){
  try{
    FILE *fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    fread(&m_type, sizeof(m_type), 1, fp);
    fread(&Rnum, sizeof(Rnum), 1, fp);
    fread(&Cnum, sizeof(Cnum), 1, fp);
    fread(&diag, sizeof(diag), 1, fp);
    fread(&ongpu, sizeof(ongpu), 1, fp);
    if(m_type == REAL){
      init(ongpu, REAL);
      if(elemNum())
        elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
    }
    if(m_type == COMPLEX){
      init(ongpu, COMPLEX);
      if(elemNum())
        elemBzero(cm_elem, elemNum() * sizeof(Complex), ongpu);
    }
    if(m_type == REAL){
      Real* elem = m_elem;
      if(ongpu)
        elem = (Real*)malloc(elemNum() * sizeof(Real));
      fread(elem, sizeof(Real), elemNum(), fp);
      if(ongpu){
        elemCopy(m_elem, elem, elemNum() * sizeof(Real), ongpu, false);
        free(elem);
      }
    }
    if(m_type == COMPLEX){
      Complex* elem = cm_elem;
      if(ongpu)
        elem = (Complex*)malloc(elemNum() * sizeof(Complex));
      fread(elem, sizeof(Complex), elemNum(), fp);
      if(ongpu){
        elemCopy(cm_elem, elem, elemNum() * sizeof(Complex), ongpu, false);
        free(elem);
      }
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::load(std::string&):");
  }
}

/********************************************************************************************/

Matrix::Matrix(size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu): Block(_Rnum, _Cnum, _diag){
/*
  try{
    init(_ongpu, REAL);
    if(elemNum())
      elemBzero(m_elem, elemNum() * sizeof(double), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, bool=false):");
  }
*/
}

Matrix& Matrix::operator*= (const Block& Mb){
  try{
    if(!ongpu)
      m_elem = (double*)mvGPU(m_elem, elemNum() * sizeof(double), ongpu);
    *this = *this * Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator*=(uni10::Matrix&):");
  }
  return *this;
}

Matrix& Matrix::operator*= (double a){
  try{
    if(!ongpu)
      m_elem = (double*)mvGPU(m_elem, elemNum() * sizeof(double), ongpu);
    vectorScal(a, m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator*=(double):");
  }
	return *this;
}


Matrix& Matrix::operator+= (const Block& Mb){
  try{
    if(!ongpu)
      m_elem = (double*)mvGPU(m_elem, elemNum() * sizeof(double), ongpu);
    if (diag && !Mb.diag) {
      Matrix Mc(Rnum,Cnum);
      setDiag(Mc.m_elem,m_elem,Mc.Rnum,Mc.Cnum,Rnum,Mc.ongpu,ongpu);
      this->resize(Rnum,Cnum);
      *this=Mc;
      vectorAdd(m_elem, Mb.m_elem, elemNum(), ongpu, Mb.ongpu);
    } else if (Mb.diag && !diag) {
      Matrix Mc(Mb.Rnum,Mb.Cnum);
      setDiag(Mc.m_elem,Mb.m_elem,Mc.Rnum,Mc.Cnum,Mb.Cnum,Mc.ongpu,Mb.ongpu);
      vectorAdd(m_elem, Mc.m_elem, Mc.elemNum(), Mc.ongpu, ongpu);
    } else {
    vectorAdd(m_elem, Mb.m_elem, elemNum(), ongpu, Mb.ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator+=(uni10::Matrix&):");
  }
  return *this;
}

double& Matrix::operator[](size_t idx){
  try{
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of the matrix elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    m_elem = (double*)mvCPU(m_elem, elemNum() * sizeof(double), ongpu);
    return m_elem[idx];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::opeartor[](size_t):");
    return m_elem[0];
  }
}

double* Matrix::getHostElem(){
  try{
    if(ongpu){
      m_elem = (double*)mvCPU(m_elem, elemNum() * sizeof(double), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::getHostElem():");
  }
	return m_elem;
}

double& Matrix::at(size_t r, size_t c){
  try{
    if(!((r < Rnum) && (c < Cnum))){
      std::ostringstream err;
      err<<"The input indices are out of range.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    m_elem = (double*)mvCPU(m_elem, elemNum() * sizeof(double), ongpu);
    if(diag){
      if(!(r == c && r < elemNum())){
        std::ostringstream err;
        err<<"The matrix is diagonal, there is no off-diagonal element.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      return m_elem[r];
    }
    else
      return m_elem[r * Cnum + c];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::at(size_t, size_t):");
    return m_elem[0];
  }
}
/**********************************************************************************/

Matrix exp(double a, const Block& mat){
  try{
    std::vector<Matrix> rets = mat.eig();
    Matrix Uinv = rets[1].inverse();
    vectorExp(a, rets[0].getComplexElem(), rets[0].row(), rets[0].isOngpu());
    return Uinv * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function exp(double, uni10::Matrix&):");
    return Matrix();
  }
}

Matrix exph(double a, const Block& mat){
  try{
    std::vector<Matrix> rets = mat.eigh();
    Matrix UT(rets[1]);
    UT.cTranspose();
    if(mat.m_type == REAL)
      vectorExp(a, rets[0].getRealElem(), rets[0].row(), rets[0].isOngpu());
    if(mat.m_type == COMPLEX)
      vectorExp(a, rets[0].getComplexElem(), rets[0].row(), rets[0].isOngpu());
    return UT * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function exph(double, uni10::Matrix&):");
    return Matrix();
  }
}

Matrix exp(const std::complex<double>& a, const Block& mat){
  try{
    std::vector<Matrix> rets = mat.eig();
    Matrix Uinv = rets[1].inverse();
    vectorExp(a, rets[0].getComplexElem(), rets[0].row(), rets[0].isOngpu());
    return Uinv * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function exp(std::complex<double>, uni10::Matrix&):");
    return Matrix();
  }
}

Matrix exp(const Block& mat){
  return exp(1.0, mat);
}

Matrix exph(const Block& mat){
  return exph(1.0, mat);
}

};	/* namespace uni10 */

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

void Matrix::init(rflag tp, bool _ongpu){

  if(elemNum()){
    if(_ongpu)	// Try to allocate GPU memory
      m_elem = (Real*)elemAlloc(elemNum() * sizeof(Real), ongpu);
    else{
      m_elem = (Real*)elemAllocForce(elemNum() * sizeof(Real), false);
      ongpu = false;
    }
  }
  cm_elem = NULL;

}

void Matrix::init(const Real* _elem, bool src_ongpu){
  init(RTYPE, ongpu);
  elemCopy(m_elem, _elem, elemNum() * sizeof(Real), ongpu, src_ongpu);
}

Matrix::Matrix(rflag tp, size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu):Block(tp, _Rnum, _Cnum, _diag){
  try{
    throwTypeError(tp);
    init(tp, _ongpu);
    if(elemNum())
      elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(uni10::rflag, size_t, size_t, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const Real* _elem, bool _diag, bool _ongpu, bool src_ongpu): Block(RTYPE, _Rnum, _Cnum, _diag){
  try{
    ongpu = _ongpu;
    init(_elem, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, double*, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const std::vector<Real>& _elem, bool _diag, bool _ongpu, bool src_ongpu): Block(RTYPE, _Rnum, _Cnum, _diag){
  try{
    ongpu = _ongpu;
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

Matrix::Matrix(rflag tp, const std::string& fname){
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
      elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
    Real* elem = m_elem;
    if(ongpu)
      elem = (Real*)malloc(elemNum() * sizeof(Real));
    fread(elem, sizeof(Real), elemNum(), fp);
    if(ongpu){
      elemCopy(m_elem, elem, elemNum() * sizeof(Real), ongpu, false);
      free(elem);
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(uni10::rflag ,std::string& )");
  }
}

void Matrix::setElem(const std::vector<Real>& elem, bool src_ongpu){
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
    propogate_exception(e, "In function Matrix::setElem(std::vector<double>&, bool=false):");
  }
}

void Matrix::setElem(const Real* elem, bool src_ongpu){
  try{
    if(typeID() == 2){
      r_flag = RTYPE;
      c_flag = CNULL;
      init(RTYPE, ongpu);
    }
    elemCopy(m_elem, elem, elemNum() * sizeof(Real), ongpu, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(const double*, bool=false):");
  }
}

void Matrix::identity(rflag tp){
  try{
    throwTypeError(tp);
    diag = true;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(Real), ongpu);

    m_elem = (Real*)elemAlloc(elemNum() * sizeof(Real), ongpu);
    Real* elemI = (Real*)malloc(elemNum() * sizeof(Real));
    
    for(size_t i = 0; i < elemNum(); i++)
      elemI[i] = 1;
    this->setElem(elemI);
    free(elemI);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::identity(uni10::rflag ):");
  }
}

void Matrix::set_zero(rflag tp){
  try{
    throwTypeError(tp);
    if(elemNum())
      elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::set_zero(uni10::flag ):");
  }
}

void Matrix::randomize(rflag tp){
  try{
    throwTypeError(tp);
    if(!ongpu)
      m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
    elemRand(m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::randomize(uni10::rflag ):");
  }
}

void Matrix::orthoRand(rflag tp){
  try{
    throwTypeError(tp);
    if(!ongpu)
      m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
    if(!diag)
      orthoRandomize(m_elem, Rnum, Cnum, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::orthoRand(uni10::rflag ):");
  }
}

Matrix& Matrix::transpose(rflag tp){
  try{
    throwTypeError(tp);
    if(!ongpu)
      m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
    if(!(diag || Rnum == 1 || Cnum == 1))
      setTranspose(m_elem, Rnum, Cnum, ongpu);
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose(uni10::rflag ):");
  }
  return *this;
}

Matrix& Matrix::cTranspose(rflag tp){
  try{
    transpose(tp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::cTranspose(uni10::rflag ):");
  }
  return *this;
}

Matrix& Matrix::conj(rflag tp){
  throwTypeError(tp);
  return *this;
}

Matrix& Matrix::resize(rflag tp, size_t row, size_t col){
  try{
    throwTypeError(tp);
    if(diag){
      size_t _elemNum = row < col ? row : col;
      if(_elemNum > elemNum()){
        bool des_ongpu;
        Real* elem = (Real*)elemAlloc(_elemNum * sizeof(Real), des_ongpu);
        elemBzero(elem, _elemNum * sizeof(Real), des_ongpu);
        elemCopy(elem, m_elem, elemNum() * sizeof(Real), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
        m_elem = elem;
        ongpu = des_ongpu;
      }
      else
        shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Real), ongpu);
      Rnum = row;
      Cnum = col;
    }
    else{
      if(col == Cnum){
        size_t _elemNum = row * col;
        if(row > Rnum){
          bool des_ongpu;
          Real* elem = (Real*)elemAlloc(_elemNum * sizeof(Real), des_ongpu);
          elemBzero(elem, _elemNum * sizeof(Real), des_ongpu);
          elemCopy(elem, m_elem, elemNum() * sizeof(Real), des_ongpu, ongpu);
          if(m_elem != NULL)
            elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
          m_elem = elem;
          ongpu = des_ongpu;
        }
        else
          shrinkWithoutFree((elemNum() - _elemNum) * sizeof(Real), ongpu);
        Rnum = row;
      }
      else{
        size_t data_row = row < Rnum ? row : Rnum;
        size_t data_col = col < Cnum ? col : Cnum;
        bool des_ongpu;
        Real* elem = (Real*)elemAlloc(row * col * sizeof(Real), des_ongpu);
        elemBzero(elem, row * col * sizeof(Real), des_ongpu);
        for(size_t r = 0; r < data_row; r++)
          elemCopy(&(elem[r * col]), &(m_elem[r * Cnum]), data_col * sizeof(Real), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
        m_elem = elem;
        ongpu = des_ongpu;
        Rnum = row;
        Cnum = col;
      }
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::resize(uni10::rflag, size_t, size_t):");
  }
  return *this;
}

void Matrix::assign(rflag tp, size_t _Rnum, size_t _Cnum){
  try{
    throwTypeError(tp);
    Matrix M(RTYPE, _Rnum, _Cnum, diag, ongpu);
    *this = M;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::assign(uni10::rflag, size_t ,size_t ):");
  }
}

bool Matrix::toGPU(rflag tp){
  throwTypeError(tp);
  if(!ongpu)
    m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
  return ongpu;
}

Real& Matrix::at(rflag tp, size_t idx){
  try{
    throwTypeError(tp);
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of the matrix elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    m_elem = (Real*)mvCPU(m_elem, elemNum() * sizeof(Real), ongpu);
    return m_elem[idx];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::at(uni10::rflag ,size_t ):");
  }
  return m_elem[0];
}

Real* Matrix::getHostElem(rflag tp){
  try{
    throwTypeError(tp);
    if(ongpu)
      m_elem = (Real*)mvCPU(m_elem, elemNum() * sizeof(Real), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::getHostElem(uni10::rflag ):");
  }
  return m_elem;
}

/*********************  developping  **********************/

Real Matrix::max(rflag tp, bool on_gpu){
  try{
    throwTypeError(tp);
    return elemMax(m_elem,  elemNum(),  on_gpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::max(uni10::rflag, bool ):");
  }
  return 0;
}

Real Matrix::absMax(rflag tp, bool on_gpu){
  try{
    throwTypeError(tp);
    return elemAbsMax(m_elem,  elemNum(),  on_gpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::absMax(uni10::rflag, bool ):");
    return 0;
  }
}

Matrix& Matrix::normalize(rflag tp){
  try{
    throwTypeError(tp);
    Real norm = vectorNorm(m_elem, elemNum(), 1, ongpu);
    vectorScal((1./norm), m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::normalize(uni10::rflag ):");
  }
  return *this;
}

Matrix& Matrix::maxNorm(rflag tp){
  try{
    throwTypeError(tp);
    Real max = elemMax(m_elem, elemNum(), ongpu);
    vectorScal((1./max), m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::maxNorm(uni10::rflag ):");
  }
  return *this;
}

Matrix& Matrix::absMaxNorm(rflag tp){
  try{
    throwTypeError(tp);
    Real absMax = elemAbsMax(m_elem, elemNum(), ongpu);
    vectorScal((1./absMax), m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::absMaxNorm(uni10::rflag ):");
  }
  return *this;
}


};	/* namespace uni10 */

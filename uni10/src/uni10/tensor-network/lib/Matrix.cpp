/****************************************************************************
*  @file Matrix.cpp
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

Matrix& Matrix::operator=(const Matrix& _m){
  try{
    r_flag = _m.r_flag;
    c_flag = _m.c_flag;
    Rnum = _m.Rnum;
    Cnum = _m.Cnum;
    diag = _m.diag;
    ongpu = _m.ongpu;
    MelemFree();
    setMelemBNULL();
    init(_m.m_elem, _m.cm_elem, _m.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator=(uni10::Matrix&):");
  }
  return *this;
}

Matrix& Matrix::operator=(const Block& _b){
  try{
    r_flag = _b.r_flag;
    c_flag = _b.c_flag;
    Rnum = _b.Rnum;
    Cnum = _b.Cnum;
    diag = _b.diag;
    MelemFree();
    setMelemBNULL();
    init(_b.m_elem, _b.cm_elem, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator=(uni10::Block&):");
  }
  return *this;
}

Matrix& Matrix::operator*= (const Block& Mb){
  try{
    if(!ongpu && typeID() == 1)
      m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
    if(!ongpu && typeID() == 2)
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    *this = *this * Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator*=(uni10::Matrix&):");
  }
  return *this;
}

Matrix& Matrix::operator*= (Real a){
  try{
    if(!ongpu && typeID() == 1){
      m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
      vectorScal(a, m_elem, elemNum(), ongpu);
    }
    if(!ongpu && typeID() == 2){
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
      vectorScal(a, cm_elem, elemNum(), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator*=(double):");
  }
  return *this;
}

Matrix& Matrix::operator*= (Complex a){
  try{
    if(a.imag() == 0)
      *this = *this * a.real();
    else{
      if(typeID() == 1)
        RtoC(*this);
      if(!ongpu){
        cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
        vectorScal(a, cm_elem, elemNum(), ongpu);
      }
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator*=(double):");
  }
  return *this;
}

Matrix& Matrix::operator+= (const Block& Mb){
  try{
    if ((Rnum != Mb.Rnum) || (Cnum != Mb.Cnum) ){
      std::ostringstream err;
      err<<"These two matrix has different shape do not match for matrix add. ";
      throw std::runtime_error(exception_msg(err.str())); ;
    };
    if(!ongpu && typeID() == 1)
      m_elem = (Real*)mvGPU(m_elem, elemNum() * sizeof(Real), ongpu);
    if(!ongpu && typeID() == 2)
      cm_elem = (Complex*)mvGPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
    if(typeID() == 1 && Mb.typeID() == 1)
      RAddR(*this, Mb);
    else if(typeID() == 2 && Mb.typeID() == 2)
      CAddC(*this, Mb);
    else if(typeID() == 1 && Mb.typeID() == 2)
      RAddC(*this,Mb);
    else if(typeID() == 2 && Mb.typeID() == 1)
      CAddR(*this, Mb);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator+=(uni10::Matrix&):");
  }
  return *this;
}

Real& Matrix::operator[](size_t idx){
  try{
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of the matrix elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 2){
      std::ostringstream err;
      err<<"This matrix is COMPLEX. Please use operator()." << std::endl << "In the file Matrix.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    else if(typeID() == 1)
      m_elem = (Real*)mvCPU(m_elem, elemNum() * sizeof(Real), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::opeartor[](size_t):");
  }
  return m_elem[idx];
}

Complex& Matrix::operator()(size_t idx){
  try{
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of the matrix elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 1){
      std::ostringstream err;
      err<<"This matrix is REAL. Please use operator[] instead." << std::endl << "In the file Matrix.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 2)
      cm_elem = (Complex*)mvCPU(cm_elem, elemNum() * sizeof(Complex), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::opeartor()(size_t):");
  }
  return cm_elem[idx];
}

/*********************  NO TYPE **************************/

void Matrix::MelemFree(){
  if(m_elem != NULL)
    elemFree(m_elem, elemNum() * sizeof(Real), ongpu);
  if(cm_elem != NULL)
    elemFree(cm_elem, elemNum() * sizeof(Complex), ongpu);
}
void Matrix::setMelemBNULL(){
  m_elem = NULL;
  cm_elem = NULL;
}

void Matrix::init(const Real* _m_elem, const Complex* _cm_elem, bool src_ongpu){
  assert(!(_m_elem != NULL && _cm_elem != NULL));
  if(_m_elem != NULL)
    init(_m_elem, src_ongpu);
  else if(_cm_elem != NULL)
    init(_cm_elem, src_ongpu);
}

Matrix::Matrix(): Block(){}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu): Block(RTYPE, _Rnum, _Cnum, _diag){
  try{
    init(RTYPE, _ongpu);
    if(elemNum())
      elemBzero(m_elem, elemNum() * sizeof(Real), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, bool=false):");
  }
}

Matrix::Matrix(std::string stp, size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu){
  try{
    if(stp == "R"){
      Matrix tmp(RTYPE, _Rnum, _Cnum, _diag, _ongpu);
      *this = tmp;
    }else if(stp == "C"){
      Matrix tmp(CTYPE, _Rnum, _Cnum, _diag, _ongpu);
      *this = tmp;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(std::string , size_t, size_t, bool = false, bool=false):");
  }
}

Matrix::Matrix(const Matrix& _m): Block(_m.Rnum, _m.Cnum, _m.diag){
  try{
    r_flag = _m.r_flag;
    c_flag = _m.c_flag;
    ongpu = _m.ongpu;
    init(_m.m_elem, _m.cm_elem, _m.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(const uni10::Matrix& ):");
  }
}

Matrix::Matrix(const Block& _b): Block(_b){
  try{
    init(_b.m_elem, _b.cm_elem, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(const uni10::Block&):");
  }
}

Matrix::Matrix(const std::string& fname):Block(){
  try{
    if(typeID() == 1){
      Matrix tmp(RTYPE, fname);
      *this = tmp;
    }else if(typeID() == 2){
      Matrix tmp(CTYPE, fname);
      *this = tmp;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(std::string& )");
  }
}

Matrix::~Matrix(){
  try{
    MelemFree();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In destructor Matrix::~Matrix():");
  }
}

void Matrix::identity(){
  try{
    if(typeID() == 1)
      identity(RTYPE);
    else if(typeID() == 2)
      identity(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::identity():");
  }
}

void Matrix::set_zero(){
  try{
    if(typeID() == 1)
      set_zero(RTYPE);
    else if(typeID() == 2)
      set_zero(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::set_zero():");
  }
}

void Matrix::randomize(){
  try{
    if(typeID() == 1)
      randomize(RTYPE);
    else if(typeID() == 2)
      randomize(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::randomize():");
  }
}

void Matrix::orthoRand(){
  try{
    if(typeID() == 1)
      orthoRand(RTYPE);
    else if(typeID() == 2)
      orthoRand(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::orthoRand():");
  }
}

Matrix& Matrix::transpose(){
  try{
    if(typeID() == 1)
      return transpose(RTYPE);
    else if(typeID() == 2)
      return transpose(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose():");
  }
  return *this;
}

Matrix& Matrix::cTranspose(){
  try{
    if(typeID() == 1)
     return cTranspose(RTYPE);
    else if(typeID() == 2)
     return cTranspose(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::cTranspose():");
  }
  return *this;
}

Matrix& Matrix::conj(){
  try{
    if(typeID() == 1)
      return  conj(RTYPE);
    else if(typeID() == 2)
      return  conj(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::conj():");
  }
  return *this;
}

Matrix& Matrix::resize(size_t row, size_t col){
  try{
    if(typeID() == 1)
      return  resize(RTYPE, row, col);
    else if(typeID() == 2)
      return  resize(CTYPE, row, col);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::resize(size_t, size_t):");
  }
  return *this;
}

void Matrix::load(const std::string& fname){
  try{
    setMelemBNULL();
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

    if(typeID() == 1){
      init(RTYPE, ongpu);
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
    }

    else if(typeID() == 2){
      init(CTYPE, ongpu);
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
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::load(std::string&):");
  }
}

void Matrix::assign(size_t _Rnum, size_t _Cnum){
  try{
    if(typeID() == 1)
      this->assign(RTYPE, _Rnum, _Cnum);
    else if(typeID() == 2)
      this->assign(CTYPE, _Rnum, _Cnum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::assign(size_t ,size_t ):");
  }
}

bool Matrix::toGPU(){
  try{
    if(typeID() == 1)
      return toGPU(RTYPE);
    else if(typeID() == 2)
      return toGPU(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::toGPU( ):");
  }
  return false;
}

Real& Matrix::at(size_t r, size_t c){
  try{
    if(!((r < Rnum) && (c < Cnum))){
      std::ostringstream err;
      err<<"The input indices are out of range.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 1)
      m_elem = (Real*)mvCPU(m_elem, elemNum() * sizeof(Real), ongpu);
    if(diag){
      if(!(r == c && r < elemNum())){
        std::ostringstream err;
        err<<"The matrix is diagonal, there is no off-diagonal element.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      return m_elem[r];
    }
    else{
      return m_elem[r * Cnum + c];
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::at(size_t, size_t):");
  }
  return m_elem[0];
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

/*********************  developping  **********************/

Real Matrix::max(bool _ongpu){
  try{
    if(elemNum() == 0){
      std::ostringstream err;
      err<<"There is no element in this matrix ";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 2){
      std::ostringstream err;
      err<< "Can't Comparison. The type of matirx is COMPLEX." << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return this->max(RTYPE, _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::max():");
  }
  return 0.;
}

Real Matrix::absMax(bool _ongpu){
  try{
    if(elemNum() == 0){
      std::ostringstream err;
      err<<"There is no element in this matrix ";
      throw std::runtime_error(exception_msg(err.str()));
    }
    else if(typeID() == 2){
      std::ostringstream err;
      err<< "Can't Comparison. The type of matirx is COMPLEX." << std::endl <<"In the file Block.cpp, line(" << __LINE__ << ")";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return this->absMax(RTYPE, _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::absMax():");
    return 0.;
  }
}

Matrix& Matrix::normalize(){
  try{
    if(typeID() == 1)
      return this->normalize(RTYPE);
    else if(typeID() == 2)
      return this->normalize(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::normalize():");
  }
  return *this;
}

Matrix& Matrix::absMaxNorm(){
  try{
    if(typeID() == 2){
      std::ostringstream err;
      err<< "Can't perform Matrix::absMaxNorm() on this matrix. The type of matirx is COMPLEX.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return this->absMaxNorm(RTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::absMaxNorm():");
  }
  return *this;
}

Matrix& Matrix::maxNorm(){
  try{
    if(typeID() == 2){
      std::ostringstream err;
      err<< "Can't perform Matrix::maxNorm() on this matrix. The type of matirx is COMPLEX.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return this->maxNorm(RTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::maxNorm():");
  }
  return *this;
}


};	/* namespace uni10 */

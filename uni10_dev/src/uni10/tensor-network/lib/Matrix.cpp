/****************************************************************************
*  @file CMakeLists.txt
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao
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
*  @brief Main specification file for CMake
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <cmath>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
namespace uni10{
Matrix::Matrix(): Block(){}
Matrix::Matrix(const Matrix& _m): Block(_m.Rnum, _m.Cnum, _m.diag){
  try{
    if(elemNum()){
      m_elem = (double*)elemAlloc(elemNum() * sizeof(double), ongpu);
      elemCopy(m_elem, _m.m_elem, elemNum() * sizeof(double), ongpu, _m.ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(uni10::Matrix&):");
  }
}
Matrix::Matrix(const Block& _b): Block(_b.Rnum, _b.Cnum, _b.diag){
  try{
    if(elemNum()){
      m_elem = (double*)elemAlloc(elemNum() * sizeof(double), ongpu);
      elemCopy(m_elem, _b.m_elem, elemNum() * sizeof(double), ongpu, _b.ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(uni10::Matrix&):");
  }
}

void Matrix::init(bool _ongpu){
	if(elemNum()){
		if(_ongpu)	// Try to allocate GPU memory
			m_elem = (double*)elemAlloc(elemNum() * sizeof(double), ongpu);
		else{
			m_elem = (double*)elemAllocForce(elemNum() * sizeof(double), false);
			ongpu = false;
		}
	}
}

void Matrix::init(const double* _elem, bool src_ongpu){
	init(true);
	elemCopy(m_elem, _elem, elemNum() * sizeof(double), ongpu, src_ongpu);
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const double* _elem, bool _diag, bool src_ongpu): Block(_Rnum, _Cnum, _diag){
  try{
	  init(_elem, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, double*, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const std::vector<double>& _elem, bool _diag, bool src_ongpu): Block(_Rnum, _Cnum, _diag){
  try{
	  init(&_elem[0], src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, std::vector<double>&, bool=false):");
  }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu): Block(_Rnum, _Cnum, _diag){
  try{
    init(_ongpu);
    if(elemNum())
      elemBzero(m_elem, elemNum() * sizeof(double), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, bool=false):");
  }
}

Matrix& Matrix::operator=(const Matrix& _m){
  try{
    Rnum = _m.Rnum;
    Cnum = _m.Cnum;
    diag = _m.diag;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(double), ongpu);
    m_elem = (double*)elemAlloc(elemNum() * sizeof(double), ongpu);
    elemCopy(m_elem, _m.m_elem, elemNum() * sizeof(double), ongpu, _m.ongpu);
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
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(double), ongpu);
    m_elem = (double*)elemAlloc(elemNum() * sizeof(double), ongpu);
    elemCopy(m_elem, _b.m_elem, elemNum() * sizeof(double), ongpu, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator=(uni10::Block&):");
  }
	return *this;
}

Matrix::~Matrix(){
  try{
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(double), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In destructor Matrix::~Matrix():");
  }
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

void Matrix::setElem(const std::vector<double>& elem, bool _ongpu){
  try{
    setElem(&elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(std::vector<double>&, bool=false):");
  }
}
void Matrix::setElem(const double* elem, bool _ongpu){
  try{
	  elemCopy(m_elem, elem, elemNum() * sizeof(double), ongpu, _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(double*, bool=false):");
  }
}

void Matrix::randomize(){
  try{
    if(!ongpu)
      m_elem = (double*)mvGPU(m_elem, elemNum() * sizeof(double), ongpu);
    elemRand(m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::randomize():");
  }
}


void Matrix::orthoRand(){
  try{
    if(!ongpu)
      m_elem = (double*)mvGPU(m_elem, elemNum() * sizeof(double), ongpu);
    if(!diag){
      orthoRandomize(m_elem, Rnum, Cnum, ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::orthoRand():");
  }
}

void Matrix::identity(){
  try{
    diag = true;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(double), ongpu);
    m_elem = (double*)elemAlloc(elemNum() * sizeof(double), ongpu);
    double* elemI = (double*)malloc(elemNum() * sizeof(double));
    for(int i = 0; i < elemNum(); i++)
      elemI[i] = 1;
    this->setElem(elemI, false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::identity():");
  }
}

void Matrix::set_zero(){
  try{
	if(elemNum())
		elemBzero(m_elem, elemNum() * sizeof(double), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::set_zero():");
  }
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
    vectorAdd(m_elem, Mb.m_elem, elemNum(), ongpu, Mb.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator+=(uni10::Matrix&):");
  }
	return *this;
}

Matrix& Matrix::transpose(){
  try{
    if(!ongpu)
      m_elem = (double*)mvGPU(m_elem, elemNum() * sizeof(double), ongpu);
    if(!diag){
      double* transElem;
      size_t memsize = elemNum() * sizeof(double);
      transElem = (double*)elemAllocForce(memsize, ongpu);
      setTranspose(m_elem, Rnum, Cnum, transElem, ongpu);
      if(m_elem != NULL)
        elemFree(m_elem, memsize, ongpu);
      m_elem = transElem;
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
        double* elem = (double*)elemAlloc(_elemNum * sizeof(double), des_ongpu);
        elemBzero(elem, _elemNum * sizeof(double), des_ongpu);
        elemCopy(elem, m_elem, elemNum() * sizeof(double), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(double), ongpu);
        m_elem = elem;
        ongpu = des_ongpu;
      }
      else
        shrinkWithoutFree((elemNum() - _elemNum) * sizeof(double), ongpu);
      Rnum = row;
      Cnum = col;
    }
    else{
      if(col == Cnum){
        size_t _elemNum = row * col;
        if(row > Rnum){
          bool des_ongpu;
          double* elem = (double*)elemAlloc(_elemNum * sizeof(double), des_ongpu);
          elemBzero(elem, _elemNum * sizeof(double), des_ongpu);
          elemCopy(elem, m_elem, elemNum() * sizeof(double), des_ongpu, ongpu);
          if(m_elem != NULL)
            elemFree(m_elem, elemNum() * sizeof(double), ongpu);
          m_elem = elem;
          ongpu = des_ongpu;
        }
        else
          shrinkWithoutFree((elemNum() - _elemNum) * sizeof(double), ongpu);
        Rnum = row;
      }
      else{
        size_t data_row = row < Rnum ? row : Rnum;
        size_t data_col = col < Cnum ? col : Cnum;
        bool des_ongpu;
        double* elem = (double*)elemAlloc(row * col * sizeof(double), des_ongpu);
        elemBzero(elem, row * col * sizeof(double), des_ongpu);
        for(size_t r = 0; r < data_row; r++)
          elemCopy(&(elem[r * col]), &(m_elem[r * Cnum]), data_col * sizeof(double), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(double), ongpu);
        m_elem = elem;
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

void Matrix::load(const std::string& fname){
  try{
    FILE *fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    double* elem = m_elem;
    if(ongpu)
      elem = (double*)malloc(elemNum() * sizeof(double));
    fread(elem, sizeof(double), elemNum(), fp);
    fclose(fp);
    if(ongpu){
      elemCopy(m_elem, elem, elemNum() * sizeof(double), ongpu, false);
      free(elem);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::load(std::string&):");
  }
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

bool Matrix::toGPU(){
	if(!ongpu)
		m_elem = (double*)mvGPU(m_elem, elemNum() * sizeof(double), ongpu);
	return ongpu;
}

Matrix takeExp(double a, const Block& mat){
  try{
    std::vector<Matrix> rets = mat.eigh();
    Matrix UT(rets[1]);
    UT.transpose();
    vectorExp(a, rets[0].getElem(), rets[0].row(), rets[0].isOngpu());
    return UT * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function takeExp(double, uni10::Matrix&):");
    return Matrix();
  }
}
};	/* namespace uni10 */

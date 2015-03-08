/****************************************************************************
*  @file CMatrix.cpp
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
*  @brief Main specification file for CMake
*  @author Yun-Da Hsieh,Ying-Jer Kao
*  @date 2015-03-06
*  @since 1.0.0
*
*****************************************************************************/
#include <uni10/tools/uni10_tools.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/tensor-network/CMatrix.h>

typedef std::complex<double> Complex;

namespace uni10{
CMatrix::CMatrix(): CBlock(){}
CMatrix::CMatrix(const CMatrix& _m): CBlock(_m.Rnum, _m.Cnum, _m.diag){
  try{
      init(_m.m_elem, _m.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor CMatrix::CMatrix(uni10::CMatrix&):");
  }
}
CMatrix::CMatrix(const CBlock& _b): CBlock(_b){
  try{
      init(_b.m_elem, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor CMatrix::CMatrix(uni10::CBlock&):");
  }
}

void CMatrix::init(bool _ongpu){
	if(elemNum()){
		if(_ongpu)	// Try to allocate GPU memory
			m_elem = (Complex*)elemAlloc(elemNum() * sizeof(Complex), ongpu);
		else{
			m_elem = (Complex*)elemAllocForce(elemNum() * sizeof(Complex), false);
			ongpu = false;
		}
	}
}

void CMatrix::init(const Complex* _elem, bool src_ongpu){
	init(true);
	elemCopy(m_elem, _elem, elemNum() * sizeof(Complex), ongpu, src_ongpu);
}

CMatrix::CMatrix(size_t _Rnum, size_t _Cnum, const Complex* _elem, bool _diag, bool src_ongpu): CBlock(_Rnum, _Cnum, _diag){
  try{
	  init(_elem, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor CMatrix::CMatrix(size_t, size_t, std::complex<double>*, bool=false):");
  }
}

CMatrix::CMatrix(size_t _Rnum, size_t _Cnum, const std::vector<Complex>& _elem, bool _diag, bool src_ongpu): CBlock(_Rnum, _Cnum, _diag){
  try{
	  init(&_elem[0], src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor CMatrix::CMatrix(size_t, size_t, std::vector<std::<double> >&, bool=false):");
  }
}

CMatrix::CMatrix(size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu): CBlock(_Rnum, _Cnum, _diag){
  try{
    init(_ongpu);
    if(elemNum())
      elemBzero(m_elem, elemNum() * sizeof(Complex), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor CMatrix::CMatrix(size_t, size_t, bool=false):");
  }
}

CMatrix& CMatrix::operator=(const CMatrix& _m){
  try{
    Rnum = _m.Rnum;
    Cnum = _m.Cnum;
    diag = _m.diag;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(Complex), ongpu);
    init(_m.m_elem, _m.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator=(uni10::CMatrix&):");
  }
	return *this;
}
CMatrix& CMatrix::operator=(const CBlock& _b){
  try{
    Rnum = _b.Rnum;
    Cnum = _b.Cnum;
    diag = _b.diag;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(Complex), ongpu);
    init(_b.m_elem, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator=(uni10::CBlock&):");
  }
	return *this;
}

CMatrix::~CMatrix(){
  try{
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(Complex), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In destructor CMatrix::~CMatrix():");
  }
}

CMatrix& CMatrix::operator*= (const CBlock& Mb){
  try{
    if(!ongpu)
      m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    *this = *this * Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator*=(uni10::CMatrix&):");
  }
  return *this;
}

void CMatrix::setElem(const std::vector<Complex>& elem, bool _ongpu){
  try{
    setElem(&elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::setElem(std::vector< std::complex<double> >&, bool=false):");
  }
}
void CMatrix::setElem(const Complex* elem, bool _ongpu){
  try{
	  elemCopy(m_elem, elem, elemNum() * sizeof(Complex), ongpu, _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::setElem(std::complex<double>*, bool=false):");
  }
}

void CMatrix::randomize(){
  try{
    if(!ongpu)
      m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    elemRand(m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::randomize():");
  }
}


void CMatrix::orthoRand(){
  try{
    if(!ongpu)
      m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    if(!diag){
      orthoRandomize(m_elem, Rnum, Cnum, ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::orthoRand():");
  }
}

void CMatrix::identity(){
  try{
    diag = true;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(Complex), ongpu);
    m_elem = (Complex*)elemAlloc(elemNum() * sizeof(Complex), ongpu);
    Complex* elemI = (Complex*)malloc(elemNum() * sizeof(Complex));
    for(int i = 0; i < elemNum(); i++)
      elemI[i] = 1;
    this->setElem(elemI, false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::identity():");
  }
}

void CMatrix::set_zero(){
  try{
	if(elemNum())
		elemBzero(m_elem, elemNum() * sizeof(Complex), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::set_zero():");
  }
}

CMatrix& CMatrix::operator*= (double a){
  try{
    if(!ongpu)
      m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    vectorScal(a, m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator*=(double):");
  }
	return *this;
}


CMatrix& CMatrix::operator+= (const CBlock& Mb){
  try{
    if(!ongpu)
      m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    vectorAdd(m_elem, Mb.m_elem, elemNum(), ongpu, Mb.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator+=(uni10::CMatrix&):");
  }
	return *this;
}

CMatrix& CMatrix::transpose(){
  try{
    if(!ongpu)
      m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    if(!(diag || Rnum == 1 || Cnum == 1))
      setTranspose(m_elem, Rnum, Cnum, ongpu);
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::transpose():");
  }
  return *this;
}

CMatrix& CMatrix::cTranspose(){
  try{
    if(!ongpu)
      m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    if(diag || Rnum == 1 || Cnum == 1)
      this->conj();
    else
      setCTranspose(m_elem, Rnum, Cnum, ongpu);
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::cTranspose():");
  }
  return *this;
}

CMatrix& CMatrix::resize(size_t row, size_t col){
  try{
    if(diag){
      size_t _elemNum = row < col ? row : col;
      if(_elemNum > elemNum()){
        bool des_ongpu;
        Complex* elem = (Complex*)elemAlloc(_elemNum * sizeof(Complex), des_ongpu);
        elemBzero(elem, _elemNum * sizeof(Complex), des_ongpu);
        elemCopy(elem, m_elem, elemNum() * sizeof(Complex), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(Complex), ongpu);
        m_elem = elem;
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
          elemCopy(elem, m_elem, elemNum() * sizeof(Complex), des_ongpu, ongpu);
          if(m_elem != NULL)
            elemFree(m_elem, elemNum() * sizeof(Complex), ongpu);
          m_elem = elem;
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
          elemCopy(&(elem[r * col]), &(m_elem[r * Cnum]), data_col * sizeof(Complex), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(Complex), ongpu);
        m_elem = elem;
        ongpu = des_ongpu;
        Rnum = row;
        Cnum = col;
      }
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::resize(size_t, size_t):");
  }
  return *this;
}

void CMatrix::load(const std::string& fname){
  try{
    FILE *fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    Complex* elem = m_elem;
    if(ongpu)
      elem = (Complex*)malloc(elemNum() * sizeof(Complex));
    fread(elem, sizeof(Complex), elemNum(), fp);
    fclose(fp);
    if(ongpu){
      elemCopy(m_elem, elem, elemNum() * sizeof(Complex), ongpu, false);
      free(elem);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::load(std::string&):");
  }
}

Complex& CMatrix::operator[](size_t idx){
  try{
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of the matrix elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    m_elem = (Complex*)mvCPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    return m_elem[idx];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::opeartor[](size_t):");
    return m_elem[0];
  }
}

Complex* CMatrix::getHostElem(){
  try{
    if(ongpu){
      m_elem = (Complex*)mvCPU(m_elem, elemNum() * sizeof(Complex), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::getHostElem():");
  }
	return m_elem;
}


Complex& CMatrix::at(size_t r, size_t c){
  try{
    if(!((r < Rnum) && (c < Cnum))){
      std::ostringstream err;
      err<<"The input indices are out of range.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    m_elem = (Complex*)mvCPU(m_elem, elemNum() * sizeof(Complex), ongpu);
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
    propogate_exception(e, "In function CMatrix::at(size_t, size_t):");
    return m_elem[0];
  }
}

bool CMatrix::toGPU(){
	if(!ongpu)
		m_elem = (Complex*)mvGPU(m_elem, elemNum() * sizeof(Complex), ongpu);
	return ongpu;
}

CMatrix exp(double a, const CBlock& mat){
  try{
    std::vector<CMatrix> rets = mat.eig();
    CMatrix Uinv = rets[1].inverse();
    vectorExp(a, rets[0].getElem(), rets[0].row(), rets[0].isOngpu());
    return Uinv * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function exp(double, uni10::CMatrix&):");
    return CMatrix();
  }
}

CMatrix exph(double a, const CBlock& mat){
  try{
    std::vector<CMatrix> rets = mat.eigh();
    CMatrix UT(rets[1]);
    UT.cTranspose();
    vectorExp(a, rets[0].getElem(), rets[0].row(), rets[0].isOngpu());
    return UT * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function exph(double, uni10::CMatrix&):");
    return CMatrix();
  }
}

CMatrix exp(const std::complex<double>& a, const CBlock& mat){
  try{
    std::vector<CMatrix> rets = mat.eig();
    CMatrix Uinv = rets[1].inverse();
    vectorExp(a, rets[0].getElem(), rets[0].row(), rets[0].isOngpu());
    return Uinv * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function exp(std::complex<double>, uni10::CMatrix&):");
    return CMatrix();
  }
}

CMatrix exp(const CBlock& mat){
  return exp(1.0, mat);
}

CMatrix exph(const CBlock& mat){
  return exph(1.0, mat);
}

};	/* namespace uni10 */

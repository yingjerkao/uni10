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
#include <uni10/tensor-network/Matrix.h>
#include <uni10/tensor-network/CMatrix.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#ifndef UNI10_DTYPE
#define UNI10_DTYPE double
#endif
#ifndef UNI10_MATRIX
#define UNI10_MATRIX Matrix
#endif
#ifndef UNI10_BLOCK
#define UNI10_BLOCK Block
#endif
namespace uni10{
UNI10_MATRIX::UNI10_MATRIX(): UNI10_BLOCK(){}
UNI10_MATRIX::UNI10_MATRIX(const UNI10_MATRIX& _m): UNI10_BLOCK(_m.Rnum, _m.Cnum, _m.diag){
  try{
    if(elemNum()){
      m_elem = (UNI10_DTYPE*)elemAlloc(elemNum() * sizeof(UNI10_DTYPE), ongpu);
      elemCopy(m_elem, _m.m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu, _m.ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(uni10::Matrix&):");
  }
}
UNI10_MATRIX::UNI10_MATRIX(const UNI10_BLOCK& _b): UNI10_BLOCK(_b){
  try{
    if(elemNum()){
      m_elem = (UNI10_DTYPE*)elemAlloc(elemNum() * sizeof(UNI10_DTYPE), ongpu);
      elemCopy(m_elem, _b.m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu, _b.diag);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In copy constructor Matrix::Matrix(uni10::Block&):");
  }
}

void UNI10_MATRIX::init(bool _ongpu){
	if(elemNum()){
		if(_ongpu)	// Try to allocate GPU memory
			m_elem = (UNI10_DTYPE*)elemAlloc(elemNum() * sizeof(UNI10_DTYPE), ongpu);
		else{
			m_elem = (UNI10_DTYPE*)elemAllocForce(elemNum() * sizeof(UNI10_DTYPE), false);
			ongpu = false;
		}
	}
}

void UNI10_MATRIX::init(const UNI10_DTYPE* _elem, bool src_ongpu){
	init(true);
	elemCopy(m_elem, _elem, elemNum() * sizeof(UNI10_DTYPE), ongpu, src_ongpu);
}

UNI10_MATRIX::UNI10_MATRIX(size_t _Rnum, size_t _Cnum, const UNI10_DTYPE* _elem, bool _diag, bool src_ongpu): UNI10_BLOCK(_Rnum, _Cnum, _diag){
  try{
	  init(_elem, src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, double*, bool=false):");
  }
}

UNI10_MATRIX::UNI10_MATRIX(size_t _Rnum, size_t _Cnum, const std::vector<UNI10_DTYPE>& _elem, bool _diag, bool src_ongpu): UNI10_BLOCK(_Rnum, _Cnum, _diag){
  try{
	  init(&_elem[0], src_ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, std::vector<double>&, bool=false):");
  }
}

UNI10_MATRIX::UNI10_MATRIX(size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu): UNI10_BLOCK(_Rnum, _Cnum, _diag){
  try{
    init(_ongpu);
    if(elemNum())
      elemBzero(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, bool=false):");
  }
}

UNI10_MATRIX& UNI10_MATRIX::operator=(const UNI10_MATRIX& _m){
  try{
    Rnum = _m.Rnum;
    Cnum = _m.Cnum;
    diag = _m.diag;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    m_elem = (UNI10_DTYPE*)elemAlloc(elemNum() * sizeof(UNI10_DTYPE), ongpu);
    elemCopy(m_elem, _m.m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu, _m.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator=(uni10::Matrix&):");
  }
	return *this;
}
UNI10_MATRIX& UNI10_MATRIX::operator=(const UNI10_BLOCK& _b){
  try{
    Rnum = _b.Rnum;
    Cnum = _b.Cnum;
    diag = _b.diag;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    m_elem = (UNI10_DTYPE*)elemAlloc(elemNum() * sizeof(UNI10_DTYPE), ongpu);
    elemCopy(m_elem, _b.m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu, _b.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator=(uni10::Block&):");
  }
	return *this;
}

UNI10_MATRIX::~UNI10_MATRIX(){
  try{
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In destructor Matrix::~Matrix():");
  }
}

UNI10_MATRIX& UNI10_MATRIX::operator*= (const UNI10_BLOCK& Mb){
  try{
    if(!ongpu)
      m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    *this = *this * Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator*=(uni10::Matrix&):");
  }
  return *this;
}

void UNI10_MATRIX::setElem(const std::vector<UNI10_DTYPE>& elem, bool _ongpu){
  try{
    setElem(&elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(std::vector<double>&, bool=false):");
  }
}
void UNI10_MATRIX::setElem(const UNI10_DTYPE* elem, bool _ongpu){
  try{
	  elemCopy(m_elem, elem, elemNum() * sizeof(UNI10_DTYPE), ongpu, _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::setElem(double*, bool=false):");
  }
}

void UNI10_MATRIX::randomize(){
  try{
    if(!ongpu)
      m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    elemRand(m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::randomize():");
  }
}


void UNI10_MATRIX::orthoRand(){
  try{
    if(!ongpu)
      m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    if(!diag){
      orthoRandomize(m_elem, Rnum, Cnum, ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::orthoRand():");
  }
}

void UNI10_MATRIX::identity(){
  try{
    diag = true;
    if(m_elem != NULL)
      elemFree(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    m_elem = (UNI10_DTYPE*)elemAlloc(elemNum() * sizeof(UNI10_DTYPE), ongpu);
    UNI10_DTYPE* elemI = (UNI10_DTYPE*)malloc(elemNum() * sizeof(UNI10_DTYPE));
    for(int i = 0; i < elemNum(); i++)
      elemI[i] = 1;
    this->setElem(elemI, false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::identity():");
  }
}

void UNI10_MATRIX::set_zero(){
  try{
	if(elemNum())
		elemBzero(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::set_zero():");
  }
}
UNI10_MATRIX& UNI10_MATRIX::operator*= (double a){
  try{
    if(!ongpu)
      m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    vectorScal(a, m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator*=(double):");
  }
	return *this;
}

UNI10_MATRIX& UNI10_MATRIX::operator+= (const UNI10_BLOCK& Mb){
  try{
    if(!ongpu)
      m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    vectorAdd(m_elem, Mb.m_elem, elemNum(), ongpu, Mb.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::operator+=(uni10::Matrix&):");
  }
	return *this;
}

UNI10_MATRIX& UNI10_MATRIX::transpose(){
  try{
    if(!ongpu)
      m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    if(!diag)
      setTranspose(m_elem, Rnum, Cnum, ongpu);
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose():");
  }
  return *this;
}

UNI10_MATRIX& UNI10_MATRIX::cTranspose(){
  try{
    if(!ongpu)
      m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    if(!diag)
      setCTranspose(m_elem, Rnum, Cnum, ongpu);
    size_t tmp = Rnum;
    Rnum = Cnum;
    Cnum = tmp;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::transpose():");
  }
  return *this;
}

UNI10_MATRIX& UNI10_MATRIX::resize(size_t row, size_t col){
  try{
    if(diag){
      size_t _elemNum = row < col ? row : col;
      if(_elemNum > elemNum()){
        bool des_ongpu;
        UNI10_DTYPE* elem = (UNI10_DTYPE*)elemAlloc(_elemNum * sizeof(UNI10_DTYPE), des_ongpu);
        elemBzero(elem, _elemNum * sizeof(UNI10_DTYPE), des_ongpu);
        elemCopy(elem, m_elem, elemNum() * sizeof(UNI10_DTYPE), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
        m_elem = elem;
        ongpu = des_ongpu;
      }
      else
        shrinkWithoutFree((elemNum() - _elemNum) * sizeof(UNI10_DTYPE), ongpu);
      Rnum = row;
      Cnum = col;
    }
    else{
      if(col == Cnum){
        size_t _elemNum = row * col;
        if(row > Rnum){
          bool des_ongpu;
          UNI10_DTYPE* elem = (UNI10_DTYPE*)elemAlloc(_elemNum * sizeof(UNI10_DTYPE), des_ongpu);
          elemBzero(elem, _elemNum * sizeof(UNI10_DTYPE), des_ongpu);
          elemCopy(elem, m_elem, elemNum() * sizeof(UNI10_DTYPE), des_ongpu, ongpu);
          if(m_elem != NULL)
            elemFree(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
          m_elem = elem;
          ongpu = des_ongpu;
        }
        else
          shrinkWithoutFree((elemNum() - _elemNum) * sizeof(UNI10_DTYPE), ongpu);
        Rnum = row;
      }
      else{
        size_t data_row = row < Rnum ? row : Rnum;
        size_t data_col = col < Cnum ? col : Cnum;
        bool des_ongpu;
        UNI10_DTYPE* elem = (UNI10_DTYPE*)elemAlloc(row * col * sizeof(UNI10_DTYPE), des_ongpu);
        elemBzero(elem, row * col * sizeof(UNI10_DTYPE), des_ongpu);
        for(size_t r = 0; r < data_row; r++)
          elemCopy(&(elem[r * col]), &(m_elem[r * Cnum]), data_col * sizeof(UNI10_DTYPE), des_ongpu, ongpu);
        if(m_elem != NULL)
          elemFree(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
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

void UNI10_MATRIX::load(const std::string& fname){
  try{
    FILE *fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UNI10_DTYPE* elem = m_elem;
    if(ongpu)
      elem = (UNI10_DTYPE*)malloc(elemNum() * sizeof(UNI10_DTYPE));
    fread(elem, sizeof(UNI10_DTYPE), elemNum(), fp);
    fclose(fp);
    if(ongpu){
      elemCopy(m_elem, elem, elemNum() * sizeof(UNI10_DTYPE), ongpu, false);
      free(elem);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::load(std::string&):");
  }
}

UNI10_DTYPE& UNI10_MATRIX::operator[](size_t idx){
  try{
    if(!(idx < elemNum())){
      std::ostringstream err;
      err<<"Index exceeds the number of the matrix elements("<<elemNum()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    m_elem = (UNI10_DTYPE*)mvCPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    return m_elem[idx];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::opeartor[](size_t):");
    return m_elem[0];
  }
}

UNI10_DTYPE* UNI10_MATRIX::getHostElem(){
  try{
    if(ongpu){
      m_elem = (UNI10_DTYPE*)mvCPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::getHostElem():");
  }
	return m_elem;
}


UNI10_DTYPE& UNI10_MATRIX::at(size_t r, size_t c){
  try{
    if(!((r < Rnum) && (c < Cnum))){
      std::ostringstream err;
      err<<"The input indices are out of range.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    m_elem = (UNI10_DTYPE*)mvCPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
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

bool UNI10_MATRIX::toGPU(){
	if(!ongpu)
		m_elem = (UNI10_DTYPE*)mvGPU(m_elem, elemNum() * sizeof(UNI10_DTYPE), ongpu);
	return ongpu;
}

#ifndef UNI10_COMPLEX //Only for real version
UNI10_MATRIX takeExp(double a, const UNI10_BLOCK& mat){
  try{
    std::vector<UNI10_MATRIX> rets = mat.eigh();
    UNI10_MATRIX UT(rets[1]);
    UT.transpose();
    vectorExp(a, rets[0].getElem(), rets[0].row(), rets[0].isOngpu());
    return UT * (rets[0] * rets[1]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function takeExp(double, uni10::Matrix&):");
    return UNI10_MATRIX();
  }
}
#endif

#ifdef UNI10_COMPLEX //Only for complex version
UNI10_MATRIX& UNI10_MATRIX::conj(){
  setConjugate(m_elem, elemNum(), ongpu);
  return *this;
}
#endif

};	/* namespace uni10 */
#ifdef UNI10_BLOCK
#undef UNI10_BLOCK
#endif
#ifdef UNI10_MATRIX
#undef UNI10_MATRIX
#endif
#ifdef UNI10_DTYPE
#undef UNI10_DTYPE
#endif

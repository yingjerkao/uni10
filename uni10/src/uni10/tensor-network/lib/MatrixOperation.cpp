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

namespace uni10{

Matrix takeExp(double a, const Block& mat){
  try{
    return exph(a, mat);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function takeExp(double, uni10::Matrix&):");
    return Matrix();
  }
}

CMatrix& CMatrix::conj(){
  setConjugate(m_elem, elemNum(), ongpu);
  return *this;
}

Matrix::Matrix(const CBlock& _cb): Block(_cb.Rnum, _cb.Cnum, _cb.diag){
  init(RTYPE, true);
  elemCast(m_elem, _cb.m_elem, elemNum(), ongpu, _cb.ongpu);
}

CMatrix::CMatrix(const Block& _b): CBlock(_b.Rnum, _b.Cnum, _b.diag){
  init(true);
  elemCast(m_elem, _b.m_elem, elemNum(), ongpu, _b.ongpu);
}


CMatrix& CMatrix::operator*= (const std::complex<double>& a){
  try{
    if(!ongpu)
      m_elem = (std::complex<double>*)mvGPU(m_elem, elemNum() * sizeof(std::complex<double>), ongpu);
    vectorScal(a, m_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator*=(std::complex<double>):");
  }
	return *this;
}

CMatrix operator*(const CBlock& Ma, const Block& Mb){
  CMatrix cMb(Mb);
  return Ma * cMb;
}

CMatrix operator*(const Block& Ma, const CBlock& Mb){
  CMatrix cMa(Ma);
  return cMa * Mb;
}

bool operator== (const Block& m1, const CBlock& m2){
  try{
    double diff;
    if(m1.elemNum() == m2.elemNum()){
      for(size_t i = 0; i < m1.elemNum(); i++){
        if(std::abs(m2[i].imag()) > 1E-12)
          return false;
        if((std::abs(m1[i] - m2[i].real()))  > 1E-12)
          return false;
      }
    }
    else
      return false;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator==(uni10::Matrix&, uni10::CMatrix&):");
  }
  return true;
}

bool operator== (const CBlock& m1, const Block& m2){
  return (m2 == m1);
}

CMatrix operator+(const CBlock& Ma, const Block& Mb){
  try{
    CMatrix Mc(Ma);
    vectorAdd(Mc.getElem(), Mb.getElem(RTYPE), Mc.elemNum(), Mc.isOngpu(), Mb.isOngpu());
    return Mc;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator+(uni10::CMatrix&, uni10::Matrix&):");
    return CMatrix();
  }
}

CMatrix operator+(const Block& Ma, const CBlock& Mb){
  try{
    return Mb + Ma;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator+(uni10::Matrix&, uni10::CMatrix&):");
    return CMatrix();
  }
}

CMatrix& CMatrix::operator*= (const Block& Mb){
  try{
    if(!ongpu)
      m_elem = (std::complex<double>*)mvGPU(m_elem, elemNum() * sizeof(std::complex<double>), ongpu);
    *this = *this * Mb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator*=(uni10::Matrix&):");
  }
  return *this;

}

CMatrix& CMatrix::operator+= (const Block& Mb){
  try{
    if(!ongpu)
      m_elem = (std::complex<double>*)mvGPU(m_elem, elemNum() * sizeof(std::complex<double>), ongpu);
    vectorAdd(m_elem, Mb.m_elem, elemNum(), ongpu, Mb.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function CMatrix::operator+=(uni10::Matrix&):");
  }
	return *this;
}
};	/* namespace uni10 */

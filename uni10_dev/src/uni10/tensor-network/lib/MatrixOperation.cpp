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

CMatrix& CMatrix::conj(){
  setConjugate(m_elem, elemNum(), ongpu);
  return *this;
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

Matrix::Matrix(const CBlock& _cb): Block(_cb.Rnum, _cb.Cnum, _cb.diag){
  init(true);
  elemCast(m_elem, _cb.m_elem, elemNum(), ongpu, _cb.ongpu);
}

CMatrix::CMatrix(const Block& _b): CBlock(_b.Rnum, _b.Cnum, _b.diag){
  init(true);
  elemCast(m_elem, _b.m_elem, elemNum(), ongpu, _b.ongpu);
}



};	/* namespace uni10 */

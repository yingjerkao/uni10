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
#include <uni10/tools/uni10_tools.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/tensor-network/CMatrix.h>
#include <uni10/tensor-network/UniTensor.h>

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
UNI10_MATRIX otimes(const UNI10_BLOCK& Ma, const UNI10_BLOCK& Mb){
  try{
    std::vector<Bond> bonds;
    Bond bdr_a(BD_IN, Ma.row());
    Bond bdc_a(BD_OUT, Ma.col());
    Bond bdr_b(BD_IN, Mb.row());
    Bond bdc_b(BD_OUT, Mb.col());
    bonds.push_back(bdr_a); bonds.push_back(bdc_a);
    UniTensor Ta(bonds);
    bonds.clear();
    bonds.push_back(bdr_b); bonds.push_back(bdc_b);
    UniTensor Tb(bonds);
    Ta.putBlock(Ma);
    Tb.putBlock(Mb);
    return otimes(Ta, Tb).getBlock();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function otimes(uni10::Matrix&, uni10::Matrix&):");
    return UNI10_MATRIX();
  }
}
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

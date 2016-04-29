/****************************************************************************
 *  @file BlockTools.cpp
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
 *  @brief Implementation file of Block class
 *  @author Yun-Da Hsieh, Yun-Hsuan Chou
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#include <uni10/numeric/lapack/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#include <uni10/tensor-network/Matrix.h>


namespace uni10{

    
  Matrix RDotR(const Block& Ma, const Block& Mb){
    try{
      if(!(Ma.Cnum == Mb.Rnum)){
        std::ostringstream err;
        err<<"The dimensions of the two matrices do not match for matrix multiplication.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if((!Ma.diag) && (!Mb.diag)){
        Matrix Mc(RTYPE, Ma.Rnum, Mb.Cnum);
        matrixMul(Ma.m_elem, Mb.m_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.m_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
        return Mc;
      }
      else if(Ma.diag && (!Mb.diag)){
        Matrix Mc(Mb);
        Mc.resize(Ma.Rnum, Mb.Cnum);
        diagRowMul(Mc.m_elem, Ma.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
        return Mc;
      }
      else if((!Ma.diag) && Mb.diag){
        Matrix Mc(Ma);
        Mc.resize(Ma.Rnum, Mb.Cnum);
        diagColMul(Mc.m_elem, Mb.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
        return Mc;
      }
      else{
        Matrix Mc(RTYPE, Ma.Rnum, Mb.Cnum, true);
        Mc.set_zero(RTYPE);
        size_t min = std::min(Ma.elemNum(), Mb.elemNum());
        elemCopy(Mc.m_elem, Ma.m_elem, min * sizeof(Real), Mc.ongpu, Ma.ongpu);
        vectorMul(Mc.m_elem, Mb.m_elem, min, Mc.ongpu, Mb.ongpu);
        return Mc;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix RDotR(const uni10::Block&, const uni10::Block&):");
      return Block();
    }
    return Matrix();
  }

  Matrix CDotC(const Block& Ma, const Block& Mb){
    if(!(Ma.Cnum == Mb.Rnum)){
      std::ostringstream err;
      err<<"The dimensions of the two matrices do not match for matrix multiplication.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((!Ma.diag) && (!Mb.diag)){
      Matrix Mc(CTYPE, Ma.Rnum, Mb.Cnum);
      matrixMul(Ma.cm_elem, Mb.cm_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.cm_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
      return Mc;
    }
    else if(Ma.diag && (!Mb.diag)){
      Matrix Mc(Mb);
      Mc.resize(Ma.Rnum, Mb.Cnum);
      diagRowMul(Mc.cm_elem, Ma.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
      return Mc;
    }
    else if((!Ma.diag) && Mb.diag){
      Matrix Mc(Ma);
      Mc.resize(Ma.Rnum, Mb.Cnum);
      diagColMul(Mc.cm_elem, Mb.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
      return Mc;
    }
    else{
      Matrix Mc(CTYPE, Ma.Rnum, Mb.Cnum, true);
      Mc.set_zero(CTYPE);
      size_t min = std::min(Ma.elemNum(), Mb.elemNum());
      elemCopy(Mc.cm_elem, Ma.cm_elem, min * sizeof(Complex), Mc.ongpu, Ma.ongpu);
      vectorMul(Mc.cm_elem, Mb.cm_elem, min, Mc.ongpu, Mb.ongpu);
      return Mc;
    }
  }

  Matrix RDotC(const Block& Ma, const Block& Mb){

    Matrix _Ma(Ma);
    RtoC(_Ma);
    return CDotC(_Ma, Mb);

  }

  Matrix CDotR(const Block& Ma, const Block& Mb){

    Matrix _Mb(Mb);
    RtoC(_Mb);
    return CDotC(Ma, _Mb);

  }

};	/* namespace uni10 */
/*
#ifdef Block
#undef Block
#endif
#ifdef Matrix
#undef Matrix
#endif
#ifdef Real
#undef Real
#endif
*/

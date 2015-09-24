/****************************************************************************
*  @file Operation.cpp
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
*  @brief Implementation of cross UniTensor and CUniTensor methods
*  @author Yun-Da Hsieh, Ying-Jer Kao
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
#include <uni10/tensor-network/CUniTensor.h>

namespace uni10{
/*  
  UniTensor::UniTensor(const CUniTensor& UniT): name(UniT.name), status(UniT.status), bonds(UniT.bonds){
    try{
      initUniT();
      setLabel(UniT.labels);
      if(status & HAVEELEM)
        elemCast(elem, UniT.elem, m_elemNum, ongpu, UniT.ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In constructor UniTensor::UniTensor(uni10::CUniTensor& )");
    }
  }
*/
  CUniTensor::CUniTensor(const UniTensor& UniT): name(UniT.name), status(UniT.status), bonds(UniT.bonds){
    try{
      initUniT();
      setLabel(UniT.labels);
      if(status & HAVEELEM)
        elemCast(elem, UniT.elem, m_elemNum, ongpu, UniT.ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In constructor CUniTensor::CUniTensor(uni10::UniTensor& )");
    }
  }
/*
  UniTensor::UniTensor(const CBlock& blk){
    try{
      Bond bdi(BD_IN, blk.Rnum);
      Bond bdo(BD_OUT, blk.Cnum);
      bonds.push_back(bdi);
      bonds.push_back(bdo);
      initUniT();
      this->putBlock(blk);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In constructor UniTensor::UniTensor(uni10::CBlock&");
    }
  }
*/
  CUniTensor::CUniTensor(const Block& blk){
    try{
      Bond bdi(BD_IN, blk.Rnum);
      Bond bdo(BD_OUT, blk.Cnum);
      bonds.push_back(bdi);
      bonds.push_back(bdo);
      initUniT();
      putBlock(blk);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In constructor CUniTensor::CUniTensor(uni10::Block&");
    }
  }
/*
  void UniTensor::putBlock(const Qnum& qnum, const CBlock& mat){
    try{
      Matrix cm(mat);
      putBlock(qnum, cm);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function UniTensor::putBlock(uni10::Qnum&, uni10::CBlock&):");
    }
  }

  void UniTensor::putBlock(const CBlock& mat){
    try{
      Qnum q0(0);
      putBlock(q0, mat);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function UniTensor::putBlock(uni10::CBlock&):");
    }
  }
*/
  void CUniTensor::putBlock(const Qnum& qnum, const Block& mat){
    try{
      CMatrix cm(mat);
      putBlock(qnum, cm);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function CUniTensor::putBlock(uni10::Qnum&, uni10::Block&):");
    }
  }

  void CUniTensor::putBlock(const Block& mat){
    try{
      Qnum q0(0);
      putBlock(q0, mat);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function CUniTensor::putBlock(uni10::Block&):");
    }
  }
/*
  void UniTensor::setRawElem(const CBlock& blk){
    try{
      Matrix m(blk);
      setRawElem(m);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function UniTensor::setRawElem(uni10::CBlock&):");
    }
  }
*/
  void CUniTensor::setRawElem(const Block& blk){
    try{
      CMatrix m(blk);
      setRawElem(m);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function UniTensor::setRawElem(uni10::Block&):");
    }
  }

  CUniTensor& CUniTensor::operator*= (const std::complex<double>& a){
    try{
      if(!(status & HAVEELEM)){
        std::ostringstream err;
        err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      vectorScal(a, elem, m_elemNum, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function CUniTensor::operator*=(complex<double>&):");
    }
    return *this;
  }

  CUniTensor& CUniTensor::operator*= (const UniTensor& Tb){
    try{
      *this = *this * Tb;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function CUniTensor::operator*=(uni10::UniTensor&):");
    }
    return *this;
  }

  CUniTensor& CUniTensor::operator+= (const UniTensor& Tb){
    try{
      if(!(status & Tb.status & HAVEELEM)){
        std::ostringstream err;
        err<<"Cannot perform addition of tensors before setting their elements.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(!(bonds == Tb.bonds)){
        std::ostringstream err;
        err<<"Cannot perform addition of two tensors having different bonds.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      vectorAdd(elem, Tb.elem, m_elemNum, ongpu, Tb.ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function CUniTensor::operator+=(uni10::UniTensor&):");
    }
    return *this;
  }
  
  CUniTensor operator*(const CUniTensor& Ta, const UniTensor& Tb){
    try{
      CUniTensor cTa(Ta);
      CUniTensor cTb(Tb);
      return contract(cTa, cTb, true);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::CUniTensor&, uni10::UniTensor&):");
      return CUniTensor();
    }
  }

  CUniTensor operator*(const UniTensor& Ta, const CUniTensor& Tb){
    try{
      CUniTensor cTa(Ta);
      CUniTensor cTb(Tb);
      return contract(cTa, cTb, true);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::UniTensor&, uni10::CUniTensor&):");
      return CUniTensor();
    }
  }

  CUniTensor operator*(const CUniTensor& Ta, const std::complex<double>& a){
    try{
      if(!(Ta.status & Ta.HAVEELEM)){
        std::ostringstream err;
        err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      CUniTensor Tb(Ta);
      vectorScal(a, Tb.elem, Tb.m_elemNum, Tb.ongpu);
      return Tb;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::CUniTensor&, complex<double>&):");
      return CUniTensor();
    }
  }
/*  
  CUniTensor operator*(const UniTensor& Ta, const std::complex<double>& a){
    try{
      if(!(Ta.status & Ta.HAVEELEM)){
        std::ostringstream err;
        err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      CUniTensor Tb(Ta);
      vectorScal(a, Tb.elem, Tb.m_elemNum, Tb.ongpu);
      return Tb;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::UniTensor&, complex<double>&):");
      return CUniTensor();
    }
  }
*/  
  CUniTensor operator*(const std::complex<double>& a, const CUniTensor& Ta){return Ta * a;}
//  CUniTensor operator*(const std::complex<double>& a, const UniTensor& Ta){return Ta * a;}
  
  CUniTensor operator+(const CUniTensor& Ta, const UniTensor& Tb){
    CUniTensor Tc(Ta);
    return Tc += Tb;
  }
  CUniTensor operator+(const UniTensor& Ta, const CUniTensor& Tb){ return Tb + Ta;}

  CUniTensor& CUniTensor::conj(){
    setConjugate(elem, m_elemNum, ongpu);
    return *this;
  }
  CUniTensor& CUniTensor::cTranspose(){
    this->transpose();
    return this->conj();
  }


};	/* namespace uni10 */


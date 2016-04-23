/****************************************************************************
 *  @file UniTensor.cpp
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
 *  @brief Implementation of UniTensor Class
 *  @author Yun-Da Hsieh, Ying-Jer Kao
 *  @date 2015-03-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#include <uni10/tools/uni10_tools.h>
#include <uni10/numeric/lapack/uni10_lapack.h>
#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/tensor-network/UniTensor.h>
#include <uni10/hdf5io/uni10_hdf5io.h>



namespace uni10{

  void RtoC(UniTensor& UniT){
    try{
      if(UniT.typeID() == 1){
        UniT.r_flag = RNULL;
        UniT.c_flag = CTYPE;
        UniT.c_elem = (Complex*)elemAlloc( UniT.m_elemNum * sizeof(Complex), UniT.ongpu);
        elemCast(UniT.c_elem, UniT.elem, UniT.m_elemNum, UniT.ongpu, UniT.ongpu);
        UniT.initBlocks(CTYPE);
        UniT.elem = NULL;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In constructor UniTensor::RtoC( )");
    }
  }

  UniTensor contract(UniTensor& _Ta, UniTensor& _Tb, bool fast){
    try{
      if(_Ta.typeID() == 0 || _Tb.typeID() == 0){
        std::ostringstream err;
        err<<"This tensor is EMPTY ";
        throw std::runtime_error(exception_msg(err.str()));
      }else if(_Ta.typeID() == 1 && _Tb.typeID() == 1)
        return contract(RTYPE, _Ta, _Tb, fast);
      else if(_Ta.typeID() == 2 && _Tb.typeID() == 2)
        return contract(CTYPE, _Ta, _Tb, fast);
      else if(_Ta.typeID() == 1 && _Tb.typeID() == 2){
        UniTensor Ta(_Ta);
        RtoC(Ta);
        return contract(CTYPE, Ta, _Tb, fast);
      }else{
        UniTensor Tb(_Tb);
        RtoC(Tb);
        return contract(CTYPE, _Ta, Tb, fast);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function contract(uni10::UniTensor&, uni10::UniTensor, bool):");
      return UniTensor();
    }
  }

  UniTensor otimes(const UniTensor & Ta, const UniTensor& Tb){
    try{
      UniTensor T1 = Ta;
      UniTensor T2 = Tb;
      std::vector<int> label1(T1.bondNum());
      std::vector<int> label2(T2.bondNum());
      for(size_t i = 0; i < T1.bondNum(); i++){
        if(i < T1.inBondNum())
          label1[i] = i;
        else
          label1[i] = T2.inBondNum() + i;
      }
      for(size_t i = 0; i < T2.bondNum(); i++){
        if(i < T2.inBondNum())
          label2[i] = i + T1.inBondNum();
        else
          label2[i] = i + T1.bondNum();
      }
      T1.setLabel(label1);
      T2.setLabel(label2);
      return contract(T1, T2, true);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function otimes(uni10::UniTensor&, uni10::UniTensor&):");
      return UniTensor();
    }
  }

  UniTensor contract(rflag tp, UniTensor& Ta, UniTensor& Tb, bool fast){
    try{
      throwTypeError(tp);
      if(!(Ta.status & Tb.status & Ta.HAVEELEM)){
        std::ostringstream err;
        err<<"Cannot perform contraction of two tensors before setting their elements.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(&Ta == &Tb){
        UniTensor Ttmp = Tb;
        return contract(Ta, Ttmp, fast);
      }

      if(Ta.status & Ta.HAVEBOND && Tb.status & Ta.HAVEBOND){
        int AbondNum = Ta.bonds.size();
        int BbondNum = Tb.bonds.size();
        std::vector<int> oldLabelA = Ta.labels;
        std::vector<int> oldLabelB = Tb.labels;
        int oldRnumA = Ta.RBondNum;
        int oldRnumB = Tb.RBondNum;
        std::vector<int> newLabelA;
        std::vector<int> interLabel;
        std::vector<int> newLabelB;
        std::vector<int> markB(BbondNum, 0);
        std::vector<int> newLabelC;
        bool match;
        for(int a = 0; a < AbondNum; a++){
          match = false;
          for(int b = 0; b < BbondNum; b++)
            if(Ta.labels[a] == Tb.labels[b]){
              markB[b] = 1;
              interLabel.push_back(Ta.labels[a]);
              newLabelB.push_back(Tb.labels[b]);
              if(!(Ta.bonds[a].dim() == Tb.bonds[b].dim())){
                std::ostringstream err;
                err<<"Cannot contract two bonds having different dimensions";
                throw std::runtime_error(exception_msg(err.str()));
              }
              match = true;
              break;
            }
          if(!match){
            newLabelA.push_back(Ta.labels[a]);
            newLabelC.push_back(Ta.labels[a]);
          }
        }
        for(size_t a = 0; a < interLabel.size(); a++)
          newLabelA.push_back(interLabel[a]);
        for(int b = 0; b < BbondNum; b++)
          if(markB[b] == 0){
            newLabelB.push_back(Tb.labels[b]);
            newLabelC.push_back(Tb.labels[b]);
          }
        int conBond = interLabel.size();
        Ta.permute(RTYPE, newLabelA, AbondNum - conBond);
        Tb.permute(RTYPE, newLabelB, conBond);
        std::vector<Bond> cBonds;
        for(int i = 0; i < AbondNum - conBond; i++)
          cBonds.push_back(Ta.bonds[i]);
        for(int i = conBond; i < BbondNum; i++)
          cBonds.push_back(Tb.bonds[i]);
        UniTensor Tc(RTYPE, cBonds);
        if(cBonds.size())
          Tc.setLabel(newLabelC);
        Block blockA, blockB, blockC;
        std::map<Qnum, Block>::iterator it;
        std::map<Qnum, Block>::iterator it2;
        for(it = Ta.blocks.begin() ; it != Ta.blocks.end(); it++){
          if((it2 = Tb.blocks.find(it->first)) != Tb.blocks.end()){
            blockA = it->second;
            blockB = it2->second;
            blockC = Tc.blocks[it->first];
            if(!(blockA.row() == blockC.row() && blockB.col() == blockC.col() && blockA.col() == blockB.row())){
              std::ostringstream err;
              err<<"The dimensions the bonds to be contracted out are different.";
              throw std::runtime_error(exception_msg(err.str()));
            }
            matrixMul(blockA.getElem(RTYPE), blockB.getElem(RTYPE), blockA.row(), blockB.col(), blockA.col(), blockC.getElem(RTYPE), Ta.ongpu, Tb.ongpu, Tc.ongpu);
          }
        }
        Tc.status |= Tc.HAVEELEM;

        if(conBond == 0){	//Outer product
          int idx = 0;
          for(int i = 0; i < oldRnumA; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(int i = 0; i < oldRnumB; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          for(int i = oldRnumA; i < AbondNum; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(int i = oldRnumB; i < BbondNum; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          Tc.permute(newLabelC, oldRnumA + oldRnumB);
        }

        if(!fast){
          Ta.permute(RTYPE, oldLabelA, oldRnumA);
          Tb.permute(RTYPE, oldLabelB, oldRnumB);
        }
        return Tc;
      }
      else if(Ta.status & Ta.HAVEBOND)
        return Ta * Tb.at(RTYPE, 0);
      else if(Tb.status & Tb.HAVEBOND)
        return Ta.at(RTYPE, 0) * Tb;
      else
        return UniTensor(Ta.at(RTYPE, 0) * Tb.at(RTYPE, 0));
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function contract(uni10::UniTensor&, uni10::UniTensor, bool):");
      return UniTensor();
    }
  }

  UniTensor otimes(rflag tp, const UniTensor & Ta, const UniTensor& Tb){
    try{
      throwTypeError(tp);
      UniTensor T1 = Ta;
      UniTensor T2 = Tb;
      std::vector<int> label1(T1.bondNum());
      std::vector<int> label2(T2.bondNum());
      for(size_t i = 0; i < T1.bondNum(); i++){
        if(i < T1.inBondNum())
          label1[i] = i;
        else
          label1[i] = T2.inBondNum() + i;
      }
      for(size_t i = 0; i < T2.bondNum(); i++){
        if(i < T2.inBondNum())
          label2[i] = i + T1.inBondNum();
        else
          label2[i] = i + T1.bondNum();
      }
      T1.setLabel(label1);
      T2.setLabel(label2);
      return contract(RTYPE, T1, T2, true);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function otimes(uni10::UniTensor&, uni10::UniTensor&):");
      return UniTensor();
    }
  }

  UniTensor contract(cflag tp, UniTensor& Ta, UniTensor& Tb, bool fast){
    try{
      throwTypeError(tp);
      if(!(Ta.status & Tb.status & Ta.HAVEELEM)){
        std::ostringstream err;
        err<<"Cannot perform contraction of two tensors before setting their elements.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(&Ta == &Tb){
        UniTensor Ttmp = Tb;
        return contract(Ta, Ttmp, fast);
      }

      if(Ta.status & Ta.HAVEBOND && Tb.status & Ta.HAVEBOND){
        int AbondNum = Ta.bonds.size();
        int BbondNum = Tb.bonds.size();
        std::vector<int> oldLabelA = Ta.labels;
        std::vector<int> oldLabelB = Tb.labels;
        int oldRnumA = Ta.RBondNum;
        int oldRnumB = Tb.RBondNum;
        std::vector<int> newLabelA;
        std::vector<int> interLabel;
        std::vector<int> newLabelB;
        std::vector<int> markB(BbondNum, 0);
        std::vector<int> newLabelC;
        bool match;
        for(int a = 0; a < AbondNum; a++){
          match = false;
          for(int b = 0; b < BbondNum; b++)
            if(Ta.labels[a] == Tb.labels[b]){
              markB[b] = 1;
              interLabel.push_back(Ta.labels[a]);
              newLabelB.push_back(Tb.labels[b]);
              if(!(Ta.bonds[a].dim() == Tb.bonds[b].dim())){
                std::ostringstream err;
                err<<"Cannot contract two bonds having different dimensions";
                throw std::runtime_error(exception_msg(err.str()));
              }
              match = true;
              break;
            }
          if(!match){
            newLabelA.push_back(Ta.labels[a]);
            newLabelC.push_back(Ta.labels[a]);
          }
        }
        for(size_t a = 0; a < interLabel.size(); a++)
          newLabelA.push_back(interLabel[a]);
        for(int b = 0; b < BbondNum; b++)
          if(markB[b] == 0){
            newLabelB.push_back(Tb.labels[b]);
            newLabelC.push_back(Tb.labels[b]);
          }
        int conBond = interLabel.size();
        Ta.permute(CTYPE, newLabelA, AbondNum - conBond);
        Tb.permute(CTYPE, newLabelB, conBond);
        std::vector<Bond> cBonds;
        for(int i = 0; i < AbondNum - conBond; i++)
          cBonds.push_back(Ta.bonds[i]);
        for(int i = conBond; i < BbondNum; i++)
          cBonds.push_back(Tb.bonds[i]);
        UniTensor Tc(CTYPE, cBonds);
        if(cBonds.size())
          Tc.setLabel(newLabelC);
        Block blockA, blockB, blockC;
        std::map<Qnum, Block>::iterator it;
        std::map<Qnum, Block>::iterator it2;
        for(it = Ta.blocks.begin() ; it != Ta.blocks.end(); it++){
          if((it2 = Tb.blocks.find(it->first)) != Tb.blocks.end()){
            blockA = it->second;
            blockB = it2->second;
            blockC = Tc.blocks[it->first];
            if(!(blockA.row() == blockC.row() && blockB.col() == blockC.col() && blockA.col() == blockB.row())){
              std::ostringstream err;
              err<<"The dimensions the bonds to be contracted out are different.";
              throw std::runtime_error(exception_msg(err.str()));
            }
            matrixMul(blockA.getElem(CTYPE), blockB.getElem(CTYPE), blockA.row(), blockB.col(), blockA.col(), blockC.getElem(CTYPE), Ta.ongpu, Tb.ongpu, Tc.ongpu);
          }
        }
        Tc.status |= Tc.HAVEELEM;

        if(conBond == 0){	//Outer product
          int idx = 0;
          for(int i = 0; i < oldRnumA; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(int i = 0; i < oldRnumB; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          for(int i = oldRnumA; i < AbondNum; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(int i = oldRnumB; i < BbondNum; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          Tc.permute(newLabelC, oldRnumA + oldRnumB);
        }

        if(!fast){
          Ta.permute(CTYPE, oldLabelA, oldRnumA);
          Tb.permute(CTYPE, oldLabelB, oldRnumB);
        }
        return Tc;
      }
      else if(Ta.status & Ta.HAVEBOND)
        return Ta * Tb.at(CTYPE, 0);
      else if(Tb.status & Tb.HAVEBOND)
        return Ta.at(CTYPE, 0) * Tb;
      else
        return UniTensor(Ta.at(CTYPE, 0) * Tb.at(CTYPE, 0));
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function contract(uni10::UniTensor&, uni10::UniTensor, bool):");
      return UniTensor();
    }
  }

  UniTensor otimes(cflag tp, const UniTensor & Ta, const UniTensor& Tb){
    try{
      throwTypeError(tp);
      UniTensor T1 = Ta;
      UniTensor T2 = Tb;
      std::vector<int> label1(T1.bondNum());
      std::vector<int> label2(T2.bondNum());
      for(size_t i = 0; i < T1.bondNum(); i++){
        if(i < T1.inBondNum())
          label1[i] = i;
        else
          label1[i] = T2.inBondNum() + i;
      }
      for(size_t i = 0; i < T2.bondNum(); i++){
        if(i < T2.inBondNum())
          label2[i] = i + T1.inBondNum();
        else
          label2[i] = i + T1.bondNum();
      }
      T1.setLabel(label1);
      T2.setLabel(label2);
      return contract(CTYPE, T1, T2, true);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function otimes(uni10::UniTensor&, uni10::UniTensor&):");
      return UniTensor();
    }
  }

}; /* namespace uni10 */

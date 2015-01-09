/****************************************************************************
*  @file CMakeLists.txt
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
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/tensor-network/UniTensor.h>
#include <uni10/tools/uni10_tools.h>
#include <uni10/numeric/uni10_lapack.h>
//using namespace uni10::datatype;
namespace uni10 {
UniTensor& UniTensor::operator+= (const UniTensor& Tb) {
    try {
        if(!(status & Tb.status & HAVEELEM)) {
            std::ostringstream err;
            err<<"Cannot perform addition of tensors before setting their elements.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if(!(bonds == Tb.bonds)) {
            std::ostringstream err;
            err<<"Cannot perform addition of two tensors having different bonds.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        vectorAdd(elem, Tb.elem, m_elemNum, ongpu, Tb.ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function UniTensor::operator+=(uni10::UniTensor&):");
    }
    return *this;
}
UniTensor operator+(const UniTensor& Ta, const UniTensor& Tb) {
    try {
        if(!(Ta.status & Tb.status & Ta.HAVEELEM)) {
            std::ostringstream err;
            err<<"Cannot perform addition of tensors before setting their elements.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if(!(Ta.bonds == Tb.bonds)) {
            std::ostringstream err;
            err<<"Cannot perform addition of two tensors having diffirent bonds.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        UniTensor Tc(Ta);
        vectorAdd(Tc.elem, Tb.elem, Tc.m_elemNum, Tc.ongpu, Tb.ongpu);
        return Tc;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator+(uni10::UniTensor&, uni10::UniTensor&):");
        return UniTensor();
    }
}

UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb) {
    try {
        UniTensor cTa = Ta;
        UniTensor cTb = Tb;
        return contract(cTa, cTb, true);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator*(uni10::UniTensor&, uni10::UniTensor&):");
        return UniTensor();
    }
}

UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast) {
    try {
        if(!(Ta.status & Tb.status & Ta.HAVEELEM)) {
            std::ostringstream err;
            err<<"Cannot perform contraction of two tensors before setting their elements.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if(&Ta == &Tb) {
            UniTensor Ttmp = Tb;
            return contract(Ta, Ttmp, fast);
        }
        if(Ta.status & Ta.HAVEBOND && Tb.status & Ta.HAVEBOND) {
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
            for(int a = 0; a < AbondNum; a++) {
                match = false;
                for(int b = 0; b < BbondNum; b++)
                    if(Ta.labels[a] == Tb.labels[b]) {
                        markB[b] = 1;
                        interLabel.push_back(Ta.labels[a]);
                        newLabelB.push_back(Tb.labels[b]);
                        if(!(Ta.bonds[a].dim() == Tb.bonds[b].dim())) {
                            std::ostringstream err;
                            err<<"Cannot contract two bonds having different dimensions";
                            throw std::runtime_error(exception_msg(err.str()));
                        }
                        match = true;
                        break;
                    }
                if(!match) {
                    newLabelA.push_back(Ta.labels[a]);
                    newLabelC.push_back(Ta.labels[a]);
                }
            }
            for(int a = 0; a < interLabel.size(); a++)
                newLabelA.push_back(interLabel[a]);
            for(int b = 0; b < BbondNum; b++)
                if(markB[b] == 0) {
                    newLabelB.push_back(Tb.labels[b]);
                    newLabelC.push_back(Tb.labels[b]);
                }
            int conBond = interLabel.size();
            Ta.permute(newLabelA, AbondNum - conBond);
            Tb.permute(newLabelB, conBond);
            std::vector<Bond> cBonds;
            for(int i = 0; i < AbondNum - conBond; i++)
                cBonds.push_back(Ta.bonds[i]);
            for(int i = conBond; i < BbondNum; i++)
                cBonds.push_back(Tb.bonds[i]);
            UniTensor Tc(cBonds);
            if(cBonds.size())
                Tc.setLabel(newLabelC);
            Block blockA, blockB, blockC;
            std::map<Qnum,Block>::iterator it;
            std::map<Qnum,Block>::iterator it2;
            for(it = Ta.blocks.begin() ; it != Ta.blocks.end(); it++) {
                if((it2 = Tb.blocks.find(it->first)) != Tb.blocks.end()) {
                    blockA = it->second;
                    blockB = it2->second;
                    blockC = Tc.blocks[it->first];
                    if(!(blockA.Rnum == blockC.Rnum && blockB.Cnum == blockC.Cnum && blockA.Cnum == blockB.Rnum)) {
                        std::ostringstream err;
                        err<<"The dimensions the bonds to be contracted out are different.";
                        throw std::runtime_error(exception_msg(err.str()));
                    }
                    matrixMul(blockA.elem, blockB.elem, blockA.Rnum, blockB.Cnum, blockA.Cnum, blockC.elem, Ta.ongpu, Tb.ongpu, Tc.ongpu);
                }
            }
            Tc.status |= Tc.HAVEELEM;

            if(conBond == 0) { //Outer product
                int idx = 0;
                for(int i = 0; i < oldRnumA; i++) {
                    newLabelC[idx] = oldLabelA[i];
                    idx++;
                }
                for(int i = 0; i < oldRnumB; i++) {
                    newLabelC[idx] = oldLabelB[i];
                    idx++;
                }
                for(int i = oldRnumA; i < AbondNum; i++) {
                    newLabelC[idx] = oldLabelA[i];
                    idx++;
                }
                for(int i = oldRnumB; i < BbondNum; i++) {
                    newLabelC[idx] = oldLabelB[i];
                    idx++;
                }
                Tc.permute(newLabelC, oldRnumA + oldRnumB);
            }

            if(!fast) {
                Ta.permute(oldLabelA, oldRnumA);
                Tb.permute(oldLabelB, oldRnumB);
            }
            return Tc;
        }
        else if(Ta.status & Ta.HAVEBOND)
            return Ta * Tb[0];
        else if(Tb.status & Tb.HAVEBOND)
            return Ta[0] * Tb;
        else
            return UniTensor(Ta[0] * Tb[0]);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function contract(uni10::UniTensor&, uni10::UniTensor, bool):");
        return UniTensor();
    }
}

UniTensor& UniTensor::operator*=(const UniTensor& uT) {
    try {
        *this = *this * uT;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function UniTensor::operator*=(uni10::UniTensor&):");
    }
    return *this;
}

UniTensor& UniTensor::operator*= (double a) {
    try {
        if(!(status & HAVEELEM)) {
            std::ostringstream err;
            err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        vectorScal(a, elem, m_elemNum, ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function UniTensor::operator*=(double):");
    }
    return *this;
}

UniTensor operator*(const UniTensor& Ta, double a) {
    try {
        if(!(Ta.status & Ta.HAVEELEM)) {
            std::ostringstream err;
            err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        UniTensor Tb(Ta);
        vectorScal(a, Tb.elem, Tb.m_elemNum, Tb.ongpu);
        return Tb;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator*(uni10::UniTensor&, double):");
        return UniTensor();
    }
}
UniTensor otimes(const UniTensor & Ta, const UniTensor& Tb) {
    try {
        UniTensor T1 = Ta;
        UniTensor T2 = Tb;
        std::vector<int> label1(T1.bondNum());
        std::vector<int> label2(T2.bondNum());
        for(int i = 0; i < T1.bondNum(); i++) {
            if(i < T1.inBondNum())
                label1[i] = i;
            else
                label1[i] = T2.inBondNum() + i;
        }
        for(int i = 0; i < T2.bondNum(); i++) {
            if(i < T2.inBondNum())
                label2[i] = i + T1.inBondNum();
            else
                label2[i] = i + T1.bondNum();
        }
        T1.setLabel(label1);
        T2.setLabel(label2);
        return contract(T1, T2, true);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function otimes(uni10::UniTensor&, uni10::UniTensor&):");
        return UniTensor();
    }
}
Matrix otimes(const Matrix& Ma, const Matrix& Mb) {
    try {
        std::vector<Bond> bonds;
        Bond bdr_a(BD_IN, Ma.row());
        Bond bdc_a(BD_OUT, Ma.col());
        Bond bdr_b(BD_IN, Mb.row());
        Bond bdc_b(BD_OUT, Mb.col());
        bonds.push_back(bdr_a);
        bonds.push_back(bdc_a);
        UniTensor Ta(bonds);
        bonds.clear();
        bonds.push_back(bdr_b);
        bonds.push_back(bdc_b);
        UniTensor Tb(bonds);
        Ta.putBlock(Ma);
        Tb.putBlock(Mb);
        return otimes(Ta, Tb).getBlock();
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function otimes(uni10::Matrix&, uni10::Matrix&):");
        return Matrix();
    }

}
};  /* namespace uni10 */

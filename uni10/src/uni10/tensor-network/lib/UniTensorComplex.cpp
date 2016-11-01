/****************************************************************************
*  @file UniTensorComplex.cpp
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


namespace uni10{

UniTensor::UniTensor(Complex val): status(0){ //GPU
  try{
    initUniT(CTYPE);
    setElemAt(0, val, c_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(Complex val):");
  }
}

UniTensor::UniTensor(cflag _tp, const std::vector<Bond>& _bonds, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    throwTypeError(_tp);
    initUniT(_tp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(uni10::cflag, std::vector<Bond>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(cflag _tp, const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    throwTypeError(_tp);
    initUniT(_tp);
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(uni10::cflag, std::vector<Bond>&, std::vector<int>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(cflag _tp, const std::vector<Bond>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    throwTypeError(_tp);
    initUniT(_tp);
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(uni10::cflag, std::vector<Bond>&, int*, std::string& = \"\"):");
  }
}

void UniTensor::setRawElem(const std::vector<Complex>& rawElem){
  try{
    setRawElem(&rawElem[0]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(std::vector<std::complex<double>>&):");
  }
}

void UniTensor::setRawElem(cflag tp, const Block& blk){
  try{
    throwTypeError(tp);
    setRawElem(blk.getElem(CTYPE));
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(uni10::cflag, uni10::Block&):");
  }
}

void UniTensor::setRawElem(const Complex* rawElem){
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"Setting elements to a tensor without bonds is not supported.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 1)
      this->assign(CTYPE, this->bond());

    int bondNum = bonds.size();
    std::vector<int> Q_idxs(bondNum, 0);
    std::vector<int> Q_Bdims(bondNum, 0);
    std::vector<int> sB_idxs(bondNum, 0);
    std::vector<int> sB_sBdims(bondNum, 0);
    std::vector<int> rAcc(bondNum, 1);
    for(int b = 0; b < bondNum; b++)
      Q_Bdims[b] = bonds[b].Qnums.size();
    for(int b = bondNum - 1; b > 0; b--)
      rAcc[b - 1] = rAcc[b] * bonds[b].dim();
    int Q_off;
    int tmp;
    int RQoff, CQoff;
    size_t sB_r, sB_c;	//sub-block of a Qidx
    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
    size_t B_cDim;
    size_t E_off;
    int R_off;
    Complex* work = c_elem;
    if(ongpu){
      work = (Complex*)malloc(m_elemNum * sizeof(Complex));
    }
    for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
      Q_off = it->first;
      tmp = Q_off;
      for(int b = bondNum - 1; b >= 0; b--){
        Q_idxs[b] = tmp % Q_Bdims[b];
        tmp /= Q_Bdims[b];
      }
      R_off = 0;
      for(int b = 0; b < bondNum; b++){
        R_off += rAcc[b] * bonds[b].offsets[Q_idxs[b]];
        sB_sBdims[b] = bonds[b].Qdegs[Q_idxs[b]];
      }
      RQoff = Q_off / CQdim;
      CQoff = Q_off % CQdim;
      B_cDim = RQidx2Blk[RQoff]->Cnum;
      E_off = (RQidx2Blk[RQoff]->cm_elem - c_elem) + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
      sB_rDim = RQidx2Dim[RQoff];
      sB_cDim = CQidx2Dim[CQoff];
      sB_idxs.assign(bondNum, 0);
      for(sB_r = 0; sB_r < sB_rDim; sB_r++)
        for(sB_c = 0; sB_c < sB_cDim; sB_c++){
          work[E_off + (sB_r * B_cDim) + sB_c] = rawElem[R_off];
          for(int bend = bondNum - 1; bend >= 0; bend--){
            sB_idxs[bend]++;
            if(sB_idxs[bend] < sB_sBdims[bend]){
              R_off += rAcc[bend];
              break;
            }
            else{
              R_off -= rAcc[bend] * (sB_idxs[bend] - 1);
              sB_idxs[bend] = 0;
            }
          }
        }
    }
    if(ongpu){
      elemCopy(c_elem, work, m_elemNum * sizeof(Complex), ongpu, false);
      free(work);
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(uni10::cflag, std::complex<double>*):");
  }
}

void UniTensor::setElem(const Complex* _elem, bool _ongpu){
  try{
    if(typeID() == 1)
      this->assign(CTYPE, this->bond());
    elemCopy(c_elem, _elem, m_elemNum * sizeof(Complex), ongpu, _ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(std::complex<double>*, bool=false):");
  }
}

void UniTensor::setElem(const std::vector<Complex>& _elem, bool _ongpu){
  try{
    setElem(&_elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(std::vector<double>&, bool=false):");
  }
}

void UniTensor::putBlock(cflag uni10_tp, const Block& mat, bool force){

  try{

    //checkUni10TypeError(uni10_tp);
    throwTypeError(uni10_tp);

    Qnum q0(0);

    putBlock(CTYPE, q0, mat, force);

  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::cflag, uni10::Block&):");
  }

}

void UniTensor::putBlock(cflag uni10_tp, const Qnum& qnum, const Block& mat, bool force){

  try{

    //checkUni10TypeError(uni10_tp);
    throwTypeError(uni10_tp);

    std::map<Qnum, Block>::iterator it;

    if( !force && mat.typeID() == 1){
      std::ostringstream err;
      err<<"\n1. Can not put a Real(RTYPE) Matrix into a Complex(CTYPE) UniTensor\n\n2. Or you can turn on the force flag, UniTensor::putBlock(qnum, mat, true) or UniTensor::putBlock(CTYPE, qnum, mat, true). \n";
      throw std::runtime_error(exception_msg(err.str()));
    }
    
    if(!((it = blocks.find(qnum)) != blocks.end())){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }

    if(!(mat.row() == it->second.Rnum && mat.col() == it->second.Cnum)){
      std::ostringstream err;
      err<<"The dimension of input matrix does not match for the dimension of the block with quantum number "<<qnum<<std::endl;
      err<<"  Hint: Use Matrix::resize(int, int)";
      throw std::runtime_error(exception_msg(err.str()));
    }

    Matrix tmp = mat;

    if( force && mat.typeID() == 1)
      RtoC(tmp);
      
    if(tmp.cm_elem != it->second.cm_elem){

      if(tmp.isDiag()){

        elemBzero(it->second.cm_elem, it->second.Rnum * it->second.Cnum * sizeof(Complex), ongpu);

        setDiag(it->second.cm_elem,tmp.getElem(CTYPE), it->second.Rnum, it->second.Cnum, tmp.elemNum(), ongpu, tmp.isOngpu());

      }
      else
        elemCopy(it->second.cm_elem, tmp.getElem(CTYPE), it->second.Rnum * it->second.Cnum * sizeof(Complex), ongpu, tmp.isOngpu());
    }

    status |= HAVEELEM;

  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::cflag, uni10::Qnum&, uni10::Block&):");
  }

}

Matrix UniTensor::getRawElem(cflag tp)const{
  try{
    throwTypeError(tp);
    if(status & HAVEBOND && status & HAVEELEM){
      int bondNum = bonds.size();
      size_t rowNum = 1;
      size_t colNum = 1;
      for(std::vector<Bond>::const_iterator it = bonds.begin(); it != bonds.end(); ++it){
        if(it->type() == BD_IN)
          rowNum *= it->dim();
        else
          colNum *= it->dim();
      }
      std::vector<size_t> idxs(bondNum, 0);
      int bend;
      std::vector<Complex> rawElem;
      while(1){
        rawElem.push_back(at(CTYPE, idxs));
        for(bend = bondNum - 1; bend >= 0; bend--){
          idxs[bend]++;
          if(idxs[bend] < bonds[bend].dim())
            break;
          else
            idxs[bend] = 0;
        }
        if(bend < 0)
          break;
      }
      return Matrix(rowNum, colNum, &rawElem[0]);
    }
    else if(status & HAVEELEM)
      return Matrix(CTYPE, 1, 1, c_elem);
    else
      return Matrix();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getRawElem(uni10::cflag ):");
    return Matrix();
  }
}

Complex* UniTensor::getElem(cflag tp){
  try{
    throwTypeError(tp);
    if(typeID() == 1){
      std::ostringstream err;
      err<<"This Tensor is REAL. Please use UniTensor::getElem(uni10::rflag ) instead";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getElem(uni10::cflag ):");
  }
  return c_elem;
}

std::map<Qnum, Matrix> UniTensor::getBlocks(cflag tp)const{
  std::map<Qnum, Matrix> mats;
  try{
    throwTypeError(tp);
    for(std::map<Qnum, Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.cm_elem, false, ongpu);
      mats.insert(std::pair<Qnum, Matrix>(it->first, mat));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlocks(uni10::cflag ):");
  }
  return mats;
}

Matrix UniTensor::getBlock(cflag tp, bool diag)const{
  try{
    throwTypeError(tp);
    Qnum q0(0);
    return getBlock(CTYPE, q0, diag);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(uni10::cflag, bool=false):");
    return Matrix();
  }
}

Matrix UniTensor::getBlock(cflag tp, const Qnum& qnum, bool diag)const{
  try{
    throwTypeError(tp);
    std::map<Qnum, Block>::const_iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag)
      return it->second.getDiag();
    else{
      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.cm_elem, false, ongpu);
      return mat;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(uni10::cflag, uni10::Qnum&):");
    return Matrix();
  }
}

void UniTensor::set_zero(cflag tp){
  try{
    throwTypeError(tp);
    elemBzero(c_elem, m_elemNum * sizeof(Complex), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero(uni10::cflag ):");
  }
}

void UniTensor::set_zero(cflag tp, const Qnum& qnum){
  try{
    throwTypeError(tp);
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    elemBzero(block.cm_elem, block.Rnum * block.Cnum * sizeof(Complex), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero(uni10::cflag, std::Qnum&):");
  }
}

void UniTensor::identity(cflag tp){
  try{
    throwTypeError(tp);
    std::map<Qnum, Block>::iterator it;
    for ( it = blocks.begin() ; it != blocks.end(); it++ )
      setIdentity(it->second.cm_elem, it->second.Rnum, it->second.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity(uni10::cflag ):");
  }
}

void UniTensor::identity(cflag tp, const Qnum& qnum){
  try{
    throwTypeError(tp);
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    setIdentity(block.cm_elem, block.Rnum, block.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity(uni10::cflag, std::Qnum&):");
  }
}

void UniTensor::randomize(cflag tp){
  try{
    throwTypeError(tp);
    elemRand(c_elem, m_elemNum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::randomize(uni10::cflag ):");
  }
}

void UniTensor::orthoRand(cflag tp){
  try{
    throwTypeError(tp);
    std::map<Qnum, Block>::iterator it;
    for ( it = blocks.begin() ; it != blocks.end(); it++ )
      orthoRandomize(it->second.cm_elem, it->second.Rnum, it->second.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand(uni10::cflag ):");
  }
}

void UniTensor::orthoRand(cflag tp, const Qnum& qnum){
  try{
    throwTypeError(tp);
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    orthoRandomize(block.cm_elem, block.Rnum, block.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand(uni10::cflag, std::Qnum&):");
  }
}

UniTensor& UniTensor::transpose(cflag tp){
  try{
    throwTypeError(tp);
    if(!(status & HAVEBOND)){
      std::ostringstream err;
      err<<"There is no bond in the tensor(scalar) to perform transposition.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<int> rsp_outin(bondNum);
    int rbondNum = 0;
    for(int b = 0; b < bondNum; b++)
      if(bonds[b].type() == BD_IN)
        rbondNum++;
      else
        break;
    int cbondNum = bondNum - rbondNum;
    for(int b = 0; b < bondNum; b++)
      if(b < cbondNum)
        rsp_outin[b] = rbondNum + b;
      else
        rsp_outin[b] = b - cbondNum;
    std::vector<int> outLabels(bondNum, 0);
    std::vector<Bond> outBonds;
    for(size_t b = 0; b < bonds.size(); b++){
      outBonds.push_back(bonds[rsp_outin[b]]);
      outLabels[b] = labels[rsp_outin[b]];
    }
    for(int b = 0; b < bondNum; b++){
      if(b < cbondNum)
        outBonds[b].m_type = BD_IN;
      else
        outBonds[b].m_type = BD_OUT;
    }
    UniTensor UniTout(CTYPE, outBonds, name);
    UniTout.setLabel(outLabels);
    if(status & HAVEELEM){
      std::map<Qnum, Block>::iterator it_in;
      std::map<Qnum, Block>::iterator it_out;
      Complex* elem_in;
      Complex* elem_out;
      size_t Rnum, Cnum;
      for ( it_in = blocks.begin() ; it_in != blocks.end(); it_in++ ){
        it_out = UniTout.blocks.find((it_in->first));
        Rnum = it_in->second.Rnum;
        Cnum = it_in->second.Cnum;
        elem_in = it_in->second.cm_elem;
        elem_out = it_out->second.cm_elem;
        setTranspose(elem_in, Rnum, Cnum, elem_out, ongpu, UniTout.ongpu);
      }
      UniTout.status |= HAVEELEM;
    }
    *this = UniTout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::transpose(uni10::cflag ):");
  }
  return *this;
}

UniTensor& UniTensor::cTranspose(cflag tp){
  try{
    throwTypeError(tp);
    if(!(status & HAVEBOND)){
      std::ostringstream err;
      err<<"There is no bond in the tensor(scalar) to perform transposition.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<int> rsp_outin(bondNum);
    int rbondNum = 0;
    for(int b = 0; b < bondNum; b++)
      if(bonds[b].type() == BD_IN)
        rbondNum++;
      else
        break;
    int cbondNum = bondNum - rbondNum;
    for(int b = 0; b < bondNum; b++)
      if(b < cbondNum)
        rsp_outin[b] = rbondNum + b;
      else
        rsp_outin[b] = b - cbondNum;
    std::vector<int> outLabels(bondNum, 0);
    std::vector<Bond> outBonds;
    for(size_t b = 0; b < bonds.size(); b++){
      outBonds.push_back(bonds[rsp_outin[b]]);
      outLabels[b] = labels[rsp_outin[b]];
    }
    for(int b = 0; b < bondNum; b++){
      if(b < cbondNum)
        outBonds[b].m_type = BD_IN;
      else
        outBonds[b].m_type = BD_OUT;
    }
    UniTensor UniTout(CTYPE, outBonds, name);
    UniTout.setLabel(outLabels);
    if(status & HAVEELEM){
      std::map<Qnum, Block>::iterator it_in;
      std::map<Qnum, Block>::iterator it_out;
      Complex* elem_in;
      Complex* elem_out;
      size_t Rnum, Cnum;
      for ( it_in = blocks.begin() ; it_in != blocks.end(); it_in++ ){
        it_out = UniTout.blocks.find((it_in->first));
        Rnum = it_in->second.Rnum;
        Cnum = it_in->second.Cnum;
        elem_in = it_in->second.cm_elem;
        elem_out = it_out->second.cm_elem;
        setCTranspose(elem_in, Rnum, Cnum, elem_out, ongpu, UniTout.ongpu);
      }
      UniTout.status |= HAVEELEM;
    }
    *this = UniTout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::cTranspose(uni10::cflag ):");
  }
  return *this;

}

UniTensor& UniTensor::permute(cflag tp, int rowBondNum){
  try{
    throwTypeError(tp);
    std::vector<int> ori_labels = labels;
    this->permute(CTYPE, ori_labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(uni10::cflag, int):");
  }
  return *this;
}

UniTensor& UniTensor::permute(cflag tp, int* newLabels, int rowBondNum){
  try{
    throwTypeError(tp);
    std::vector<int> _labels(newLabels, newLabels + bonds.size());
    this->permute(CTYPE, _labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(uni10::clfag, int*, int):");
  }
  return *this;
}

UniTensor& UniTensor::permute(cflag tp, const std::vector<int>& newLabels, int rowBondNum){
  try{
    throwTypeError(tp);
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"There is no bond in the tensor(scalar) to permute.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((labels.size() == newLabels.size()) == 0){
      std::ostringstream err;
      err<<"The size of the input new labels does not match for the number of bonds.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<int> rsp_outin(bondNum);
    int cnt = 0;
    for(int i = 0; i < bondNum; i++)
      for(int j = 0; j < bondNum; j++)
        if(labels[i] == newLabels[j]){
          rsp_outin[j] = i;
          cnt++;
        }
    if((cnt == newLabels.size()) == 0){
      std::ostringstream err;
      err<<"The input new labels do not 1-1 correspond to the labels of the tensor.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    bool inorder = true;
    for(int i = 1; i < bondNum; i++)
      if(rsp_outin[i] != i){
        inorder = false;
        break;
      }
    if(inorder && RBondNum == rowBondNum)	//do nothing
      return *this;
    else{
      std::vector<Bond> outBonds;
      bool withoutSymmetry = true;
      for(size_t b = 0; b < bonds.size(); b++){
        outBonds.push_back(bonds[rsp_outin[b]]);
        if(bonds[b].Qnums.size() != 1)
          withoutSymmetry = false;
      }
      for(size_t b = 0; b < bonds.size(); b++){
        if(b < rowBondNum)
          outBonds[b].change(BD_IN);
        else
          outBonds[b].change(BD_OUT);
      }
      UniTensor UniTout(CTYPE, outBonds, name);
      if(status & HAVEELEM){
        if(withoutSymmetry){
          if(!inorder){
            if(ongpu && UniTout.ongpu){
              size_t* perInfo = (size_t*)malloc(bondNum * 2 * sizeof(size_t));
              std::vector<size_t> newAcc(bondNum);
              newAcc[bondNum - 1] = 1;
              perInfo[bondNum - 1] = 1;
              for(int b = bondNum - 1; b > 0; b--){
                newAcc[b - 1] = newAcc[b] * UniTout.bonds[b].Qdegs[0];
                perInfo[b - 1] = perInfo[b] * bonds[b].Qdegs[0];
              }
              for(int b = 0; b < bondNum; b++)
                perInfo[bondNum + rsp_outin[b]] = newAcc[b];
              Complex* des_elem = UniTout.c_elem;
              Complex* src_elem = c_elem;
              reshapeElem(src_elem, bondNum, m_elemNum, perInfo, des_elem);
              free(perInfo);
            }
            else{
              Complex* des_elem = UniTout.c_elem;
              Complex* src_elem = c_elem;
              size_t memsize = m_elemNum * sizeof(Complex);
              if(ongpu){
                src_elem = (Complex*)elemAllocForce(memsize, false);
                elemCopy(src_elem, c_elem, memsize, false, ongpu);
              }
              if(UniTout.ongpu)
                des_elem = (Complex*)elemAllocForce(memsize, false);

              std::vector<size_t> transAcc(bondNum);
              std::vector<size_t> newAcc(bondNum);
              transAcc[bondNum - 1] = 1;
              newAcc[bondNum - 1] = 1;
              for(int b = bondNum - 1; b > 0; b--)
                newAcc[b - 1] = newAcc[b] * UniTout.bonds[b].Qdegs[0];
              std::vector<int> bondDims(bondNum);
              std::vector<int> idxs(bondNum);
              for(int b = 0; b < bondNum; b++){
                transAcc[rsp_outin[b]] = newAcc[b];
                bondDims[b] = bonds[b].Qdegs[0];
                idxs[b] = 0;
              }
              size_t cnt_ot = 0;
              for(size_t i = 0; i < m_elemNum; i++){
                des_elem[cnt_ot] = src_elem[i];
                for(int bend = bondNum - 1; bend >= 0; bend--){
                  idxs[bend]++;
                  if(idxs[bend] < bondDims[bend]){
                    cnt_ot += transAcc[bend];
                    break;
                  }
                  else{
                    cnt_ot -= transAcc[bend] * (idxs[bend] - 1);
                    idxs[bend] = 0;
                  }
                }
              }
              if(ongpu)
                elemFree(src_elem, memsize, false);
              if(UniTout.ongpu){
                elemCopy(UniTout.c_elem, des_elem, memsize, UniTout.ongpu, false);
                elemFree(des_elem, memsize, false);
              }
            }
          }
          else{  //non-symmetry inorder
            size_t memsize = m_elemNum * sizeof(Complex);
            elemCopy(UniTout.c_elem, c_elem, memsize, UniTout.ongpu, ongpu);
          }
        }
        else{
          Real sign = 1.0;
          //For Fermionic system
          std::vector<_Swap> swaps;
          if(Qnum::isFermionic()){
            std::vector<int> inLabelF(bondNum);
            std::vector<int> outLabelF(bondNum);
            std::vector<int> ordF(bondNum);

            for(int b = 0; b < RBondNum; b++){
              inLabelF[b] = labels[b];
              ordF[b] = b;
            }
            for(int b = 0; b < UniTout.RBondNum; b++)
              outLabelF[b] = newLabels[b];
            for(int b = bondNum - 1; b >= RBondNum; b--){
              ordF[b] = bondNum - b + RBondNum - 1;
              inLabelF[ordF[b]] = labels[b];
            }
            for(int b = bondNum - 1; b >= UniTout.RBondNum; b--)
              outLabelF[bondNum - b + UniTout.RBondNum - 1] = newLabels[b];

            std::vector<int> rspF_outin(bondNum);
            for(int i = 0; i < bondNum; i++)
              for(int j = 0; j < bondNum; j++)
                if(inLabelF[i] == outLabelF[j])
                  rspF_outin[j] = i;
            swaps = recSwap(rspF_outin, ordF);
          }
          //End Fermionic system
          std::vector<int> Qin_idxs(bondNum, 0);
          std::vector<int> Qot_idxs(bondNum, 0);
          int Qin_off, Qot_off;
          int tmp;
          int Qin_RQoff, Qin_CQoff;
          int Qot_CQoff, Qot_RQoff;
          size_t sBin_r, sBin_c;	//sub-block of a Qidx
          size_t sBin_rDim, sBin_cDim;	//sub-block of a Qidx
          size_t sBot_cDim;	//sub-block of a Qidx
          size_t sBot_r, sBot_c;
          size_t Bin_cDim, Bot_cDim;
          Complex* Ein_ptr;
          Complex* Eot_ptr;
          std::vector<int> sBin_idxs(bondNum, 0);
          std::vector<int> sBin_sBdims(bondNum, 0);
          std::vector<int> Qot_acc(bondNum, 1);
          std::vector<int> sBot_acc(bondNum, 1);
          for(int b = bondNum	- 1; b > 0; b--)
            Qot_acc[b - 1] = Qot_acc[b] * UniTout.bonds[b].Qnums.size();

          for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
            Qin_off = it->first;
            tmp = Qin_off;
            int qdim;
            for(int b = bondNum - 1; b >= 0; b--){
              qdim = bonds[b].Qnums.size();
              Qin_idxs[b] = tmp % qdim;
              sBin_sBdims[b] = bonds[b].Qdegs[Qin_idxs[b]];
              tmp /= qdim;
            }
            Qot_off = 0;
            for(int b = 0; b < bondNum; b++){
              Qot_idxs[b] = Qin_idxs[rsp_outin[b]];
              Qot_off += Qot_idxs[b] * Qot_acc[b];
            }
            for(int b = bondNum - 1; b > 0; b--)
              sBot_acc[rsp_outin[b-1]] = sBot_acc[rsp_outin[b]] * bonds[rsp_outin[b]].Qdegs[Qot_idxs[b]];
            Qin_RQoff = Qin_off / CQdim;
            Qin_CQoff = Qin_off % CQdim;
            Qot_RQoff = Qot_off / UniTout.CQdim;
            Qot_CQoff = Qot_off % UniTout.CQdim;
            Bin_cDim = RQidx2Blk[Qin_RQoff]->Cnum;
            Bot_cDim = UniTout.RQidx2Blk[Qot_RQoff]->Cnum;
            Ein_ptr = RQidx2Blk[Qin_RQoff]->cm_elem + (RQidx2Off[Qin_RQoff] * Bin_cDim) + CQidx2Off[Qin_CQoff];
            Eot_ptr = UniTout.RQidx2Blk[Qot_RQoff]->cm_elem + (UniTout.RQidx2Off[Qot_RQoff] * Bot_cDim) + UniTout.CQidx2Off[Qot_CQoff];
            sBin_rDim = RQidx2Dim[Qin_RQoff];
            sBin_cDim = CQidx2Dim[Qin_CQoff];
            sBot_cDim = UniTout.CQidx2Dim[Qot_CQoff];
            int cnt_ot = 0;
            sBin_idxs.assign(bondNum, 0);
            if(Qnum::isFermionic()){
              int sign01 = 0;
              for(size_t i = 0; i < swaps.size(); i++)
                sign01 ^= (bonds[swaps[i].b1].Qnums[Qin_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Qin_idxs[swaps[i].b2]].prtF());
              sign = sign01 ? -1.0 : 1.0;
            }
            for(sBin_r = 0; sBin_r < sBin_rDim; sBin_r++)
              for(sBin_c = 0; sBin_c < sBin_cDim; sBin_c++){
                sBot_r = cnt_ot / sBot_cDim;
                sBot_c = cnt_ot % sBot_cDim;
                Eot_ptr[(sBot_r * Bot_cDim) + sBot_c] = sign * Ein_ptr[(sBin_r * Bin_cDim) + sBin_c];
                for(int bend = bondNum - 1; bend >= 0; bend--){
                  sBin_idxs[bend]++;
                  if(sBin_idxs[bend] < sBin_sBdims[bend]){
                    cnt_ot += sBot_acc[bend];
                    break;
                  }
                  else{
                    cnt_ot -= sBot_acc[bend] * (sBin_idxs[bend] - 1);
                    sBin_idxs[bend] = 0;
                  }
                }
              }
          }
        }
        UniTout.status |= HAVEELEM;
      }
      *this = UniTout;
      this->setLabel(newLabels);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(uni10::cflag, std::vector<int>&, int):");
  }
  return *this;
}

Complex UniTensor::at(cflag tp, size_t idx)const{
  try{
    throwTypeError(tp);
    if(!(idx < m_elemNum)){
      std::ostringstream err;
      err<<"Index exceeds the number of elements("<<m_elemNum<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }else if(typeID() == 1){
      std::ostringstream err;
      err<<"This Tensor is REAL. Please use UniTensor::at(rflag, size_t) or UniTensor::at(size_t) instead ";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return getElemAt(idx, c_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(uni10::cflag, size_t):");
    return 0;
  }
}

UniTensor& UniTensor::combineBond(cflag tp, const std::vector<int>&cmbLabels){
  try{
    throwTypeError(tp);
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"There is no bond in the tensor to be combined.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(cmbLabels.size() > 1)){
      return *this;
    }
    std::vector<int> rsp_labels(labels.size(), 0);
    std::vector<int> reduced_labels(labels.size() - cmbLabels.size() + 1, 0);

    std::vector<int> marked(labels.size(), 0);
    std::vector<int> picked(cmbLabels.size(), 0);
    for(size_t p = 0; p < cmbLabels.size(); p++){
      for(size_t l = 0; l < labels.size(); l++){
        if(cmbLabels[p] == labels[l]){
          picked[p] = l;
          marked[l] = 1;
          break;
        }
      }
    }
    int mark = 0;
    for(size_t m = 0; m < marked.size(); m++)
      if(marked[m])
        mark++;
    if(!(mark == cmbLabels.size())){
      std::ostringstream err;
      err<<"The input labels do not match for the labels of the tensor.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int enc = 0;
    int enc_r = 0;
    std::vector<Bond> newBonds;
    int RBnum = 0;
    for(size_t l = 0; l < labels.size(); l++){
      if(marked[l] && l == picked[0]){
        for(size_t ll = 0; ll < cmbLabels.size(); ll++){
          rsp_labels[enc] = cmbLabels[ll];
          enc++;
        }
        std::vector<Bond> tmpBonds;
        for(size_t p = 0; p < picked.size(); p++)
          tmpBonds.push_back(bonds[picked[p]]);
        if(bonds[picked[0]].type() == BD_IN)
          RBnum += picked.size();
        newBonds.push_back(combine(tmpBonds));
        reduced_labels[enc_r] = labels[l];
        enc_r++;
      }
      else if(marked[l] == 0){
        rsp_labels[enc] = labels[l];
        reduced_labels[enc_r] = labels[l];
        if(bonds[l].type() == BD_IN)
          RBnum++;
        newBonds.push_back(bonds[l]);
        enc_r++;
        enc++;
      }
    }
    this->permute(CTYPE, rsp_labels, RBnum);
    UniTensor Tout(CTYPE, newBonds, reduced_labels);

    if(status & HAVEELEM)
      Tout.setElem(c_elem, ongpu);

    *this = Tout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::combineBond(uni10::cflag, std::vector<int>&):");
  }
  return *this;
}

void UniTensor::addGate(cflag tp, const std::vector<_Swap>& swaps){
  try{
    throwTypeError(tp);
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"Adding swap gates to a tensor without bonds(scalar).";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((status & HAVEELEM) == 0){
      std::ostringstream err;
      err<<"Cannot add swap gates to a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int sign = 1;
    int bondNum = bonds.size();
    std::vector<int> Q_idxs(bondNum, 0);
    std::vector<int> Q_Bdims(bondNum, 0);
    for(int b = 0; b < bondNum; b++)
      Q_Bdims[b] = bonds[b].Qnums.size();
    int Q_off;
    int tmp;
    int RQoff, CQoff;
    size_t sB_r, sB_c;	//sub-block of a Qidx
    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
    size_t B_cDim;
    Complex* Eptr;
    for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
      Q_off = it->first;
      tmp = Q_off;
      for(int b = bondNum - 1; b >= 0; b--){
        Q_idxs[b] = tmp % Q_Bdims[b];
        tmp /= Q_Bdims[b];
      }
      RQoff = Q_off / CQdim;
      CQoff = Q_off % CQdim;
      B_cDim = RQidx2Blk[RQoff]->Cnum;
      Eptr = RQidx2Blk[RQoff]->cm_elem + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
      sB_rDim = RQidx2Dim[RQoff];
      sB_cDim = CQidx2Dim[CQoff];

      int sign01 = 0;
      for(size_t i = 0; i < swaps.size(); i++)
        sign01 ^= (bonds[swaps[i].b1].Qnums[Q_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Q_idxs[swaps[i].b2]].prtF());
      sign = sign01 ? -1 : 1;

      for(sB_r = 0; sB_r < sB_rDim; sB_r++)
        for(sB_c = 0; sB_c < sB_cDim; sB_c++)
          Eptr[(sB_r * B_cDim) + sB_c] *= sign;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::addGate(uni10::cflag, std::vector<_Swap>&):");
  }
}

Complex UniTensor::trace(cflag tp)const{
  try{
    throwTypeError(tp);
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot trace a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(status & HAVEBOND){
      Complex trVal(0, 0);
      for(std::map<Qnum, Block>::const_iterator it = blocks.begin() ; it != blocks.end(); it++ ){
        if(!(it->second.Rnum == it->second.Cnum)){
          std::ostringstream err;
          err<<"Cannot trace a non-square block.";
          throw std::runtime_error(exception_msg(err.str()));
        }
        trVal += it->second.trace(CTYPE);
      }
      return trVal;
    }
    else{
      return getElemAt(0, c_elem, ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::trace(uni10::cflag ):");
    return 0;
  }
}

UniTensor& UniTensor::partialTrace(cflag tp, int la, int lb){
  try{
    throwTypeError(tp);
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot trace bonds of a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(bonds.size() > 2)){
      std::ostringstream err;
      err<<"The number of bonds must larger than 2 for performing partialTrace.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int bondNum = bonds.size();
    std::vector<Bond> newBonds;
    std::vector<int>newLabels(bondNum - 2, 0);
    std::vector<int>rsp_labels(bondNum);
    int ia, ib;
    int enc = 0;
    for(size_t l = 0; l < labels.size(); l++){
      if(labels[l] == la)
        ia = l;
      else if(labels[l] == lb)
        ib = l;
      else{
        newBonds.push_back(bonds[l]);
        newLabels[enc] = labels[l];
        rsp_labels[enc] = labels[l];
        enc++;
      }
    }
    if(!(enc == newLabels.size())){
      std::ostringstream err;
      err<<"Cannot find the two bonds with the given two labels.";
      throw std::runtime_error(exception_msg(err.str()));
    }

    UniTensor Tt(CTYPE, newBonds, newLabels);
    rsp_labels[bondNum - 2] = labels[ia];
    rsp_labels[bondNum - 1] = labels[ib];
    ia = bondNum - 2;
    ib = bondNum - 1;
    this->permute(CTYPE, rsp_labels, Tt.RBondNum);
    std::vector<int> Q_acc(bondNum, 1);
    for(int b = bondNum - 1; b > 0; b--)
      Q_acc[b - 1] = Q_acc[b] * bonds[b].Qnums.size();
    int tQdim = bonds[ia].Qnums.size();
    /*Sanity Check*/
    if(tQdim == bonds[ib].Qnums.size()){
      std::ostringstream err;
      err<<"The bonds of the given two labels does not match for trace.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    Qnum q0(0, PRT_EVEN);
    for(int q = 0; q < tQdim; q++){
      if(!((bonds[ia].Qnums[q] * bonds[ib].Qnums[q] == q0) && (bonds[ia].Qdegs[q] == bonds[ib].Qdegs[q]))){
        std::ostringstream err;
        err<<"The bonds of the given two labels does not match for trace.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    }
    /*END*/
    int tBnum = Tt.bonds.size();
    std::vector<int> Qt_Bdims(tBnum, 0);
    for(int b = 0; b < tBnum; b++)
      Qt_Bdims[b] = Tt.bonds[b].Qnums.size();

    int Qt_off;
    int Q_off;
    int Qt_RQoff, Qt_CQoff;
    int Q_RQoff, Q_CQoff;
    size_t sBt_rDim, sBt_cDim;	//sub-block of a Qidx of Tt
    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
    size_t Bt_cDim;
    Complex* Et_ptr;
    std::vector<Complex*> E_offs(tQdim);
    std::vector<size_t> B_cDims(tQdim);
    int tQdim2 = tQdim * tQdim;
    int Qenc = Q_acc[ia] + Q_acc[ib];
    for(std::map<int, size_t>::iterator it = Tt.QidxEnc.begin(); it != Tt.QidxEnc.end(); it++){
      Qt_off = it->first;
      Qt_RQoff = Qt_off / Tt.CQdim;
      Qt_CQoff = Qt_off % Tt.CQdim;
      Bt_cDim = Tt.RQidx2Blk[Qt_RQoff]->Cnum;
      Et_ptr = Tt.RQidx2Blk[Qt_RQoff]->cm_elem + (Tt.RQidx2Off[Qt_RQoff] * Bt_cDim) + Tt.CQidx2Off[Qt_CQoff];
      sBt_rDim = Tt.RQidx2Dim[Qt_RQoff];
      sBt_cDim = Tt.CQidx2Dim[Qt_CQoff];

      for(int q = 0; q < tQdim; q++){
        Q_off = Qt_off * tQdim2 + q * Qenc;
        Q_RQoff = Q_off / CQdim;
        Q_CQoff = Q_off % CQdim;
        B_cDims[q] = RQidx2Blk[Q_RQoff]->Cnum;
        E_offs[q] = RQidx2Blk[Q_RQoff]->cm_elem + (RQidx2Off[Q_RQoff] * B_cDims[q]) + CQidx2Off[Q_CQoff];
      }
      int tQdeg, sB_c_off;
      Complex trVal;
      for(size_t sB_r = 0; sB_r < sBt_rDim; sB_r++)
        for(size_t sB_c = 0; sB_c < sBt_cDim; sB_c++){
          trVal = 0;
          for(int q = 0; q < tQdim; q++){
            tQdeg = bonds[ia].Qdegs[q];
            sB_c_off = sB_c * (tQdeg * tQdeg);
            for(int t = 0; t < tQdeg; t++){
              trVal += E_offs[q][(sB_r * B_cDims[q]) + sB_c_off + t * (tQdeg + 1)];
            }
          }
          Et_ptr[sB_r * Bt_cDim + sB_c] = trVal;
        }
      Tt.status |= HAVEELEM;
    }
    *this = Tt;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::partialTrace(uni10::cflag, int, int):");
  }
  return *this;
}

UniTensor& UniTensor::assign(cflag tp, const std::vector<Bond>& _bond){
  try{
    throwTypeError(tp);
    UniTensor T(CTYPE, _bond);
    *this = T;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::assign(uni10::cflag, std::vector<Bond>&):");
  }
  return *this;
}

Complex UniTensor::at(cflag tp, const std::vector<int>& idxs)const{
  try{
    throwTypeError(tp);
    if(typeID() == 1){
      std::ostringstream err;
      err<<"This Tensor is REAL. Please use UniTensor::at(cflag, size_t) or UniTensor::at(size_t) instead ";
      throw std::runtime_error(exception_msg(err.str()));
    }
    std::vector<size_t> _idxs(idxs.size());
    for(size_t i = 0; i < idxs.size(); i++)
      _idxs[i] = idxs[i];
    return at(CTYPE, _idxs);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::cflag, std::vector<int>&):");
    return 0;
  }
}

Complex UniTensor::at(cflag tp, const std::vector<size_t>& idxs)const{
  try{
    throwTypeError(tp);
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"The tensor is a scalar. Use UniTensor::operator() instead.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(idxs.size() == bonds.size())){
      std::ostringstream err;
      err<<"The size of input indices array does not match with the number of the bonds.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 1){
      std::ostringstream err;
      err<<"This Tensor is REAL. Please use UniTensor::at(rflag, const vector<size_t>&) or UniTensor::at(const vector<size_t>& ) instead ";
      throw std::runtime_error(exception_msg(err.str()));
    }

    int bondNum = bonds.size();
    std::vector<int> Qidxs(bondNum, 0);
    for(int b = 0; b < bondNum; b++){
      if(!(idxs[b] < bonds[b].dim())){
        std::ostringstream err;
        err<<"The input indices are out of range.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      for(int q = bonds[b].offsets.size() - 1; q >= 0; q--){
        if(idxs[b] < bonds[b].offsets[q])
          continue;
        Qidxs[b] = q;
        break;
      }
    }
    std::vector<int> Q_acc(bondNum, 1);
    for(int b = bondNum	- 1; b > 0; b--)
      Q_acc[b - 1] = Q_acc[b] * bonds[b].Qnums.size();
    int Qoff = 0;
    for(int b = 0; b < bondNum; b++)
      Qoff += Q_acc[b] * Qidxs[b];

    if(QidxEnc.find(Qoff) != QidxEnc.end()){
      int Q_RQoff = Qoff / CQdim;
      int Q_CQoff = Qoff % CQdim;
      Block* blk = RQidx2Blk.find(Q_RQoff)->second;
      size_t B_cDim = blk->Cnum;
      size_t sB_cDim = CQidx2Dim.find(Q_CQoff)->second;
      size_t blkRoff = RQidx2Off.find(Q_RQoff)->second;
      size_t blkCoff = CQidx2Off.find(Q_CQoff)->second;
      Complex* boff = blk->cm_elem + (blkRoff * B_cDim) + blkCoff;
      int cnt = 0;
      std::vector<int> D_acc(bondNum, 1);
      for(int b = bondNum	- 1; b > 0; b--)
        D_acc[b - 1] = D_acc[b] * bonds[b].Qdegs[Qidxs[b]];
      for(int b = 0; b < bondNum; b++)
        cnt += (idxs[b] - bonds[b].offsets[Qidxs[b]]) * D_acc[b];
      return boff[(cnt / sB_cDim) * B_cDim + cnt % sB_cDim];
    }
    else{
      return 0.0;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(uni10::cflag, const std::vector<size_t>&):");
    return 0;
  }
}

/*********************  Private **********************/

void UniTensor::initUniT(cflag tp){ //GPU
  r_flag= RNULL;
  c_flag= CTYPE;
  if(bonds.size()){
    m_elemNum = grouping(CTYPE);
    if(!(blocks.size() > 0)){ //No block in Tensor, Error!
      std::ostringstream err;
      err<<"There is no symmetry block with the given bonds:\n";
      for(size_t b = 0; b < bonds.size(); b++)
        err<<"    "<<bonds[b];
      throw std::runtime_error(exception_msg(err.str()));
    }
    labels.assign(bonds.size(), 0);
    for(size_t b = 0; b < bonds.size(); b++)
      labels[b] = b;
    status |= HAVEBOND;
  }
  else{
    Qnum q0(0);
    blocks[q0] = Block(CTYPE, 1, 1);
    RBondNum = 0;
    RQdim = 0;
    CQdim = 0;
    m_elemNum = 1;
    status |= HAVEELEM;
  }
  elem = NULL;
  c_elem = NULL;

  ELEMNUM += m_elemNum;
  COUNTER++;
  if(ELEMNUM > MAXELEMNUM)
    MAXELEMNUM = ELEMNUM;
  if(m_elemNum > MAXELEMTEN)
    MAXELEMTEN = m_elemNum;

  TelemAlloc(CTYPE);
  initBlocks(CTYPE);
  TelemBzero(CTYPE);
}

size_t UniTensor::grouping(cflag tp){
  blocks.clear();
  int row_bondNum = 0;
  int col_bondNum = 0;
  RQdim = 1;
  CQdim = 1;
  bool IN_BONDS_BEFORE_OUT_BONDS = true;
  for(size_t i = 0; i < bonds.size(); i++){
    if(bonds[i].type() == BD_IN){
      if(!(IN_BONDS_BEFORE_OUT_BONDS == true)){
        std::ostringstream err;
        err<<"Error in the input bond array: BD_OUT bonds must be placed after all BD_IN bonds.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      RQdim *= bonds[i].Qnums.size();
      row_bondNum++;
    }
    else{
      CQdim *= bonds[i].Qnums.size();
      col_bondNum++;
      IN_BONDS_BEFORE_OUT_BONDS = false;
    }
  }
  RBondNum = row_bondNum;
  std::map<Qnum,size_t> row_QnumMdim;
  std::vector<int> row_offs(row_bondNum, 0);
  std::map<Qnum,std::vector<int> > row_Qnum2Qidx;
  Qnum qnum;
  size_t dim;
  int boff = 0;
  std::vector<size_t>tmpRQidx2Dim(RQdim, 1);
  std::vector<size_t>tmpCQidx2Dim(CQdim, 1);
  std::vector<size_t>tmpRQidx2Off(RQdim, 0);
  std::vector<size_t>tmpCQidx2Off(CQdim, 0);
  if(row_bondNum){
    while(1){
      qnum.assign();
      dim = 1;
      for(int b = 0; b < row_bondNum; b++){
        qnum = qnum * bonds[b].Qnums[row_offs[b]];
        dim *= bonds[b].Qdegs[row_offs[b]];
      }
      if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
        tmpRQidx2Off[boff] = row_QnumMdim[qnum];
        tmpRQidx2Dim[boff] = dim;
        row_QnumMdim[qnum] += dim;
      }
      else{
        tmpRQidx2Off[boff] = 0;
        tmpRQidx2Dim[boff] = dim;
        row_QnumMdim[qnum] = dim;
      }
      row_Qnum2Qidx[qnum].push_back(boff);
      boff++;
      int bidx;
      for(bidx = row_bondNum - 1; bidx >= 0; bidx--){
        row_offs[bidx]++;
        if(row_offs[bidx] < bonds[bidx].Qnums.size())
          break;
        else
          row_offs[bidx] = 0;
      }
      if(bidx < 0)	//run over all row_bond offsets
        break;
    }
  }
  else{
    qnum.assign();
    row_QnumMdim[qnum] = 1;
    row_Qnum2Qidx[qnum].push_back(0);
  }
  std::map<Qnum,size_t> col_QnumMdim;
  std::vector<int> col_offs(col_bondNum, 0);
  std::map<Qnum,std::vector<int> > col_Qnum2Qidx;
  boff = 0;
  if(col_bondNum){
    while(1){
      qnum.assign();
      dim = 1;
      for(int b = 0; b < col_bondNum; b++){
        qnum = qnum * bonds[b + row_bondNum].Qnums[col_offs[b]];
        dim *= bonds[b + row_bondNum].Qdegs[col_offs[b]];
      }
      if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
        if(col_QnumMdim.find(qnum) != col_QnumMdim.end()){
          tmpCQidx2Off[boff] = col_QnumMdim[qnum];
          tmpCQidx2Dim[boff] = dim;
          col_QnumMdim[qnum] += dim;
        }
        else{
          tmpCQidx2Off[boff] = 0;
          tmpCQidx2Dim[boff] = dim;
          col_QnumMdim[qnum] = dim;
        }
        col_Qnum2Qidx[qnum].push_back(boff);
      }
      boff++;
      int bidx;
      for(bidx = col_bondNum - 1; bidx >= 0; bidx--){
        col_offs[bidx]++;
        if(col_offs[bidx] < bonds[bidx + row_bondNum].Qnums.size())
          break;
        else
          col_offs[bidx] = 0;
      }
      if(bidx < 0)	//run over all row_bond offsets
        break;
    }
  }
  else{
    qnum.assign();
    if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
      col_QnumMdim[qnum] = 1;
      col_Qnum2Qidx[qnum].push_back(0);
    }
  }

  std::map<Qnum,size_t>::iterator it;
  std::map<Qnum,size_t>::iterator it2;
  std::set<int> Qidx;
  int qidx;
  size_t off = 0;
  for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
    it = row_QnumMdim.find(it2->first);
    Block blk(CTYPE, it->second, it2->second); // blk(Rnum, Cnum);
    off += blk.Rnum * blk.Cnum;
    blocks[it->first] = blk;
    Block* blkptr = &(blocks[it->first]);
    std::vector<int>& tmpRQidx = row_Qnum2Qidx[it->first];
    std::vector<int>& tmpCQidx = col_Qnum2Qidx[it->first];
    for(size_t i = 0; i < tmpRQidx.size(); i++){
      RQidx2Blk[tmpRQidx[i]] = blkptr;
      for(size_t j = 0; j < tmpCQidx.size(); j++){
        RQidx2Dim[tmpRQidx[i]] = tmpRQidx2Dim[tmpRQidx[i]];
        RQidx2Off[tmpRQidx[i]] = tmpRQidx2Off[tmpRQidx[i]];
        CQidx2Dim[tmpCQidx[j]] = tmpCQidx2Dim[tmpCQidx[j]];
        CQidx2Off[tmpCQidx[j]] = tmpCQidx2Off[tmpCQidx[j]];
        qidx = tmpRQidx[i] * CQdim + tmpCQidx[j];
        Qidx.insert(qidx);
      }
    }
  }
  size_t elemEnc = 0;
  for(std::map<int, size_t>::iterator itr = RQidx2Dim.begin(); itr != RQidx2Dim.end(); itr++)
    for(std::map<int, size_t>::iterator itc = CQidx2Dim.begin(); itc != CQidx2Dim.end(); itc++){
      qidx = itr->first * CQdim + itc->first;
      if(Qidx.find(qidx) != Qidx.end()){
        QidxEnc[qidx] = elemEnc;
        elemEnc += RQidx2Dim[itr->first] * CQidx2Dim[itc->first];
      }
    }
  return off;
}

void UniTensor::initBlocks(cflag tp){
  size_t offset = 0;
  for(std::map<Qnum, Block>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
    it->second.r_flag = RNULL;
    it->second.c_flag = CTYPE;
    it->second.cm_elem = &(c_elem[offset]);
    it->second.ongpu = ongpu;
    offset += it->second.Rnum * it->second.Cnum;
  }
}

void UniTensor::TelemAlloc(cflag tp){
  c_elem = (Complex*)elemAlloc(sizeof(Complex) * m_elemNum, ongpu);
}

void UniTensor::TelemBzero(cflag tp){
  elemBzero(c_elem, sizeof(Complex) * m_elemNum, ongpu);
}

/************* developping *************/

Real UniTensor::norm(cflag tp) const{
  try{
    throwTypeError(tp);
    return vectorNorm(c_elem, elemNum(), 1, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::norm(uni10::cflag ):");
    return 0;
  }
}

UniTensor& UniTensor::normalize(cflag tp){
  try{
    throwTypeError(tp);
    Real norm = vectorNorm(c_elem, elemNum(), 1, ongpu);
    vectorScal((1./norm), c_elem, elemNum(), ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::normalize(uni10::cflag ):");
  }
  return *this;
}


std::vector<UniTensor> UniTensor::hosvd(cflag tp, int* _group_labels, int* _groups, size_t _groupsSize, std::vector<Matrix>& Ls)const{
  try{
    std::vector<int> group_labels(_group_labels, _group_labels+this->bondNum());
    std::vector<int> groups(_groups, _groups+_groupsSize);
    return hosvd(tp, group_labels, groups, Ls);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(uni10::cflag ,int* ,int* ,size_t ,std::vector<uni10::Matrix>& Ls)const;");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(cflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<Matrix>& Ls)const{
  try{
    bool withoutSymmetry = true;
    for(size_t b = 0; b < bonds.size(); b++){
      if(bonds[b].Qnums.size() != 1)
        withoutSymmetry = false;
    }
    if(!withoutSymmetry){
      std::ostringstream err;
      err<<"The tensor has symmetry quantum numbers. Cannot use non-symmetry version hosvd(size_t, std::vector<uni10::Matrix>&)";
      err<<"\n  Hint: Use UniTensor::hosvd(uni10::cflag, std::vector<int>&, std::vector<int>&, std::vector<uni10::Matrix> >&)";
      throw std::runtime_error(exception_msg(err.str()));
    }
    std::vector<std::map<Qnum, Matrix> > symLs;
    const std::vector<UniTensor>& outs = hosvd(tp , group_labels, groups, symLs, true);
    Ls.clear();
    Qnum q0(0);
    for(size_t i = 0; i < symLs.size(); i++)
      Ls.push_back(symLs[i][q0]);
    return outs;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(uni10::cflag, std::vector<int>, std::vector<int>, std::vector<Matrix>&):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(cflag tp, int* _group_labels, int* _groups, size_t _groupsSize, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const{
  try{
    std::vector<int> group_labels(_group_labels, _group_labels+this->bondNum());
    std::vector<int> groups(_groups, _groups+_groupsSize);
    return hosvd(tp, group_labels, groups, Ls, returnL);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(uni10::cflag ,int* ,int* ,size_t ,std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const;");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(cflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const{
  throwTypeError(tp);
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"Cannot perform higher order SVD on a tensor without bonds(scalar).";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((status & HAVEELEM) == 0){
      std::ostringstream err;
      err<<"Cannot perform higher order SVD on a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(group_labels.size() != this->bondNum()){
      std::ostringstream err;
      err<<"The size of Group labels is not match to the number of tensor's bond.";
      throw std::runtime_error(exception_msg(err.str()));
    }

    UniTensor T(*this);
    size_t groupElemNum=0;
    for(size_t n = 0; n < groups.size(); n++)
      groupElemNum+=groups[n];

    if(returnL)
      Ls.assign(groups.size(), std::map<Qnum, Matrix>());
    std::vector<UniTensor> Us;
    UniTensor S(T);

    std::vector<int> lrsp_labels = group_labels;
    std::vector<int> rsp_labels = group_labels;

    int min = *std::min_element(rsp_labels.begin(), rsp_labels.end());

    for(size_t m = 0; m < groups.size(); m++){
      int pos=0;
      for(size_t l = 0; l < groupElemNum; l++){
        if(l >= groupElemNum-groups[m])
          rsp_labels[pos] = lrsp_labels[l-(groupElemNum-groups[m])];
        else
          rsp_labels[pos] = lrsp_labels[l+groups[m]];
        pos++;
      }
      T.permute(CTYPE, lrsp_labels, groups[m]);
      std::vector<Bond> bonds(T.bonds.begin(), T.bonds.begin() + groups[m]);
      bonds.push_back(combine(bonds).dummy_change(BD_OUT));
      Us.push_back(UniTensor(CTYPE, bonds));
      for(std::map<Qnum, Block>::iterator it = T.blocks.begin(); it != T.blocks.end(); it++){
        std::vector<Matrix> svd = it->second.svd(CTYPE);
        Us[m].putBlock(it->first, svd[0]);
        if(returnL)
          Ls[m][it->first] = svd[1];
      }
      for(int c = 0; c < groups[m]; c++)
        Us[m].labels[c] = lrsp_labels[c];
      Us[m].labels[groups[m]] = min -m - 1;
      UniTensor UT = Us[m];
      S *= UT.cTranspose(CTYPE);
      lrsp_labels = rsp_labels;
    } 
    Us.push_back(S);
    return Us;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(uni10::cflag , std::vector<int>& , std::vector<int>& , std::vector<std::map<Qnum, Matrix> >& , bool returnL)const:");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(cflag tp, size_t modeNum, size_t fixedNum)const{
  try{
    std::vector<std::map<Qnum, Matrix> > symLs;
    std::vector<int> group_labels=this->labels;
    std::vector<int> groups(modeNum, (this->bondNum()-fixedNum)/modeNum);
    return hosvd(tp, group_labels, groups, symLs, false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(uni10::cflag, size_t, size_t = 0):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(cflag tp, size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const{
  try{
    std::vector<std::map<Qnum, Matrix> > symLs;
    std::vector<int> group_labels=this->labels;
    std::vector<int> groups(modeNum, (this->bondNum()-fixedNum)/modeNum );
    return hosvd(tp, group_labels, groups, Ls, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(uni10::cflag, size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(cflag tp, size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const{
  try{
    bool withoutSymmetry = true;
    for(size_t b = 0; b < bonds.size(); b++){
      if(bonds[b].Qnums.size() != 1)
        withoutSymmetry = false;
    }
    if(!withoutSymmetry){
      std::ostringstream err;
      err<<"The tensor has symmetry quantum numbers. Cannot use non-symmetry version hosvd(size_t, std::vector<uni10::Matrix>&)";
      err<<"\n  Hint: Use UniTensor::hosvd(uni10::cflag, size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&)";
      throw std::runtime_error(exception_msg(err.str()));
    }
    std::vector<std::map<Qnum, Matrix> > symLs;
    const std::vector<UniTensor>& outs = hosvd(tp , modeNum, fixedNum, symLs);
    Ls.clear();
    Qnum q0(0);
    for(size_t i = 0; i < symLs.size(); i++)
      Ls.push_back(symLs[i][q0]);
    return outs;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(uni10::cflag, size_t, size_t, std::vector<Matrix>&):");
    return std::vector<UniTensor>();
  }
}

}; /* namespace uni10 */

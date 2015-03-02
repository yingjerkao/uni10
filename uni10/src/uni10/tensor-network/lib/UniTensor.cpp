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
#ifndef UNI10_PURE_REAL
#include <uni10/tensor-network/CMatrix.h>
#include <uni10/tensor-network/CUniTensor.h>
#endif
#include <uni10/tensor-network/UniTensor.h>

#ifndef UNI10_DTYPE
#define UNI10_DTYPE double
#endif
#ifndef UNI10_TENSOR
#define UNI10_TENSOR UniTensor
#endif
#ifndef UNI10_MATRIX
#define UNI10_MATRIX Matrix
#endif
#ifndef UNI10_BLOCK
#define UNI10_BLOCK Block
#endif

namespace uni10{

int64_t UNI10_TENSOR::ELEMNUM = 0;
int UNI10_TENSOR::COUNTER = 0;
size_t UNI10_TENSOR::MAXELEMNUM = 0;
size_t UNI10_TENSOR::MAXELEMTEN = 0;

UNI10_TENSOR::UNI10_TENSOR(): status(0){
  try{
    initUniT();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor():");
  }
}

UNI10_TENSOR::UNI10_TENSOR(UNI10_DTYPE val): status(0){ //GPU
  try{
    initUniT();
    setElemAt(0, val, elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(double val):");
  }
}

UNI10_TENSOR::UNI10_TENSOR(const std::string& fname): status(0){ //GPU
  try{
    int namemax = 32;
    if(fname.size() > namemax)
      name = fname.substr(fname.size() - namemax);
    else
      name = fname;
    FILE* fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int st;
    fread(&st, 1, sizeof(int), fp);
    int bondNum;
    fread(&bondNum, 1, sizeof(bondNum), fp);  //OUT: bondNum(4 bytes)
    size_t qnum_sz;
    fread(&qnum_sz, 1, sizeof(size_t), fp);	//OUT: sizeof(Qnum)
    if(!(qnum_sz == sizeof(Qnum))){
      std::ostringstream err;
      err<<"Error in reading file '"<<fname<<"' in.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    for(int b = 0; b < bondNum; b++){
      int num_q;
      bondType tp;
      fread(&tp, 1, sizeof(bondType), fp);	//OUT: Number of Qnums in the bond(4 bytes)
      fread(&num_q, 1, sizeof(int), fp);	//OUT: Number of Qnums in the bond(4 bytes)
      Qnum q0;
      std::vector<Qnum> qnums(num_q, q0);
      fread(&(qnums[0]), num_q, qnum_sz, fp);
      std::vector<int> qdegs(num_q, 0);
      fread(&(qdegs[0]), num_q, sizeof(int), fp);
      std::vector<Qnum> tot_qnums;
      for(int q = 0; q < num_q; q++)
        for(int d = 0; d < qdegs[q]; d++)
          tot_qnums.push_back(qnums[q]);
      Bond bd(tp, tot_qnums);
      bonds.push_back(bd);
    }
    initUniT();
    int num_l;
    fread(&num_l, 1, sizeof(int), fp);	//OUT: Number of Labels in the Tensor(4 bytes)
    if(!(num_l == bonds.size())){
      std::ostringstream err;
      err<<"Error in reading file '"<<fname<<"' in.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    labels.assign(num_l, 0);
    fread(&(labels[0]), num_l, sizeof(int), fp);
    if(st & HAVEELEM){
      UNI10_DTYPE *tmp_elem = elem;
      size_t memsize = m_elemNum * sizeof(UNI10_DTYPE);
      if(ongpu)
        tmp_elem = (UNI10_DTYPE*)malloc(memsize);
      size_t num_el;
      fread(&num_el, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
      fread(tmp_elem, m_elemNum, sizeof(UNI10_DTYPE), fp);
      if(ongpu){
        elemCopy(elem, tmp_elem, memsize, ongpu, false);
        free(tmp_elem);
      }
      status |= HAVEELEM;
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::string&):");
  }
}

UNI10_TENSOR::UNI10_TENSOR(const std::vector<Bond>& _bonds, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::string& = \"\"):");
  }
}

UNI10_TENSOR::UNI10_TENSOR(const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT();
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::vector<int>&, std::string& = \"\"):");
  }
}
UNI10_TENSOR::UNI10_TENSOR(const std::vector<Bond>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT();
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, int*, std::string& = \"\"):");
  }
}

UNI10_TENSOR::UNI10_TENSOR(const UNI10_TENSOR& UniT): //GPU
  status(UniT.status), bonds(UniT.bonds), blocks(UniT.blocks), labels(UniT.labels), name(UniT.name),
  RBondNum(UniT.RBondNum), RQdim(UniT.RQdim), CQdim(UniT.CQdim), m_elemNum(UniT.m_elemNum), elem(NULL),
  QidxEnc(UniT.QidxEnc), RQidx2Off(UniT.RQidx2Off), CQidx2Off(UniT.CQidx2Off), RQidx2Dim(UniT.RQidx2Dim), RQidx2Blk(UniT.RQidx2Blk), CQidx2Dim(UniT.CQidx2Dim){
    try{
      elem = (UNI10_DTYPE*)elemAlloc(sizeof(UNI10_DTYPE) * m_elemNum, ongpu);
      std::map<Qnum, UNI10_BLOCK>::const_iterator it2;
      std::map<const UNI10_BLOCK* , UNI10_BLOCK*> blkmap;
      for (std::map<Qnum, UNI10_BLOCK>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
        it->second.m_elem = &(elem[(it->second.m_elem - UniT.elem)]);

        it2 = UniT.blocks.find(it->first);
        blkmap[(&(it2->second))] = &(it->second);
      }
      if(UniT.status & HAVEBOND){
        for(std::map<int, UNI10_BLOCK*>::iterator it = RQidx2Blk.begin(); it != RQidx2Blk.end(); it++)
          it->second = blkmap[it->second];
      }
      ELEMNUM += m_elemNum;
      COUNTER++;
      if(ELEMNUM > MAXELEMNUM)
        MAXELEMNUM = ELEMNUM;
      if(m_elemNum > MAXELEMTEN)
        MAXELEMTEN = m_elemNum;
      elemCopy(elem, UniT.elem, sizeof(UNI10_DTYPE) * UniT.m_elemNum, ongpu, UniT.ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In copy constructor UniTensor::UniTensor(uni10::UniTensor&):");
    }
  }

UNI10_TENSOR::UNI10_TENSOR(const UNI10_BLOCK& blk): status(0){
  try{
    Bond bdi(BD_IN, blk.Rnum);
    Bond bdo(BD_OUT, blk.Cnum);
    bonds.push_back(bdi);
    bonds.push_back(bdo);
    initUniT();
    this->putBlock(blk);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(uni10::Block&");
  }
}

UNI10_TENSOR::~UNI10_TENSOR(){
	elemFree(elem, sizeof(UNI10_DTYPE) * m_elemNum, ongpu);
	ELEMNUM -= m_elemNum;
	COUNTER--;
}

UNI10_TENSOR& UNI10_TENSOR::operator=(const UNI10_TENSOR& UniT){ //GPU
  try{
    bonds = UniT.bonds;
    blocks = UniT.blocks;
    labels = UniT.labels;
    name = UniT.name;
    RBondNum = UniT.RBondNum;
    RQdim = UniT.RQdim;
    CQdim = UniT.CQdim;
    QidxEnc = UniT.QidxEnc;
    RQidx2Off = UniT.RQidx2Off;
    CQidx2Off = UniT.CQidx2Off;
    RQidx2Dim = UniT.RQidx2Dim;
    CQidx2Dim = UniT.CQidx2Dim;
    RQidx2Blk = UniT.RQidx2Blk;

    ELEMNUM -= m_elemNum;	//free original memory
    if(elem != NULL)
      elemFree(elem, sizeof(UNI10_DTYPE) * m_elemNum, ongpu);
    status = UniT.status;
    m_elemNum = UniT.m_elemNum;
    elem = (UNI10_DTYPE*)elemAlloc(sizeof(UNI10_DTYPE) * m_elemNum, ongpu);
    std::map<Qnum, UNI10_BLOCK>::const_iterator it2;
    std::map< const UNI10_BLOCK* , UNI10_BLOCK*> blkmap;
    for (std::map<Qnum, UNI10_BLOCK>::iterator it = blocks.begin(); it != blocks.end(); it++ ){ // blocks here is UniT.blocks
      it->second.m_elem = &(elem[it->second.m_elem - UniT.elem]);

      it2 = UniT.blocks.find(it->first);
      blkmap[&(it2->second)] = &(it->second);
    }

    if(UniT.status & HAVEBOND){
      for(std::map<int, UNI10_BLOCK*>::iterator it = RQidx2Blk.begin(); it != RQidx2Blk.end(); it++)
        it->second = blkmap[it->second];
    }

    ELEMNUM += m_elemNum;
    if(ELEMNUM > MAXELEMNUM)
      MAXELEMNUM = ELEMNUM;
    if(m_elemNum > MAXELEMTEN)
      MAXELEMTEN = m_elemNum;
    elemCopy(elem, UniT.elem, sizeof(UNI10_DTYPE) * UniT.m_elemNum, ongpu, UniT.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator=(uni10::UniTensor&):");
  }
  return *this;
}

UNI10_TENSOR& UNI10_TENSOR::assign(const std::vector<Bond>& _bond){
  try{
    UNI10_TENSOR T(_bond);
	  *this = T;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::assign(std::vector<Bond>&):");
  }
  return *this;
}

void UNI10_TENSOR::setLabel(const std::vector<int>& newLabels){
  try{
    std::set<int> labelS(newLabels.begin(), newLabels.end());
    if(!(bonds.size() == labelS.size())){
      throw std::runtime_error(exception_msg("The size of input vector(labels) does not match for the number of bonds."));
    }
    labels = newLabels;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setLabel(std::vector<int>&):");
  }
}

void UNI10_TENSOR::setLabel(int* newLabels){
  try{
    std::vector<int> labels(newLabels, newLabels + bonds.size());
    setLabel(labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setLabel(int*):");
  }
}

std::vector<int> UNI10_TENSOR::label()const{
	return labels;
}

int UNI10_TENSOR::label(size_t idx)const{
  try{
    if(!(idx < labels.size())){
      std::ostringstream err;
      err<<"Index exceeds the number of the bonds("<<bonds.size()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::label(size_t):");
  }
	return labels[idx];
}

std::string UNI10_TENSOR::getName(){
	return name;
}

void UNI10_TENSOR::setName(const std::string& _name){
	name = _name;
}


size_t UNI10_TENSOR::bondNum()const{return bonds.size();}

size_t UNI10_TENSOR::inBondNum()const{return RBondNum;}

std::vector<Bond> UNI10_TENSOR::bond()const{
	return bonds;
}

Bond UNI10_TENSOR::bond(size_t idx)const{
  try{
    if(!(idx < bonds.size())){
      std::ostringstream err;
      err<<"Index exceeds the number of the bonds("<<bonds.size()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::bond(size_t):");
  }
  return bonds[idx];
}

size_t UNI10_TENSOR::elemNum()const{return m_elemNum;}

UNI10_MATRIX UNI10_TENSOR::getRawElem()const{
  try{
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
      std::vector<UNI10_DTYPE > rawElem;
      while(1){
        rawElem.push_back(at(idxs));
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
      return UNI10_MATRIX(rowNum, colNum, &rawElem[0]);
    }
    else if(status & HAVEELEM)
      return UNI10_MATRIX(1, 1, elem);
    else
      return UNI10_MATRIX();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getRawElem():");
    return UNI10_MATRIX();
  }
}

void UNI10_TENSOR::setRawElem(const UNI10_BLOCK& blk){
  try{
    setRawElem(blk.getElem());
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(uni10::Block&):");
  }
}

void UNI10_TENSOR::setRawElem(const std::vector<UNI10_DTYPE >& rawElem){
  try{
    setRawElem(&rawElem[0]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(std::vector<double>&):");
  }
}
void UNI10_TENSOR::setRawElem(const UNI10_DTYPE* rawElem){
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"Setting elements to a tensor without bonds is not supported.";
      throw std::runtime_error(exception_msg(err.str()));
    }
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
    UNI10_DTYPE* work = elem;
    if(ongpu){
      work = (UNI10_DTYPE*)malloc(m_elemNum * sizeof(UNI10_DTYPE));
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
      E_off = (RQidx2Blk[RQoff]->m_elem - elem) + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
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
      elemCopy(elem, work, m_elemNum * sizeof(UNI10_DTYPE), ongpu, false);
      free(work);
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(double*):");
  }
}

UNI10_DTYPE UNI10_TENSOR::at(const std::vector<int>& idxs)const{
  try{
    std::vector<size_t> _idxs(idxs.size());
    for(int i = 0; i < idxs.size(); i++)
      _idxs[i] = idxs[i];
    return at(_idxs);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::vector<int>&):");
    return 0;
  }
}
UNI10_DTYPE UNI10_TENSOR::at(const std::vector<size_t>& idxs)const{
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"The tensor is a scalar. Use UniTensor::operator[] instead.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(idxs.size() == bonds.size())){
      std::ostringstream err;
      err<<"The size of input indices array does not match with the number of the bonds.";
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
      UNI10_BLOCK* blk = RQidx2Blk.find(Q_RQoff)->second;
      size_t B_cDim = blk->Cnum;
      size_t sB_cDim = CQidx2Dim.find(Q_CQoff)->second;
      size_t blkRoff = RQidx2Off.find(Q_RQoff)->second;
      size_t blkCoff = CQidx2Off.find(Q_CQoff)->second;
      UNI10_DTYPE* boff = blk->m_elem + (blkRoff * B_cDim) + blkCoff;
      int cnt = 0;
      std::vector<int> D_acc(bondNum, 1);
      for(int b = bondNum	- 1; b > 0; b--)
        D_acc[b - 1] = D_acc[b] * bonds[b].Qdegs[Qidxs[b]];
      for(int b = 0; b < bondNum; b++)
        cnt += (idxs[b] - bonds[b].offsets[Qidxs[b]]) * D_acc[b];
      return boff[(cnt / sB_cDim) * B_cDim + cnt % sB_cDim];
    }
    else{
      return 0.0;;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::vector<size_t>&):");
    return 0;
  }
}

size_t UNI10_TENSOR::blockNum()const{
	return blocks.size();
}

std::vector<Qnum> UNI10_TENSOR::blockQnum()const{
	std::vector<Qnum> keys;
	for(std::map<Qnum, UNI10_BLOCK>::const_iterator it = blocks.begin(); it != blocks.end(); it++)
		keys.push_back(it->first);
	return keys;
}

Qnum UNI10_TENSOR::blockQnum(size_t idx)const{
  try{
    if(!(idx < blocks.size())){
      std::ostringstream err;
      err<<"Index exceeds the number of the blocks("<<blocks.size()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    for(std::map<Qnum, UNI10_BLOCK>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
      if(idx == 0)
        return it->first;
      idx--;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::blockQnum(size_t):");
  }
  return Qnum(0);
}

const std::map<Qnum, UNI10_BLOCK>& UNI10_TENSOR::const_getBlocks()const{
  return blocks;
}

const UNI10_BLOCK& UNI10_TENSOR::const_getBlock()const{
  try{
    Qnum q0(0);
    return const_getBlock(q0);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(bool=false):");
    return blocks.end()->second;
  }
}

const UNI10_BLOCK& UNI10_TENSOR::const_getBlock(const Qnum& qnum)const{
  try{
    std::map<Qnum, UNI10_BLOCK>::const_iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    return it->second;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(uni10::Qnum&):");
    return blocks.end()->second;
  }
}

std::map<Qnum, UNI10_MATRIX> UNI10_TENSOR::getBlocks()const{
	std::map<Qnum, UNI10_MATRIX> mats;
  try{
    for(std::map<Qnum, UNI10_BLOCK>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
      UNI10_MATRIX mat(it->second.Rnum, it->second.Cnum, it->second.m_elem, false, ongpu);
      mats.insert(std::pair<Qnum, UNI10_MATRIX>(it->first, mat));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlocks():");
  }
	return mats;
}

UNI10_MATRIX UNI10_TENSOR::getBlock(bool diag)const{
  try{
    Qnum q0(0);
    return getBlock(q0, diag);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(bool=false):");
    return UNI10_MATRIX();
  }
}

UNI10_MATRIX UNI10_TENSOR::getBlock(const Qnum& qnum, bool diag)const{
  try{
    std::map<Qnum, UNI10_BLOCK>::const_iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag){
      return it->second.getDiag();
    }
    else{
      UNI10_MATRIX mat(it->second.Rnum, it->second.Cnum, it->second.m_elem, false, ongpu);
      return mat;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(uni10::Qnum&):");
    return UNI10_MATRIX(0, 0);
  }
}

void UNI10_TENSOR::putBlock(const UNI10_BLOCK& mat){
  try{
    Qnum q0(0);
    putBlock(q0, mat);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Block&):");
  }
}

void UNI10_TENSOR::putBlock(const Qnum& qnum, const UNI10_BLOCK& mat){
  try{
    std::map<Qnum, UNI10_BLOCK>::iterator it;
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
    if(mat.m_elem != it->second.m_elem){
      if(mat.isDiag()){
        elemBzero(it->second.m_elem, it->second.Rnum * it->second.Cnum * sizeof(UNI10_DTYPE), ongpu);
        setDiag(it->second.m_elem, mat.getElem(), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
      }
      else
        elemCopy(it->second.m_elem, mat.getElem(), it->second.Rnum * it->second.Cnum * sizeof(UNI10_DTYPE), ongpu, mat.isOngpu());
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Qnum&, uni10::Block&):");
  }
}

UNI10_DTYPE* UNI10_TENSOR::getElem(){
  return elem;
}

void UNI10_TENSOR::setElem(const UNI10_DTYPE* _elem, bool _ongpu){
  try{
    elemCopy(elem, _elem, m_elemNum * sizeof(UNI10_DTYPE), ongpu, _ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(double*, bool=false):");
  }
}

void UNI10_TENSOR::setElem(const std::vector<UNI10_DTYPE >& _elem, bool _ongpu){
  try{
    setElem(&_elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(std::vector<double>&, bool=false):");
  }
}

UNI10_DTYPE UNI10_TENSOR::operator[](size_t idx)const{
  try{
    if(!(idx < m_elemNum)){
      std::ostringstream err;
      err<<"Index exceeds the number of elements("<<m_elemNum<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return getElemAt(idx, elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator[](size_t):");
    return 0;
  }
}

void UNI10_TENSOR::set_zero(){
  try{
    elemBzero(elem, m_elemNum * sizeof(UNI10_DTYPE), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero():");
  }
}

void UNI10_TENSOR::set_zero(const Qnum& qnum){
  try{
    std::map<Qnum, UNI10_BLOCK>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    UNI10_BLOCK& block = it->second;
    elemBzero(block.m_elem, block.Rnum * block.Cnum * sizeof(UNI10_DTYPE), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero(std::Qnum&):");
  }
}

void UNI10_TENSOR::identity(const Qnum& qnum){
  try{
    std::map<Qnum, UNI10_BLOCK>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    UNI10_BLOCK& block = it->second;
    setIdentity(block.m_elem, block.Rnum, block.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity(std::Qnum&):");
  }
}

void UNI10_TENSOR::identity(){
	std::map<Qnum, UNI10_BLOCK>::iterator it;
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		setIdentity(it->second.m_elem, it->second.Rnum, it->second.Cnum, ongpu);
	status |= HAVEELEM;
}

void UNI10_TENSOR::randomize(){
	elemRand(elem, m_elemNum, ongpu);
	status |= HAVEELEM;
}

void UNI10_TENSOR::orthoRand(const Qnum& qnum){
  try{
  std::map<Qnum, UNI10_BLOCK>::iterator it = blocks.find(qnum);
  if(it == blocks.end()){
    std::ostringstream err;
    err<<"There is no block with the given quantum number "<<qnum;
    throw std::runtime_error(exception_msg(err.str()));
  }
  UNI10_BLOCK& block = it->second;
	orthoRandomize(block.m_elem, block.Rnum, block.Cnum, ongpu);
	status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand(std::Qnum&):");
  }
}

void UNI10_TENSOR::orthoRand(){
	std::map<Qnum, UNI10_BLOCK>::iterator it;
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		orthoRandomize(it->second.m_elem, it->second.Rnum, it->second.Cnum, ongpu);
	status |= HAVEELEM;
}

void UNI10_TENSOR::clear(){
	status &= ~HAVEELEM;
}

void UNI10_TENSOR::save(const std::string& fname){
  try{
    if((status & HAVEBOND) == 0){   //If not INIT, NO NEED to write out to file
      throw std::runtime_error(exception_msg("Saving a tensor without bonds(scalar) is not supported."));
    }
    FILE* fp = fopen(fname.c_str(), "w");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in writing to file '"<<fname<<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    fwrite(&status, 1, sizeof(status), fp);	//OUT: status(4 bytes)
    int bondNum = bonds.size();
    fwrite(&bondNum, 1, sizeof(bondNum), fp);  //OUT: bondNum(4 bytes)
    size_t qnum_sz = sizeof(Qnum);
    fwrite(&qnum_sz, 1, sizeof(size_t), fp);	//OUT: sizeof(Qnum)
    for(int b = 0; b < bondNum; b++){
      int num_q = bonds[b].Qnums.size();
      fwrite(&(bonds[b].m_type), 1, sizeof(bondType), fp);	//OUT: Number of Qnums in the bond(4 bytes)
      fwrite(&num_q, 1, sizeof(int), fp);	//OUT: Number of Qnums in the bond(4 bytes)
      fwrite(&(bonds[b].Qnums[0]), num_q, qnum_sz, fp);
      fwrite(&(bonds[b].Qdegs[0]), num_q, sizeof(int), fp);
    }
    int num_l = labels.size();
    fwrite(&num_l, 1, sizeof(int), fp);	//OUT: Number of Labels in the Tensor(4 bytes)
    fwrite(&(labels[0]), num_l, sizeof(int), fp);
    if(status & HAVEELEM){
      fwrite(&m_elemNum, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
      size_t memsize = m_elemNum * sizeof(UNI10_DTYPE);
      UNI10_DTYPE* tmp_elem = elem;
      if(ongpu){
        tmp_elem = (UNI10_DTYPE*)malloc(memsize);
        elemCopy(tmp_elem, elem, memsize, false, ongpu);
      }
      fwrite(tmp_elem, m_elemNum, sizeof(UNI10_DTYPE), fp);
      if(ongpu)
        free(tmp_elem);
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::save(std::string&):");
  }
}

UNI10_TENSOR& UNI10_TENSOR::permute(int rowBondNum){
  try{
    std::vector<int> ori_labels = labels;
	  this->permute(ori_labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int):");
  }
	return *this;
}
UNI10_TENSOR& UNI10_TENSOR::permute(int* newLabels, int rowBondNum){
  try{
    std::vector<int> _labels(newLabels, newLabels + bonds.size());
    this->permute(_labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int*, int):");
  }
	return *this;
}

UNI10_TENSOR& UNI10_TENSOR::permute(const std::vector<int>& newLabels, int rowBondNum){
  try{
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
      for(int b = 0; b < bonds.size(); b++){
        outBonds.push_back(bonds[rsp_outin[b]]);
        if(bonds[b].Qnums.size() != 1)
          withoutSymmetry = false;
      }
      for(int b = 0; b < bonds.size(); b++){
        if(b < rowBondNum)
          outBonds[b].change(BD_IN);
        else
          outBonds[b].change(BD_OUT);
      }
      UNI10_TENSOR UniTout(outBonds, name);
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
              UNI10_DTYPE* des_elem = UniTout.elem;
              UNI10_DTYPE* src_elem = elem;
              reshapeElem(src_elem, bondNum, m_elemNum, perInfo, des_elem);
              free(perInfo);
            }
            else{
              UNI10_DTYPE* des_elem = UniTout.elem;
              UNI10_DTYPE* src_elem = elem;
              size_t memsize = m_elemNum * sizeof(UNI10_DTYPE);
              if(ongpu){
                src_elem = (UNI10_DTYPE*)elemAllocForce(memsize, false);
                elemCopy(src_elem, elem, memsize, false, ongpu);
              }
              if(UniTout.ongpu)
                des_elem = (UNI10_DTYPE*)elemAllocForce(memsize, false);

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
              for(int i = 0; i < m_elemNum; i++){
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
                elemCopy(UniTout.elem, des_elem, memsize, UniTout.ongpu, false);
                elemFree(des_elem, memsize, false);
              }
            }
          }
          else{  //non-symmetry inorder
            size_t memsize = m_elemNum * sizeof(UNI10_DTYPE);
            elemCopy(UniTout.elem, elem, memsize, UniTout.ongpu, ongpu);
          }
        }
        else{
          double sign = 1.0;
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
          UNI10_DTYPE* Ein_ptr;
          UNI10_DTYPE* Eot_ptr;
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
            for(int b = bondNum	- 1; b > 0; b--)
              sBot_acc[rsp_outin[b-1]] = sBot_acc[rsp_outin[b]] * bonds[rsp_outin[b]].Qdegs[Qot_idxs[b]];
            Qin_RQoff = Qin_off / CQdim;
            Qin_CQoff = Qin_off % CQdim;
            Qot_RQoff = Qot_off / UniTout.CQdim;
            Qot_CQoff = Qot_off % UniTout.CQdim;
            Bin_cDim = RQidx2Blk[Qin_RQoff]->Cnum;
            Bot_cDim = UniTout.RQidx2Blk[Qot_RQoff]->Cnum;
            Ein_ptr = RQidx2Blk[Qin_RQoff]->m_elem + (RQidx2Off[Qin_RQoff] * Bin_cDim) + CQidx2Off[Qin_CQoff];
            Eot_ptr = UniTout.RQidx2Blk[Qot_RQoff]->m_elem + (UniTout.RQidx2Off[Qot_RQoff] * Bot_cDim) + UniTout.CQidx2Off[Qot_CQoff];
            sBin_rDim = RQidx2Dim[Qin_RQoff];
            sBin_cDim = CQidx2Dim[Qin_CQoff];
            sBot_cDim = UniTout.CQidx2Dim[Qot_CQoff];
            int cnt_ot = 0;
            sBin_idxs.assign(bondNum, 0);
            if(Qnum::isFermionic()){
              int sign01 = 0;
              for(int i = 0; i < swaps.size(); i++)
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
    propogate_exception(e, "In function UniTensor::permute(std::vector<int>&, int):");
  }
  return *this;
}

UNI10_TENSOR& UNI10_TENSOR::transpose(){
  try{
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
    for(int b = 0; b < bonds.size(); b++){
      outBonds.push_back(bonds[rsp_outin[b]]);
      outLabels[b] = labels[rsp_outin[b]];
    }
    for(int b = 0; b < bondNum; b++){
      if(b < cbondNum)
        outBonds[b].m_type = BD_IN;
      else
        outBonds[b].m_type = BD_OUT;
    }
    UNI10_TENSOR UniTout(outBonds, name);
    UniTout.setLabel(outLabels);
    if(status & HAVEELEM){
      std::map<Qnum, UNI10_BLOCK>::iterator it_in;
      std::map<Qnum, UNI10_BLOCK>::iterator it_out;
      UNI10_DTYPE* elem_in;
      UNI10_DTYPE* elem_out;
      size_t Rnum, Cnum;
      for ( it_in = blocks.begin() ; it_in != blocks.end(); it_in++ ){
        it_out = UniTout.blocks.find((it_in->first));
        Rnum = it_in->second.Rnum;
        Cnum = it_in->second.Cnum;
        elem_in = it_in->second.m_elem;
        elem_out = it_out->second.m_elem;
        setTranspose(elem_in, Rnum, Cnum, elem_out, ongpu, UniTout.ongpu);
      }
      UniTout.status |= HAVEELEM;
    }
    *this = UniTout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::transpose():");
  }
  return *this;
}

UNI10_TENSOR& UNI10_TENSOR::combineBond(const std::vector<int>&cmbLabels){
  try{
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
	for(int p = 0; p < cmbLabels.size(); p++){
		for(int l = 0; l < labels.size(); l++){
			if(cmbLabels[p] == labels[l]){
				picked[p] = l;
				marked[l] = 1;
				break;
			}
		}
	}
	int mark = 0;
	for(int m = 0; m < marked.size(); m++)
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
	for(int l = 0; l < labels.size(); l++){
		if(marked[l] && l == picked[0]){
			for(int ll = 0; ll < cmbLabels.size(); ll++){
				rsp_labels[enc] = cmbLabels[ll];
				enc++;
			}
			std::vector<Bond> tmpBonds;
			for(int p = 0; p < picked.size(); p++)
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
	this->permute(rsp_labels, RBnum);
	UNI10_TENSOR Tout(newBonds, reduced_labels);
  if(status & HAVEELEM)
    Tout.setElem(elem, ongpu);
	*this = Tout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::combineBond(std::vector<int>&):");
  }
  return *this;
}

UNI10_TENSOR& UNI10_TENSOR::partialTrace(int la, int lb){
  try{
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
    for(int l = 0; l < labels.size(); l++){
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

    UNI10_TENSOR Tt(newBonds, newLabels);
    rsp_labels[bondNum - 2] = labels[ia];
    rsp_labels[bondNum - 1] = labels[ib];
    ia = bondNum - 2;
    ib = bondNum - 1;
    this->permute(rsp_labels, Tt.RBondNum);
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
    UNI10_DTYPE* Et_ptr;
    std::vector<UNI10_DTYPE*> E_offs(tQdim);
    std::vector<size_t> B_cDims(tQdim);
    int tQdim2 = tQdim * tQdim;
    int Qenc = Q_acc[ia] + Q_acc[ib];
    for(std::map<int, size_t>::iterator it = Tt.QidxEnc.begin(); it != Tt.QidxEnc.end(); it++){
      Qt_off = it->first;
      Qt_RQoff = Qt_off / Tt.CQdim;
      Qt_CQoff = Qt_off % Tt.CQdim;
      Bt_cDim = Tt.RQidx2Blk[Qt_RQoff]->Cnum;
      Et_ptr = Tt.RQidx2Blk[Qt_RQoff]->m_elem + (Tt.RQidx2Off[Qt_RQoff] * Bt_cDim) + Tt.CQidx2Off[Qt_CQoff];
      sBt_rDim = Tt.RQidx2Dim[Qt_RQoff];
      sBt_cDim = Tt.CQidx2Dim[Qt_CQoff];

      for(int q = 0; q < tQdim; q++){
        Q_off = Qt_off * tQdim2 + q * Qenc;
        Q_RQoff = Q_off / CQdim;
        Q_CQoff = Q_off % CQdim;
        B_cDims[q] = RQidx2Blk[Q_RQoff]->Cnum;
        E_offs[q] = RQidx2Blk[Q_RQoff]->m_elem + (RQidx2Off[Q_RQoff] * B_cDims[q]) + CQidx2Off[Q_CQoff];
      }
      int tQdeg, sB_c_off;
      UNI10_DTYPE trVal;
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
    propogate_exception(e, "In function UniTensor::partialTrace(int, int):");
  }
	return *this;
}

UNI10_DTYPE UNI10_TENSOR::trace()const{
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot trace a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(status & HAVEBOND){
      size_t Rnum;
      UNI10_DTYPE trVal = 0;
      for(std::map<Qnum, UNI10_BLOCK>::const_iterator it = blocks.begin() ; it != blocks.end(); it++ ){
        if(!(it->second.Rnum == it->second.Cnum)){
          std::ostringstream err;
          err<<"Cannot trace a non-square block.";
          throw std::runtime_error(exception_msg(err.str()));
        }
        trVal += it->second.trace();
      }
      return trVal;
    }
    else
      return getElemAt(0, elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::trace():");
    return 0;
  }
}

std::vector<UNI10_TENSOR> UNI10_TENSOR::hosvd(size_t modeNum, size_t fixedNum)const{
  try{
    std::vector<std::map<Qnum, UNI10_MATRIX> > symLs;
    return _hosvd(modeNum, fixedNum, symLs, false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t = 0):");
    return std::vector<UNI10_TENSOR>();
  }
}

std::vector<UNI10_TENSOR> UNI10_TENSOR::hosvd(size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, UNI10_MATRIX> >& Ls)const{
  try{
    return _hosvd(modeNum, fixedNum, Ls, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
    return std::vector<UNI10_TENSOR>();
  }
}

std::vector<UNI10_TENSOR> UNI10_TENSOR::hosvd(size_t modeNum, std::vector<std::map<Qnum, UNI10_MATRIX> >& Ls)const{
  try{
    return _hosvd(modeNum, 0, Ls, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
    return std::vector<UNI10_TENSOR>();
  }
}

std::vector<UNI10_TENSOR> UNI10_TENSOR::hosvd(size_t modeNum, size_t fixedNum, std::vector<UNI10_MATRIX>& Ls)const{
  try{
    bool withoutSymmetry = true;
    for(int b = 0; b < bonds.size(); b++){
      if(bonds[b].Qnums.size() != 1)
        withoutSymmetry = false;
    }
    if(!withoutSymmetry){
      std::ostringstream err;
      err<<"The tensor has symmetry quantum numbers. Cannot use non-symmetry version hosvd(size_t, std::vector<uni10::Matrix>&)";
      err<<"\n  Hint: Use UniTensor::hosvd(size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&)";
      throw std::runtime_error(exception_msg(err.str()));
    }
    std::vector<std::map<Qnum, UNI10_MATRIX> > symLs;
    const std::vector<UNI10_TENSOR>& outs = _hosvd(modeNum, fixedNum, symLs, true);
    Ls.clear();
    Qnum q0(0);
    for(int i = 0; i < symLs.size(); i++)
      Ls.push_back(symLs[i][q0]);
    return outs;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t, std::vector<Matrix>&):");
    return std::vector<UNI10_TENSOR>();
  }
}
std::vector<UNI10_TENSOR> UNI10_TENSOR::hosvd(size_t modeNum, std::vector<UNI10_MATRIX>& Ls)const{
  try{
    return hosvd(modeNum, 0, Ls);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, std::vector<Matrix>&):");
    return std::vector<UNI10_TENSOR>();
  }
}

bool UNI10_TENSOR::similar(const UNI10_TENSOR& Tb)const{
  try{
    if(bonds.size() != Tb.bonds.size())
      return false;
    for(int b = 0; b < bonds.size(); b++){
      if(bonds[b] == Tb.bonds[b]);
      else return false;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::similar(uni10::UniTensor&):");
  }
	return true;
}

bool UNI10_TENSOR::elemCmp(const UNI10_TENSOR& UniT)const{
  try{
    double diff;
    if(m_elemNum == UniT.m_elemNum){
      for(size_t i = 0; i < m_elemNum; i++){
        diff = std::abs(elem[i] - UniT.elem[i]);
        if(diff > 1E-12)
          return false;
      }
    }
    else
      return false;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::elemCmp(uni10::UniTensor&):");
  }
	return true;
}

std::string UNI10_TENSOR::printRawElem(bool print)const{
  try{
    std::ostringstream os;
    if(status & HAVEBOND && status & HAVEELEM){
      int bondNum = bonds.size();
      std::vector<Bond> ins;
      std::vector<Bond> outs;
      for(std::vector<Bond>::const_iterator it = bonds.begin(); it != bonds.end(); ++it){
        if(it->type() == BD_IN)
          ins.push_back(*it);
        else
          outs.push_back(*it);
      }
      if(ins.size() == 0 || outs.size() == 0)
        os<<getRawElem();
      else{
        Bond rBond = combine(ins);
        Bond cBond = combine(outs);
        std::vector<Qnum> rowQ = rBond.Qlist();
        std::vector<Qnum> colQ = cBond.Qlist();
        size_t rowNum = rBond.dim();
        size_t colNum = cBond.dim();
        std::vector<size_t> idxs(bondNum, 0);

        os<< "     ";
        for(int q = 0; q < colQ.size(); q++)
          os<< "   " << std::setw(2) << colQ[q].U1() << "," << colQ[q].prt();
        os<< std::endl << std::setw(5) << "" << std::setw(colQ.size() * 7 + 2) <<std::setfill('-')<<"";
        os<<std::setfill(' ');
        int cnt = 0;
        int r = 0;
        int bend;
        while(1){
          if(cnt % colNum == 0){
            os<<"\n    |\n" << std::setw(2) << rowQ[r].U1() << "," << rowQ[r].prt() << "|";
            r++;
          }
          os<< std::setw(7) << std::fixed << std::setprecision(3) << at(idxs);
          for(bend = bondNum - 1; bend >= 0; bend--){
            idxs[bend]++;
            if(idxs[bend] < bonds[bend].dim())
              break;
            else
              idxs[bend] = 0;
          }
          cnt++;
          if(bend < 0)
            break;
        }
        os <<"\n    |\n";
      }
    }
    else if(status & HAVEELEM){
      os<<"\nScalar: " << elem[0]<<"\n\n";
    }
    else{
      os<<"NO ELEMENT IN THE TENSOR!!!\n";
    }
    if(print){
      std::cout<<os.str();
      return "";
    }
    return os.str();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::printRawElem():");
    return "";
  }
}

std::string UNI10_TENSOR::profile(bool print){
  std::ostringstream os;
  os<<"\n===== Tensor profile =====\n";
	os<<"Existing Tensors: " << COUNTER << std::endl;
	os<<"Allocated Elements: " << ELEMNUM << std::endl;
	os<<"Max Allocated Elements: " << MAXELEMNUM << std::endl;
	os<<"Max Allocated Elements for a Tensor: " << MAXELEMTEN << std::endl;
  os<<"============================\n\n";
  if(print){
    std::cout<<os.str();
    return "";
  }
  return os.str();
}

UNI10_TENSOR& UNI10_TENSOR::operator*= (double a){
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    vectorScal(a, elem, m_elemNum, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator*=(double):");
  }
	return *this;
}

UNI10_TENSOR& UNI10_TENSOR::operator*=(const UNI10_TENSOR& uT){
  try{
	  *this = *this * uT;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator*=(uni10::UniTensor&):");
  }
  return *this;
}

UNI10_TENSOR& UNI10_TENSOR::operator+= (const UNI10_TENSOR& Tb){
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
    propogate_exception(e, "In function UniTensor::operator+=(uni10::UniTensor&):");
  }
  return *this;
}
std::vector<_Swap> UNI10_TENSOR::exSwap(const UNI10_TENSOR& Tb) const{
	std::vector<_Swap> swaps;
  try{
    if(status & Tb.status & HAVEBOND){
      int bondNumA = labels.size();
      int bondNumB = Tb.labels.size();
      std::vector<int> intersect;
      std::vector<int> left;
      for(int a = 0; a < bondNumA; a++){
        bool found = false;
        for(int b = 0; b < bondNumB; b++)
          if(labels[a] == Tb.labels[b])
            found = true;
        if(found)
          intersect.push_back(a);
        else
          left.push_back(a);
      }
      _Swap sp;
      for(int i = 0; i < intersect.size(); i++)
        for(int j = 0; j < left.size(); j++){
          sp.b1 = intersect[i];
          sp.b2 = left[j];
          swaps.push_back(sp);
        }
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::exSwap(uni10::UniTensor&):");
  }
	return swaps;
}

void UNI10_TENSOR::addGate(const std::vector<_Swap>& swaps){
  try{
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
	UNI10_DTYPE* Eptr;
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
		Eptr = RQidx2Blk[RQoff]->m_elem + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
		sB_rDim = RQidx2Dim[RQoff];
		sB_cDim = CQidx2Dim[CQoff];

		int sign01 = 0;
		for(int i = 0; i < swaps.size(); i++)
			sign01 ^= (bonds[swaps[i].b1].Qnums[Q_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Q_idxs[swaps[i].b2]].prtF());
		sign = sign01 ? -1 : 1;

		for(sB_r = 0; sB_r < sB_rDim; sB_r++)
			for(sB_c = 0; sB_c < sB_cDim; sB_c++)
				Eptr[(sB_r * B_cDim) + sB_c] *= sign;
	}
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::addGate(std::vector<_Swap>&):");
  }
}

UNI10_TENSOR contract(UNI10_TENSOR& Ta, UNI10_TENSOR& Tb, bool fast){
  try{
    if(!(Ta.status & Tb.status & Ta.HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform contraction of two tensors before setting their elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(&Ta == &Tb){
      UNI10_TENSOR Ttmp = Tb;
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
      for(int a = 0; a < interLabel.size(); a++)
        newLabelA.push_back(interLabel[a]);
      for(int b = 0; b < BbondNum; b++)
        if(markB[b] == 0){
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
      UNI10_TENSOR Tc(cBonds);
      if(cBonds.size())
        Tc.setLabel(newLabelC);
      UNI10_BLOCK blockA, blockB, blockC;
      std::map<Qnum, UNI10_BLOCK>::iterator it;
      std::map<Qnum, UNI10_BLOCK>::iterator it2;
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
          matrixMul(blockA.getElem(), blockB.getElem(), blockA.row(), blockB.col(), blockA.col(), blockC.getElem(), Ta.ongpu, Tb.ongpu, Tc.ongpu);
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
      return UNI10_TENSOR(Ta[0] * Tb[0]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function contract(uni10::UniTensor&, uni10::UniTensor, bool):");
    return UNI10_TENSOR();
  }
}

UNI10_TENSOR otimes(const UNI10_TENSOR & Ta, const UNI10_TENSOR& Tb){
  try{
    UNI10_TENSOR T1 = Ta;
    UNI10_TENSOR T2 = Tb;
    std::vector<int> label1(T1.bondNum());
    std::vector<int> label2(T2.bondNum());
    for(int i = 0; i < T1.bondNum(); i++){
      if(i < T1.inBondNum())
        label1[i] = i;
      else
        label1[i] = T2.inBondNum() + i;
    }
    for(int i = 0; i < T2.bondNum(); i++){
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
    return UNI10_TENSOR();
  }
}

UNI10_MATRIX otimes(const UNI10_BLOCK& Ma, const UNI10_BLOCK& Mb){
  try{
    UNI10_TENSOR Ta(Ma);
    UNI10_TENSOR Tb(Mb);
    return otimes(Ta, Tb).getBlock();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function otimes(uni10::Matrix&, uni10::Matrix&):");
    return UNI10_MATRIX();
  }
}

UNI10_TENSOR operator*(const UNI10_TENSOR& Ta, const UNI10_TENSOR& Tb){
  try{
    UNI10_TENSOR cTa = Ta;
    UNI10_TENSOR cTb = Tb;
    return contract(cTa, cTb, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::UniTensor&, uni10::UniTensor&):");
    return UNI10_TENSOR();
  }
}

UNI10_TENSOR operator*(const UNI10_TENSOR& Ta, double a){
  try{
    if(!(Ta.status & Ta.HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UNI10_TENSOR Tb(Ta);
    vectorScal(a, Tb.elem, Tb.m_elemNum, Tb.ongpu);
    return Tb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::UniTensor&, double):");
    return UNI10_TENSOR();
  }
}

UNI10_TENSOR operator*(double a, const UNI10_TENSOR& Ta){return Ta * a;};

UNI10_TENSOR operator+(const UNI10_TENSOR& Ta, const UNI10_TENSOR& Tb){
  try{
    if(!(Ta.status & Tb.status & Ta.HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform addition of tensors before setting their elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(!(Ta.bonds == Tb.bonds)){
      std::ostringstream err;
      err<<"Cannot perform addition of two tensors having diffirent bonds.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UNI10_TENSOR Tc(Ta);
    vectorAdd(Tc.elem, Tb.elem, Tc.m_elemNum, Tc.ongpu, Tb.ongpu);
    return Tc;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator+(uni10::UniTensor&, uni10::UniTensor&):");
    return UNI10_TENSOR();
  }
}

std::ostream& operator<< (std::ostream& os, const UNI10_TENSOR& UniT){
  try{
    if(!(UniT.status & UniT.HAVEBOND)){
      if(UniT.ongpu){
        os<<"\nScalar: " << getElemAt(0, UniT.elem, UniT.ongpu);
        os<<", onGPU";
      }
      else
        os<<"\nScalar: " << UniT.elem[0];
      os<<"\n\n";
      return os;
    }
    int row = 0;
    int col = 0;
    std::vector<Bond>bonds = UniT.bond();
    for(int i = 0; i < bonds.size(); i++)
      if(bonds[i].type() == BD_IN)
        row++;
      else
        col++;
    int layer = std::max(row, col);
    int nmlen = UniT.name.length() + 2;
    int star = 12 + (14 - nmlen) / 2;
    os<<std::endl;
    for(int s = 0; s < star; s++)
      os << "*";
    if(UniT.name.length() > 0)
      os << " " << UniT.name << " ";
    for(int s = 0; s < star; s++)
      os<<"*";
    if(UniT.ongpu)
      os<<"\n                 onGPU";
    os << "\n             ____________\n";
    os << "            |            |\n";
    int llab = 0;
    int rlab = 0;
    char buf[128];
    for(int l = 0; l < layer; l++){
      if(l < row && l < col){
        llab = UniT.labels[l];
        rlab = UniT.labels[row + l];
        sprintf(buf, "    %5d___|%-4d    %4d|___%-5d\n", llab, bonds[l].dim(), bonds[row + l].dim(), rlab);
        os<<buf;
      }
      else if(l < row){
        llab = UniT.labels[l];
        sprintf(buf, "    %5d___|%-4d    %4s|\n", llab, bonds[l].dim(), "");
        os<<buf;
      }
      else if(l < col){
        rlab = UniT.labels[row + l];
        sprintf(buf, "    %5s   |%4s    %4d|___%-5d\n", "", "", bonds[row + l].dim(), rlab);
        os << buf;
      }
      os << "            |            |   \n";
    }
    os << "            |____________|\n";

    os << "\n================BONDS===============\n";
    for(int b = 0; b < bonds.size(); b++){
      os << bonds[b];
    }
    os<<"\n===============BLOCKS===============\n";
    std::map<Qnum, UNI10_BLOCK> blocks = UniT.const_getBlocks();
    bool printElem = true;
    for (std::map<Qnum, UNI10_BLOCK>::const_iterator  it = blocks.begin() ; it != blocks.end(); it++ ){
      os << "--- " << it->first << ": ";// << Rnum << " x " << Cnum << " = " << Rnum * Cnum << " ---\n\n";
      if((UniT.status & UniT.HAVEELEM) && printElem)
        os<<it->second;
      else
        os<<it->second.row() << " x "<<it->second.col()<<": "<<it->second.elemNum()<<std::endl<<std::endl;
    }
    os << "Total elemNum: "<<UniT.m_elemNum<<std::endl;
    os << "***************** END ****************\n\n";
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator<<(std::ostream&, uni10::UniTensor&):");
  }
	return os;
}


/**************** Private Functions ***********************/

void UNI10_TENSOR::initUniT(){ //GPU
	if(bonds.size()){
		m_elemNum = grouping();
		if(!(blocks.size() > 0)){ //No block in Tensor, Error!
      std::ostringstream err;
      err<<"There is no symmetry block with the given bonds:\n";
      for(int b = 0; b < bonds.size(); b++)
        err<<"    "<<bonds[b];
      throw std::runtime_error(exception_msg(err.str()));
    }
		labels.assign(bonds.size(), 0);
		for(int b = 0; b < bonds.size(); b++)
			labels[b] = b;
		status |= HAVEBOND;
	}
	else{
		Qnum q0(0);
		blocks[q0] = UNI10_BLOCK(1, 1);
		RBondNum = 0;
		RQdim = 0;
		CQdim = 0;
		m_elemNum = 1;
		status |= HAVEELEM;
	}
	elem = NULL;
	elem = (UNI10_DTYPE*)elemAlloc(sizeof(UNI10_DTYPE) * m_elemNum, ongpu);
  size_t offset = 0;
	for (std::map<Qnum, UNI10_BLOCK>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
		it->second.m_elem = &(elem[offset]);
    it->second.ongpu = ongpu;
    offset += it->second.Rnum * it->second.Cnum;
  }
	elemBzero(elem, sizeof(UNI10_DTYPE) * m_elemNum, ongpu);
	ELEMNUM += m_elemNum;
	COUNTER++;
	if(ELEMNUM > MAXELEMNUM)
		MAXELEMNUM = ELEMNUM;
	if(m_elemNum > MAXELEMTEN)
		MAXELEMTEN = m_elemNum;
}

size_t UNI10_TENSOR::grouping(){
	blocks.clear();
	int row_bondNum = 0;
	int col_bondNum = 0;
	RQdim = 1;
	CQdim = 1;
        bool IN_BONDS_BEFORE_OUT_BONDS = true;
	for(int i = 0; i < bonds.size(); i++){
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
    UNI10_BLOCK blk(it->second, it2->second); // blk(Rnum, Cnum);
		off += blk.Rnum * blk.Cnum;
		blocks[it->first] = blk;
	  UNI10_BLOCK* blkptr = &(blocks[it->first]);
		std::vector<int>& tmpRQidx = row_Qnum2Qidx[it->first];
		std::vector<int>& tmpCQidx = col_Qnum2Qidx[it->first];
		for(int i = 0; i < tmpRQidx.size(); i++){
			RQidx2Blk[tmpRQidx[i]] = blkptr;
			for(int j = 0; j < tmpCQidx.size(); j++){
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

std::vector<UNI10_TENSOR> UNI10_TENSOR::_hosvd(size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, UNI10_MATRIX> >& Ls, bool returnL)const{
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
  int bondNum = bonds.size();
  if((bondNum - fixedNum) % modeNum != 0){
    std::ostringstream err;
    err<<"Bond number cannot be divided by the input mode number("<<modeNum<<").";
    throw std::runtime_error(exception_msg(err.str()));
  }
  int combNum = (bondNum - fixedNum) / modeNum;
  UNI10_TENSOR T(*this);
  for(int t = 0; t < T.labels.size(); t++)
    T.labels[t] = t;
  std::vector<int>ori_labels = T.labels;
  std::vector<int>rsp_labels = T.labels;
  if(returnL)
    Ls.assign(modeNum, std::map<Qnum, UNI10_MATRIX>());
  std::vector<UNI10_TENSOR> Us;
  UNI10_TENSOR S(T);
  std::vector<int>out_labels(S.labels.begin(), S.labels.begin() + fixedNum + modeNum);
  for(int m = 0; m < modeNum; m++){
    for(int l = 0; l < rsp_labels.size(); l++){
      if(l < combNum)
        rsp_labels[l] = ori_labels[fixedNum + (((m) * combNum + l) % (bondNum - fixedNum))];
      else if(l < combNum + fixedNum)
        rsp_labels[l] = ori_labels[l - combNum];
      else
        rsp_labels[l] = ori_labels[fixedNum + (((m) * combNum + l - fixedNum) % (bondNum - fixedNum))];
    }
    T.permute(rsp_labels, combNum);
    std::vector<Bond> bonds(T.bonds.begin(), T.bonds.begin() + combNum);
    bonds.push_back(combine(bonds).dummy_change(BD_OUT));
    Us.push_back(UNI10_TENSOR(bonds));
    for(std::map<Qnum, UNI10_BLOCK>::iterator it = T.blocks.begin(); it != T.blocks.end(); it++){
      std::vector<UNI10_MATRIX> svd = it->second.svd();
      Us[m].putBlock(it->first, svd[0]);
      if(returnL)
        Ls[m][it->first] = svd[1];
    }
    for(int c = 0; c < combNum; c++)
      Us[m].labels[c] = fixedNum + m * combNum + c;
    Us[m].labels[combNum] = -m - 1;
    out_labels[fixedNum + m] = -m -1;
    UNI10_TENSOR UT = Us[m];
    S *= UT.transpose();
  }
  S.permute(out_labels, fixedNum);
  Us.push_back(S);
  return Us;
}
}; /* namespace uni10 */
#ifdef UNI10_TENSOR
#undef UNI10_TENSOR
#endif
#ifdef UNI10_BLOCK
#undef UNI10_BLOCK
#endif
#ifdef UNI10_MATRIX
#undef UNI10_MATRIX
#endif
#ifdef UNI10_DTYPE
#undef UNI10_DTYPE
#endif

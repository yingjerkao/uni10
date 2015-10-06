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
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/tensor-network/CMatrix.h>
#include <uni10/tensor-network/CUniTensor.h>
#include <uni10/tensor-network/UniTensor.h>

typedef double Real;
typedef std::complex<double> Complex;
namespace uni10{

int64_t UniTensor::ELEMNUM = 0;
int UniTensor::COUNTER = 0;
size_t UniTensor::MAXELEMNUM = 0;
size_t UniTensor::MAXELEMTEN = 0;

/*********************  OPERATOR **************************/	    

std::ostream& operator<< (std::ostream& os, const UniTensor& UniT){
  try{
    if(!(UniT.status & UniT.HAVEBOND)){
      if(UniT.ongpu){
        if(UniT.typeID() == 0)
          os<<"This Tensor is EMPTY";
        if(UniT.typeID() == 1)
          os<<"\nScalar: " << getElemAt(0, UniT.elem, UniT.ongpu);
        if(UniT.typeID() == 2)
          os<<"\nScalar: " << getElemAt(0, UniT.c_elem, UniT.ongpu);
        os<<", onGPU";
      }
      else{
        if(UniT.typeID() == 0)
          os<<"This Tensor is EMPTY";
        if(UniT.typeID() == 1)
          os<<"\nScalar: " << UniT.elem[0];
        if(UniT.typeID() == 2)
          os<<"\nScalar: " << UniT.c_elem[0];
      }
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
    for(int s = 0; s < star; s++)
      os << "*";
    if(UniT.name.length() > 0)
      os << " " << UniT.name << " ";
    for(int s = 0; s < star; s++)
      os<<"*";
    os<<std::endl;
    if(UniT.typeID() == 0)
      os << "EMPTY" << std::endl; 
    if(UniT.typeID() == 1)
      os << "REAL" << std::endl; 
    if(UniT.typeID() == 2)
      os << "COMPLEX" << std::endl;
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
    os<<"\n===============BlockS===============\n";
    std::map<Qnum, Block> blocks = UniT.const_getBlocks();
    bool printElem = true;
    for (std::map<Qnum, Block>::const_iterator  it = blocks.begin() ; it != blocks.end(); it++ ){
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

UniTensor& UniTensor::operator*=(const UniTensor& uT){
  try{
    *this = *this * uT;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator*=(uni10::UniTensor&):");
  }
  return *this;
}

UniTensor& UniTensor::operator*= (double a){
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    *this = *this * a;
//    vectorScal(a, elem, m_elemNum, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator*=(double):");
  }
  return *this;
}

UniTensor& UniTensor::operator*= (Complex a){
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    *this = *this * a;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator*=(std::complex<double>):");
  }
  return *this;
}

UniTensor operator*(const std::complex<double>& a, const UniTensor& Ta){
  try{
    if(a.imag() == 0)
      return a.real()*Ta;
    if(!(Ta.status & Ta.HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UniTensor Tb(Ta);
    if(Tb.typeID() == 1)
      RtoC(Tb);
    vectorScal(a, Tb.c_elem, Tb.m_elemNum, Tb.ongpu);
    return Tb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::UniTensor&, complex<double>&):");
    return UniTensor();
  }
}
        
UniTensor operator*(const UniTensor& Ta, const std::complex<double>& a){return a * Ta;};

UniTensor operator*(const UniTensor& Ta, double a){
  try{
    if(!(Ta.status & Ta.HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot perform scalar multiplication on a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UniTensor Tb(Ta);
    if(Tb.typeID() == 1)
      vectorScal(a, Tb.elem, Tb.m_elemNum, Tb.ongpu);
    else if(Tb.typeID() == 2)
      vectorScal(a, Tb.c_elem, Tb.m_elemNum, Tb.ongpu);
    return Tb;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::UniTensor&, double):");
    return UniTensor();
  }
}

UniTensor operator*(double a, const UniTensor& Ta){return Ta * a;};

UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb){
  try{
    UniTensor cTa = Ta;
    UniTensor cTb = Tb;
    return contract(cTa, cTb, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator*(uni10::UniTensor&, uni10::UniTensor&):");
    return UniTensor();
  }
}

UniTensor& UniTensor::operator=(const UniTensor& UniT){ //GPU
  try{
    
    r_flag = UniT.r_flag;
    c_flag = UniT.c_flag;
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
    
    
    uelemFree(); 
    
    elem = NULL;
    c_elem = NULL;
    status = UniT.status;
    m_elemNum = UniT.m_elemNum;
 

    std::map<Qnum, Block>::const_iterator it2;
    std::map< const Block* , Block*> blkmap;
    
    if(typeID() == 0){
      for (std::map<Qnum, Block>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
        it2 = UniT.blocks.find(it->first);
        blkmap[(&(it2->second))] = &(it->second);
      }
    }else if(typeID() == 1){
      uelemAlloc(RTYPE);
      for (std::map<Qnum, Block>::iterator it = blocks.begin(); it != blocks.end(); it++ ){ // blocks here is UniT.blocks
        it->second.m_elem = &(elem[it->second.m_elem - UniT.elem]);

        it2 = UniT.blocks.find(it->first);
        blkmap[&(it2->second)] = &(it->second);
      }
    }else if(typeID() == 2){
      uelemAlloc(CTYPE);
      for (std::map<Qnum, Block>::iterator it = blocks.begin(); it != blocks.end(); it++ ){ // blocks here is UniT.blocks
        it->second.cm_elem = &(c_elem[it->second.cm_elem - UniT.c_elem]);

        it2 = UniT.blocks.find(it->first);
        blkmap[&(it2->second)] = &(it->second);
      }
    }

    if(UniT.status & HAVEBOND){
      for(std::map<int, Block*>::iterator it = RQidx2Blk.begin(); it != RQidx2Blk.end(); it++)
        it->second = blkmap[it->second];
    }

    ELEMNUM += m_elemNum;
    if(ELEMNUM > MAXELEMNUM)
      MAXELEMNUM = ELEMNUM;
    if(m_elemNum > MAXELEMTEN)
      MAXELEMTEN = m_elemNum;

    if(typeID() == 1)
      elemCopy(elem, UniT.elem, sizeof(Real) * UniT.m_elemNum, ongpu, UniT.ongpu);
    else if(typeID() == 2)
      elemCopy(c_elem, UniT.c_elem, sizeof(Complex) * UniT.m_elemNum, ongpu, UniT.ongpu);

  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator=(uni10::UniTensor&):");
  }
  return *this;
}

UniTensor operator+(const UniTensor& _Ta, const UniTensor& _Tb){
  try{
    
    UniTensor Ta(_Ta);
    UniTensor Tb(_Tb);
    if(Ta.typeID() == 1 && Tb.typeID() == 2)
      RtoC(Ta);
    if(Ta.typeID() == 2 && Tb.typeID() == 1)
      RtoC(Tb);
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
    
    UniTensor Tc(Ta);
    Tc.typeID() == 1 ? vectorAdd(Tc.elem, Tb.elem, Tc.m_elemNum, Tc.ongpu, Tb.ongpu) :vectorAdd(Tc.c_elem, Tb.c_elem, Tc.m_elemNum, Tc.ongpu, Tb.ongpu);
    return Tc;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator+(uni10::UniTensor&, uni10::UniTensor&):");
    return UniTensor();
  }
}

UniTensor& UniTensor::operator+= (const UniTensor& _Tb){
  try{
    
    UniTensor Tb(_Tb);
    if(typeID() == 1 && Tb.typeID() == 2)
      RtoC(*this);
    if(typeID() == 2 && Tb.typeID() == 1)
      RtoC(Tb);
    
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
    
    typeID() == 1 ? vectorAdd(elem, Tb.elem, m_elemNum, ongpu, Tb.ongpu) : vectorAdd(c_elem, Tb.c_elem, m_elemNum, ongpu, Tb.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator+=(uni10::UniTensor&):");
  }
  return *this;
}

Matrix otimes(const Block& Ma, const Block& Mb){
  try{
    UniTensor Ta(Ma);
    UniTensor Tb(Mb);
    return otimes(Ta, Tb).getBlock();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function otimes(uni10::Matrix&, uni10::Matrix&):");
    return Matrix();
  }
}


/*********************  NO TYPE **************************/	    

UniTensor::UniTensor(): status(0){
  try{
    initUniT();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor():");
  }
}


UniTensor::UniTensor(const std::vector<Bond>& _bonds, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT();
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::vector<int>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(const std::vector<Bond>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT();
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, int*, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(const UniTensor& UniT): //GPU
  r_flag(UniT.r_flag), c_flag(UniT.c_flag),status(UniT.status), bonds(UniT.bonds), blocks(UniT.blocks), labels(UniT.labels), name(UniT.name),
  RBondNum(UniT.RBondNum), RQdim(UniT.RQdim), CQdim(UniT.CQdim), m_elemNum(UniT.m_elemNum), elem(NULL), c_elem(NULL), 
  QidxEnc(UniT.QidxEnc), RQidx2Off(UniT.RQidx2Off), CQidx2Off(UniT.CQidx2Off), RQidx2Dim(UniT.RQidx2Dim), RQidx2Blk(UniT.RQidx2Blk), CQidx2Dim(UniT.CQidx2Dim){
    try{
      
      std::map<Qnum, Block>::const_iterator it2;
      std::map<const Block* , Block*> blkmap;
      
      if(typeID() == 0){
        for (std::map<Qnum, Block>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
          it2 = UniT.blocks.find(it->first);
          blkmap[(&(it2->second))] = &(it->second);
        }
      } 
      
      if(typeID() == 1){
        uelemAlloc(RTYPE);
        for (std::map<Qnum, Block>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
          it->second.m_elem = &(elem[(it->second.m_elem - UniT.elem)]);
          it2 = UniT.blocks.find(it->first);
          blkmap[(&(it2->second))] = &(it->second);
        }
      } 

      if(typeID() == 2){
        uelemAlloc(CTYPE);
        for (std::map<Qnum, Block>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
          it->second.cm_elem = &(c_elem[(it->second.cm_elem - UniT.c_elem)]);

          it2 = UniT.blocks.find(it->first);
          blkmap[(&(it2->second))] = &(it->second);
        }
      } 

      if(UniT.status & HAVEBOND){
        for(std::map<int, Block*>::iterator it = RQidx2Blk.begin(); it != RQidx2Blk.end(); it++)
          it->second = blkmap[it->second];
      }
      
      ELEMNUM += m_elemNum;
      COUNTER++;
      if(ELEMNUM > MAXELEMNUM)
        MAXELEMNUM = ELEMNUM;
      if(m_elemNum > MAXELEMTEN)
        MAXELEMTEN = m_elemNum;
      
      if(typeID() == 1)
        elemCopy(elem, UniT.elem, sizeof(Real) * UniT.m_elemNum, ongpu, UniT.ongpu);
      if(typeID() == 2)
        elemCopy(c_elem, UniT.c_elem, sizeof(Complex) * UniT.m_elemNum, ongpu, UniT.ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In copy constructor UniTensor::UniTensor(uni10::UniTensor&):");
    }
  }

UniTensor::UniTensor(const Block& blk): status(0){
  try{
    Bond bdi(BD_IN, blk.Rnum);
    Bond bdo(BD_OUT, blk.Cnum);
    bonds.push_back(bdi);
    bonds.push_back(bdo);
    initUniT(blk.typeID());
    this->putBlock(blk);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(uni10::Block&");
  }
}

UniTensor::~UniTensor(){
  uelemFree();
  ELEMNUM -= m_elemNum;
  COUNTER--;
}

UniTensor::UniTensor(const std::string& fname): status(0){ //GPU
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
    fread(&r_flag, 1, sizeof(r_flag), fp);
    fread(&c_flag, 1, sizeof(c_flag), fp);
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
    initUniT(typeID());
    int num_l;
    fread(&num_l, 1, sizeof(int), fp);	//OUT: Number of Labels in the Tensor(4 bytes)
    if(!(num_l == bonds.size())){
      std::ostringstream err;
      err<<"Error in reading file '"<<fname<<"' in.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    labels.assign(num_l, 0);
    fread(&(labels[0]), num_l, sizeof(int), fp);
    if(typeID() == 1){
      if(st & HAVEELEM){
        Real *tmp_elem = elem;
        size_t memsize = m_elemNum * sizeof(Real);
        if(ongpu)
          tmp_elem = (Real*)malloc(memsize);
        size_t num_el;
        fread(&num_el, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
        fread(tmp_elem, m_elemNum, sizeof(Real), fp);
        if(ongpu){
          elemCopy(elem, tmp_elem, memsize, ongpu, false);
          free(tmp_elem);
        }
        status |= HAVEELEM;
      }
    }
    if(typeID() == 2){
      if(st & HAVEELEM){
        Complex *tmp_elem = c_elem;
        size_t memsize = m_elemNum * sizeof(Complex);
        if(ongpu)
          tmp_elem = (Complex*)malloc(memsize);
        size_t num_el;
        fread(&num_el, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
        fread(tmp_elem, m_elemNum, sizeof(Complex), fp);
        if(ongpu){
          elemCopy(c_elem, tmp_elem, memsize, ongpu, false);
          free(tmp_elem);
        }
        status |= HAVEELEM;
      }
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::string&):");
  }
}


int UniTensor::typeID()const{
  return r_flag + c_flag;
}

void UniTensor::setLabel(const std::vector<int>& newLabels){
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

void UniTensor::setLabel(int* newLabels){
  try{
    std::vector<int> labels(newLabels, newLabels + bonds.size());
    setLabel(labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setLabel(int*):");
  }
}

std::vector<int> UniTensor::label()const{
  return labels;
}

int UniTensor::label(size_t idx)const{
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

std::string UniTensor::getName(){
  return name;
}

void UniTensor::setName(const std::string& _name){
  name = _name;
}

size_t UniTensor::bondNum()const{return bonds.size();}

size_t UniTensor::inBondNum()const{return RBondNum;}

std::vector<Bond> UniTensor::bond()const{
  return bonds;
}

Bond UniTensor::bond(size_t idx)const{
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

size_t UniTensor::elemNum()const{return m_elemNum;}

size_t UniTensor::blockNum()const{
	return blocks.size();
}

std::vector<Qnum> UniTensor::blockQnum()const{
  std::vector<Qnum> keys;
  for(std::map<Qnum, Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++)
    keys.push_back(it->first);
  return keys;
}

Qnum UniTensor::blockQnum(size_t idx)const{
  try{
    if(!(idx < blocks.size())){
      std::ostringstream err;
      err<<"Index exceeds the number of the blocks("<<blocks.size()<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    for(std::map<Qnum, Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
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
const std::map<Qnum, Block>& UniTensor::const_getBlocks()const{
  return blocks;
}

const Block& UniTensor::const_getBlock()const{
  try{
    Qnum q0(0);
    return const_getBlock(q0);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(bool=false):");
    return blocks.end()->second;
  }
}

const Block& UniTensor::const_getBlock(const Qnum& qnum)const{
  try{
    std::map<Qnum, Block>::const_iterator it = blocks.find(qnum);
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


std::map<Qnum, Matrix> UniTensor::getBlocks()const{
  try{
    if(typeID() == 1)
      return getBlocks(RTYPE);
    else if(typeID() == 2)
      return getBlocks(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"This tensor is EMPTY ";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlocks():");
  }
}

Matrix UniTensor::getBlock(bool diag)const{
  try{
    if(typeID() == 1)
      return getBlock(RTYPE, diag);
    else if(typeID() == 2)
      return getBlock(CTYPE, diag);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"This tensor is EMPTY ";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(bool=false):");
    return Matrix();
  }
}

Matrix UniTensor::getBlock(const Qnum& qnum, bool diag)const{
  try{
    if(typeID() == 1)
      return getBlock(RTYPE, qnum, diag);
    else if(typeID() == 2)
      return getBlock(CTYPE, qnum, diag);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"This tensor is EMPTY ";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(uni10::Qnum&):");
    return Matrix(0, 0);
  }
}

void UniTensor::set_zero(){
  try{
    if(typeID() == 1)         
      set_zero(RTYPE);
    else if(typeID() == 2)         
      set_zero(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set zeros in EMPTY tensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero():");
  }
}

void UniTensor::set_zero(const Qnum& qnum){
  try{
    if(typeID() == 1)         
      set_zero(RTYPE, qnum);
    else if(typeID() == 2)         
      set_zero(CTYPE, qnum);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set zeros in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero(std::Qnum&):");
  }
}

void UniTensor::identity(){
  try{
    if(typeID() == 1)         
      identity(RTYPE);
    else if(typeID() == 2)         
      identity(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set identity in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity():");
  }
}

void UniTensor::identity(const Qnum& qnum){
  try{
    if(typeID() == 1)         
      identity(RTYPE, qnum);
    else if(typeID() == 2)         
      identity(CTYPE, qnum);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set identity in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity(std::Qnum&):");
  }
}

void UniTensor::randomize(){
  try{
    if(typeID() == 1)         
      randomize(RTYPE);
    else if(typeID() == 2)         
      randomize(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set randomize in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::randomize():");
  }
}

void UniTensor::orthoRand(){
  try{
    if(typeID() == 1)         
      orthoRand(RTYPE);
    else if(typeID() == 2)         
      orthoRand(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set orthoRand in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand():");
  }
}

void UniTensor::orthoRand(const Qnum& qnum){
  try{
    if(typeID() == 1)         
      orthoRand(RTYPE, qnum);
    else if(typeID() == 2)         
      orthoRand(CTYPE, qnum);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set orthoRand in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand(std::Qnum&):");
  }
}

void UniTensor::save(const std::string& fname){
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
    fwrite(&r_flag, 1, sizeof(r_flag), fp);
    fwrite(&c_flag, 1, sizeof(c_flag), fp);
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
    if(typeID() == 1){
      if(status & HAVEELEM){
        fwrite(&m_elemNum, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
        size_t memsize = m_elemNum * sizeof(Real);
        Real* tmp_elem = elem;
        if(ongpu){
          tmp_elem = (Real*)malloc(memsize);
          elemCopy(tmp_elem, elem, memsize, false, ongpu);
        }
        fwrite(tmp_elem, m_elemNum, sizeof(Real), fp);
        if(ongpu)
          free(tmp_elem);
      }
    }
    if(typeID() == 2){
      if(status & HAVEELEM){
        fwrite(&m_elemNum, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
        size_t memsize = m_elemNum * sizeof(Complex);
        Complex* tmp_elem = c_elem;
        if(ongpu){
          tmp_elem = (Complex*)malloc(memsize);
          elemCopy(tmp_elem, c_elem, memsize, false, ongpu);
        }
        fwrite(tmp_elem, m_elemNum, sizeof(Complex), fp);
        if(ongpu)
          free(tmp_elem);
      }
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::save(std::string&):");
  }
}


UniTensor& UniTensor::transpose(){
  try{
    if(typeID() == 1)         
      return transpose(RTYPE);
    else if(typeID() == 2)         
      return transpose(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't set orthoRand in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::transpose():");
  }
}

UniTensor& UniTensor::permute(int rowBondNum){
  try{
    if(typeID() == 1)         
      return permute(RTYPE, rowBondNum);
    else if(typeID() == 2)         
      return permute(CTYPE, rowBondNum);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't permute in EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int):");
  }
  return *this;
}


UniTensor& UniTensor::permute(int* newLabels, int rowBondNum){
  try{
    if(typeID() == 1)         
      return permute(RTYPE, newLabels, rowBondNum);
    else if(typeID() == 2)         
      return permute(CTYPE, newLabels, rowBondNum);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't permute an EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int*, int):");
  }
	return *this;
}

UniTensor& UniTensor::permute(const std::vector<int>& newLabels, int rowBondNum){
  try{
    if(typeID() == 1)         
      return permute(RTYPE, newLabels, rowBondNum);
    else if(typeID() == 2)         
      return permute(CTYPE, newLabels, rowBondNum);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't permute an EMPTY block";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(std::vector<int>&, int):");
  }
  return *this;
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
    return UniTensor();
  }
}


std::vector<UniTensor> UniTensor::hosvd(size_t modeNum, size_t fixedNum)const{
  try{
    std::vector<std::map<Qnum, Matrix> > symLs;
    return _hosvd(modeNum, fixedNum, symLs, false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t = 0):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const{
  try{
    return _hosvd(modeNum, fixedNum, Ls, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(size_t modeNum, std::vector<std::map<Qnum, Matrix> >& Ls)const{
  try{
    return _hosvd(modeNum, 0, Ls, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const{
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
    std::vector<std::map<Qnum, Matrix> > symLs;
    const std::vector<UniTensor>& outs = _hosvd(modeNum, fixedNum, symLs, true);
    Ls.clear();
    Qnum q0(0);
    for(int i = 0; i < symLs.size(); i++)
      Ls.push_back(symLs[i][q0]);
    return outs;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t, std::vector<Matrix>&):");
    return std::vector<UniTensor>();
  }
}
std::vector<UniTensor> UniTensor::hosvd(size_t modeNum, std::vector<Matrix>& Ls)const{
  try{
    return hosvd(modeNum, 0, Ls);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, std::vector<Matrix>&):");
    return std::vector<UniTensor>();
  }
}

UniTensor& UniTensor::combineBond(const std::vector<int>&cmbLabels){
  try{
    if(typeID() == 1)         
      return combineBond(RTYPE, cmbLabels);
    else if(typeID() == 2)         
      return combineBond(CTYPE, cmbLabels);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't perform combineBond on and EMPTY tensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::combineBond(std::vector<int>&):");
  }
}

std::string UniTensor::printRawElem(bool print)const{
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
          if(typeID() == 1) 
            os<< std::setw(7) << std::fixed << std::setprecision(3) << at(idxs).real();
          else if(typeID() == 2) 
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
      if(typeID() == 1)
        os<<"\nScalar: " << elem[0]<<"\n\n";
      if(typeID() == 2)
        os<<"\nScalar: " << c_elem[0]<<"\n\n";
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

std::string UniTensor::profile(bool print){
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


void UniTensor::setRawElem(const Block& blk){
  try{
    if(blk.typeID() == 1) 
      setRawElem(RTYPE, blk);
    else if(blk.typeID() == 2) 
      setRawElem(CTYPE, blk);
    else if(blk.typeID() == 0){
      std::ostringstream err;
      err<<"Can't put an EMPTY　block in UniTensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(uni10::Block&):");
  }
}

void UniTensor::putBlock(const Block& mat){
  try{
    if(mat.typeID() == 1) 
      putBlock(RTYPE, mat);
    else if(mat.typeID() == 2) 
      putBlock(CTYPE, mat);
    else if(mat.typeID() == 0){
      std::ostringstream err;
      err<<"Can't put an EMPTY　block in UniTensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Block&):");
  }
}

void UniTensor::putBlock(const Qnum& qnum, const Block& mat){
  try{
    std::map<Qnum, Block>::iterator it;
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
    
    if(typeID() == 0){
      if(mat.typeID() == 1){
        r_flag = RTYPE;
        c_flag = CNULL;
        uelemAlloc(RTYPE);
        initBlocks(RTYPE);
      }
      else if(mat.typeID() == 2){
        r_flag = RNULL;
        c_flag = CTYPE;
        uelemAlloc(CTYPE);
        initBlocks(CTYPE);
      }
    }
    
    if(typeID() == 1){
      if(mat.m_elem != it->second.m_elem){
        if(mat.isDiag()){
          elemBzero(it->second.m_elem, it->second.Rnum * it->second.Cnum * sizeof(Real), ongpu);
          setDiag(it->second.m_elem, mat.getElem(RTYPE), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
        }
        else
          elemCopy(it->second.m_elem, mat.getElem(RTYPE), it->second.Rnum * it->second.Cnum * sizeof(Real), ongpu, mat.isOngpu());
      }
    }
    else if(typeID() == 2){
      if(mat.cm_elem != it->second.cm_elem){
        if(mat.isDiag()){
          elemBzero(it->second.cm_elem, it->second.Rnum * it->second.Cnum * sizeof(Complex), ongpu);
          setDiag(it->second.cm_elem, mat.getElem(CTYPE), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
        }
        else
          elemCopy(it->second.cm_elem, mat.getElem(CTYPE), it->second.Rnum * it->second.Cnum * sizeof(Complex), ongpu, mat.isOngpu());
      }
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Qnum&, uni10::Block&):");
  }
}

std::vector<_Swap> UniTensor::exSwap(const UniTensor& Tb) const{
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
void UniTensor::addGate(const std::vector<_Swap>& swaps){
  try{
    if(typeID() == 1) 
      addGate(RTYPE, swaps);
    else if(typeID() == 2) 
      addGate(CTYPE, swaps);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't addGate on an EMPTY block in UniTensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::addGate(std::vector<_Swap>&):");
  }
}

Complex UniTensor::trace()const{
  try{
    if(typeID() == 1) 
      return Complex(trace(RTYPE), 0);
    else if(typeID() == 2) 
      return trace(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't perform trace on an EMPTY block in UniTensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::trace():");
    return 0;
  }
}

UniTensor& UniTensor::partialTrace(int la, int lb){
  try{
    if(typeID() == 1) 
      return partialTrace(RTYPE, la, lb);
    else if(typeID() == 2) 
      return partialTrace(CTYPE, la, lb);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"Can't perform partialTrace on an EMPTY block in UniTensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::partialTrace(int, int):");
  }
}

Complex UniTensor::operator[](size_t idx)const{
  try{
    if(!(idx < m_elemNum)){
      std::ostringstream err;
      err<<"Index exceeds the number of elements("<<m_elemNum<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 1)
      return Complex(getElemAt(idx, elem, ongpu), 0);
    else if(typeID() == 2)
      return getElemAt(idx, c_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator[](size_t):");
    return 0;
  }
}

Matrix UniTensor::getRawElem()const{
  try{
    if(typeID() == 1) 
      return getRawElem(RTYPE);
    else if(typeID() == 2) 
      return getRawElem(CTYPE);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"This Tensor is EMPTY.";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getRawElem():");
    return Matrix();
  }
}

UniTensor& UniTensor::assign(const std::vector<Bond>& _bond){
  try{
    UniTensor T(_bond);
    *this = T;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::assign(std::vector<Bond>&):");
  }
  return *this;
}

bool UniTensor::isCelemEmpty(){
  return c_elem == NULL;
}
bool UniTensor::isElemEmpty(){
  return elem == NULL;
}

bool UniTensor::similar(const UniTensor& Tb)const{
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

bool UniTensor::elemCmp(const UniTensor& _UniT)const{
  try{
    if(typeID() == 0 || _UniT.typeID() == 0){
      std::ostringstream err;
      err<<"This Tensor is EMPTY.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UniTensor Ta(*this);
    UniTensor UniT(_UniT);
    if(Ta.typeID() != UniT.typeID()) 
      Ta.typeID() == 1 ? RtoC(Ta) : RtoC(UniT);

    double diff;
    if(m_elemNum == UniT.m_elemNum){
      if(Ta.typeID() == 1 && UniT.typeID() == 1){
        for(size_t i = 0; i < m_elemNum; i++){
          diff = std::abs(elem[i] - UniT.elem[i]);
          if(diff > 1E-12)
            return false;
        }
      }else{
        for(size_t i = 0; i < m_elemNum; i++){
          diff = std::abs(c_elem[i] - UniT.c_elem[i]);
          if(diff > 1E-12)
            return false;
        }
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

void UniTensor::clear(){
  status &= ~HAVEELEM;
}

Complex UniTensor::at(const std::vector<int>& idxs)const{
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

Complex UniTensor::at(const std::vector<size_t>& idxs)const{
  try{
    if(typeID() == 1) 
      return Complex(at(RTYPE, idxs), 0.0);
    else if(typeID() == 2) 
      return at(CTYPE, idxs);
    else if(typeID() == 0){
      std::ostringstream err;
      err<<"This Tensor is EMPTY.";
      throw std::runtime_error(exception_msg(err.str()));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::vector<size_t>&):");
    return 0;
  }
}
/*********************  neeee **********************/
/*********************  REAL **********************/

UniTensor::UniTensor(Real val): status(0){ //GPU
  try{
    initUniT(RTYPE);
    setElemAt(0, val, elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(double val):");
  }
}

UniTensor::UniTensor(rflag _tp, const std::vector<Bond>& _bonds, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT(RTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(rflag _tp, const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT(RTYPE);
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::vector<int>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(rflag _tp, const std::vector<Bond>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT(RTYPE);
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, int*, std::string& = \"\"):");
  }
}

void UniTensor::setRawElem(const std::vector<Real>& rawElem){
  try{
    setRawElem(&rawElem[0]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(std::vector<double>&):");
  }
}

void UniTensor::setRawElem(const Real* rawElem){
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"Setting elements to a tensor without bonds is not supported.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 2){
      std::ostringstream err;
      err<<"Can't set double* in a complex tensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 0){
      r_flag = RTYPE;
      c_flag = CNULL;
      uelemAlloc(RTYPE);
      initBlocks(RTYPE);
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
    Real* work = elem;
    if(ongpu){
      work = (Real*)malloc(m_elemNum * sizeof(Real));
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
      elemCopy(elem, work, m_elemNum * sizeof(Real), ongpu, false);
      free(work);
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(double*):");
  }
}

void UniTensor::setElem(const Real* _elem, bool _ongpu){
  try{
    if(typeID() == 2){
      std::ostringstream err;
      err<<"Can't set double* in a complex tensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
    else if(typeID() == 0){
      r_flag = RTYPE;
      c_flag = CNULL;
      uelemAlloc(RTYPE);
      initBlocks(RTYPE);
    }
    elemCopy(elem, _elem, m_elemNum * sizeof(Real), ongpu, _ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(double*, bool=false):");
  }
}

void UniTensor::setElem(const std::vector<Real>& _elem, bool _ongpu){
  try{
    setElem(&_elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(std::vector<double>&, bool=false):");
  }
}

std::map<Qnum, Matrix> UniTensor::getBlocks(rflag _tp)const{
  std::map<Qnum, Matrix> mats;
  try{
    for(std::map<Qnum, Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.m_elem, false, ongpu);
      mats.insert(std::pair<Qnum, Matrix>(it->first, mat));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlocks():");
  }
  return mats;
}

Matrix UniTensor::getBlock(rflag _tp, bool diag)const{
  try{
    Qnum q0(0);
    return getBlock(RTYPE, q0, diag);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(bool=false):");
    return Matrix();
  }
}

Matrix UniTensor::getBlock(rflag _tp, const Qnum& qnum, bool diag)const{
  try{
    std::map<Qnum, Block>::const_iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag)
      return it->second.getDiag();
    else{
      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.m_elem, false, ongpu);
      return mat;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(uni10::Qnum&):");
    return Matrix(0, 0);
  }
}

void UniTensor::set_zero(rflag _tp){
  try{
    elemBzero(elem, m_elemNum * sizeof(Real), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero():");
  }
}

void UniTensor::set_zero(rflag _tp, const Qnum& qnum){
  try{
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    elemBzero(block.m_elem, block.Rnum * block.Cnum * sizeof(Real), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero(std::Qnum&):");
  }
}

void UniTensor::identity(rflag _tp){
  try{
    std::map<Qnum, Block>::iterator it;
    for ( it = blocks.begin() ; it != blocks.end(); it++ )
      setIdentity(it->second.m_elem, it->second.Rnum, it->second.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity():");
  }
}

void UniTensor::identity(rflag _tp, const Qnum& qnum){
  try{
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    setIdentity(block.m_elem, block.Rnum, block.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity(std::Qnum&):");
  }
}

void UniTensor::randomize(rflag _tp){
  try{
    elemRand(elem, m_elemNum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::randomize():");
  }
}

void UniTensor::orthoRand(rflag _tp){
  try{
    std::map<Qnum, Block>::iterator it;
    for ( it = blocks.begin() ; it != blocks.end(); it++ )
      orthoRandomize(it->second.m_elem, it->second.Rnum, it->second.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand():");
  }
}

void UniTensor::orthoRand(rflag _tp, const Qnum& qnum){
  try{
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    orthoRandomize(block.m_elem, block.Rnum, block.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand(std::Qnum&):");
  }
}

UniTensor& UniTensor::transpose(rflag _tp){
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
    UniTensor UniTout(RTYPE, outBonds, name);
    UniTout.setLabel(outLabels);
    if(status & HAVEELEM){
      std::map<Qnum, Block>::iterator it_in;
      std::map<Qnum, Block>::iterator it_out;
      Real* elem_in;
      Real* elem_out;
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

UniTensor& UniTensor::permute(rflag _tp, int rowBondNum){
  try{
    std::vector<int> ori_labels = labels;
    this->permute(RTYPE, ori_labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int):");
  }
  return *this;
}

UniTensor& UniTensor::permute(rflag _tp, int* newLabels, int rowBondNum){
  try{
    std::vector<int> _labels(newLabels, newLabels + bonds.size());
    this->permute(RTYPE, _labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int*, int):");
  }
  return *this;
}

UniTensor& UniTensor::permute(rflag _tp, const std::vector<int>& newLabels, int rowBondNum){
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
      UniTensor UniTout(RTYPE, outBonds, name);
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
              Real* des_elem = UniTout.elem;
              Real* src_elem = elem;
              reshapeElem(src_elem, bondNum, m_elemNum, perInfo, des_elem);
              free(perInfo);
            }
            else{
              Real* des_elem = UniTout.elem;
              Real* src_elem = elem;
              size_t memsize = m_elemNum * sizeof(Real);
              if(ongpu){
                src_elem = (Real*)elemAllocForce(memsize, false);
                elemCopy(src_elem, elem, memsize, false, ongpu);
              }
              if(UniTout.ongpu)
                des_elem = (Real*)elemAllocForce(memsize, false);

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
            size_t memsize = m_elemNum * sizeof(Real);
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
          Real* Ein_ptr;
          Real* Eot_ptr;
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

Real UniTensor::at(rflag _tp, size_t idx)const{
  try{
    if(!(idx < m_elemNum)){
      std::ostringstream err;
      err<<"Index exceeds the number of elements("<<m_elemNum<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return getElemAt(idx, elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(rflag, size_t):");
    return 0;
  }
}

UniTensor contract(rflag _tp, UniTensor& Ta, UniTensor& Tb, bool fast){
  try{
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
      for(int a = 0; a < interLabel.size(); a++)
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

UniTensor otimes(rflag _tp, const UniTensor & Ta, const UniTensor& Tb){
  try{
    UniTensor T1 = Ta;
    UniTensor T2 = Tb;
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
    return contract(RTYPE, T1, T2, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function otimes(uni10::UniTensor&, uni10::UniTensor&):");
    return UniTensor();
  }
}

std::vector<UniTensor> UniTensor::hosvd(rflag _tp, size_t modeNum, size_t fixedNum)const{
  try{
    std::vector<std::map<Qnum, Matrix> > symLs;
    return _hosvd(RTYPE, modeNum, fixedNum, symLs, false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t = 0):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(rflag _tp, size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const{
  try{
    return _hosvd(RTYPE, modeNum, fixedNum, Ls, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(rflag _tp, size_t modeNum, std::vector<std::map<Qnum, Matrix> >& Ls)const{
  try{
    return _hosvd(RTYPE, modeNum, 0, Ls, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(rflag _tp, size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const{
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
    std::vector<std::map<Qnum, Matrix> > symLs;
    const std::vector<UniTensor>& outs = _hosvd(RTYPE , modeNum, fixedNum, symLs, true);
    Ls.clear();
    Qnum q0(0);
    for(int i = 0; i < symLs.size(); i++)
      Ls.push_back(symLs[i][q0]);
    return outs;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, size_t, std::vector<Matrix>&):");
    return std::vector<UniTensor>();
  }
}

std::vector<UniTensor> UniTensor::hosvd(rflag _tp, size_t modeNum, std::vector<Matrix>& Ls)const{
  try{
    return hosvd(RTYPE, modeNum, 0, Ls);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::hosvd(size_t, std::vector<Matrix>&):");
    return std::vector<UniTensor>();
  }
}

UniTensor& UniTensor::combineBond(rflag _tp, const std::vector<int>&cmbLabels){
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
    this->permute(RTYPE, rsp_labels, RBnum);
    UniTensor Tout(RTYPE, newBonds, reduced_labels);
    
    if(status & HAVEELEM)
      Tout.setElem(elem, ongpu);
    
    *this = Tout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::combineBond(std::vector<int>&):");
  }
  return *this;
}

void UniTensor::setRawElem(rflag _tp, const Block& blk){
  try{
    setRawElem(blk.getElem(RTYPE));
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(uni10::Block&):");
  }
}

void UniTensor::putBlock(rflag _tp, const Block& mat){
  try{
    Qnum q0(0);
    putBlock(RTYPE, q0, mat);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Block&):");
  }
}

void UniTensor::putBlock(rflag _tp, const Qnum& qnum, const Block& mat){
  try{
    std::map<Qnum, Block>::iterator it;
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
    if(typeID() == 2){
      std::ostringstream err;
      err<<"Can't put a REAL Block in a COMPLEX tensor."<<qnum<<std::endl;
      throw std::runtime_error(exception_msg(err.str()));
    }
    
    if(typeID() == 0){
      r_flag = RTYPE;
      c_flag = CNULL;
      uelemAlloc(RTYPE);
      initBlocks(RTYPE);
    }

    if(mat.m_elem != it->second.m_elem){
      if(mat.isDiag()){
        elemBzero(it->second.m_elem, it->second.Rnum * it->second.Cnum * sizeof(Real), ongpu);
        setDiag(it->second.m_elem, mat.getElem(RTYPE), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
      }
      else
        elemCopy(it->second.m_elem, mat.getElem(RTYPE), it->second.Rnum * it->second.Cnum * sizeof(Real), ongpu, mat.isOngpu());
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Qnum&, uni10::Block&):");
  }
}

void UniTensor::addGate(rflag _tp, const std::vector<_Swap>& swaps){
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
    Real* Eptr;
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

Real* UniTensor::getElem(rflag _tp){
  return elem;
}

Real UniTensor::trace(rflag _tp)const{
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot trace a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(status & HAVEBOND){
      size_t Rnum;
      Real trVal = 0;
      for(std::map<Qnum, Block>::const_iterator it = blocks.begin() ; it != blocks.end(); it++ ){
        if(!(it->second.Rnum == it->second.Cnum)){
          std::ostringstream err;
          err<<"Cannot trace a non-square block.";
          throw std::runtime_error(exception_msg(err.str()));
        }
        trVal += it->second.trace(RTYPE);
      }
      return trVal;
    }
    else{
        return getElemAt(0, elem, ongpu);
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::trace():");
    return 0;
  }
}

UniTensor& UniTensor::partialTrace(rflag _tp, int la, int lb){
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

    UniTensor Tt(RTYPE, newBonds, newLabels);
    rsp_labels[bondNum - 2] = labels[ia];
    rsp_labels[bondNum - 1] = labels[ib];
    ia = bondNum - 2;
    ib = bondNum - 1;
    this->permute(RTYPE, rsp_labels, Tt.RBondNum);
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
    Real* Et_ptr;
    std::vector<Real*> E_offs(tQdim);
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
      Real trVal;
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

Matrix UniTensor::getRawElem(rflag _tp)const{
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
      std::vector<Real> rawElem;
      while(1){
        rawElem.push_back(at(RTYPE, idxs));
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
      return Matrix(RTYPE, 1, 1, elem);
    else
      return Matrix();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getRawElem():");
    return Matrix();
  }
}

UniTensor& UniTensor::assign(rflag _tp, const std::vector<Bond>& _bond){
  try{
    UniTensor T(RTYPE, _bond);
    *this = T;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::assign(std::vector<Bond>&):");
  }
  return *this;
}

Real UniTensor::at(rflag _tp, const std::vector<int>& idxs)const{
  try{
    std::vector<size_t> _idxs(idxs.size());
    for(int i = 0; i < idxs.size(); i++)
      _idxs[i] = idxs[i];
    return at(RTYPE, _idxs);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::vector<int>&):");
    return 0;
  }
}

Real UniTensor::at(rflag _tp, const std::vector<size_t>& idxs)const{
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
      Block* blk = RQidx2Blk.find(Q_RQoff)->second;
      size_t B_cDim = blk->Cnum;
      size_t sB_cDim = CQidx2Dim.find(Q_CQoff)->second;
      size_t blkRoff = RQidx2Off.find(Q_RQoff)->second;
      size_t blkCoff = CQidx2Off.find(Q_CQoff)->second;
      Real* boff = blk->m_elem + (blkRoff * B_cDim) + blkCoff;
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
    propogate_exception(e, "In function UniTensor::at(std::vector<size_t>&):");
    return 0;
  }
}
/********************* reee **********************/
/*********************  COMPLEX **********************/

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
    initUniT(CTYPE);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(cflag _tp, const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT(CTYPE);
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::vector<int>&, std::string& = \"\"):");
  }
}

UniTensor::UniTensor(cflag _tp, const std::vector<Bond>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  try{
    initUniT(CTYPE);
    setLabel(_labels);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, int*, std::string& = \"\"):");
  }
}

void UniTensor::setRawElem(const std::vector<Complex>& rawElem){
  try{
    setRawElem(&rawElem[0]);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(std::vector<double>&):");
  }
}

void UniTensor::setRawElem(const Complex* rawElem){
  try{
    if((status & HAVEBOND) == 0){
      std::ostringstream err;
      err<<"Setting elements to a tensor without bonds is not supported.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 1){
      std::ostringstream err;
      err<<"Can't set complex<double>* in a real tensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 0){
      r_flag = RNULL; 
      c_flag = CTYPE;
      uelemAlloc(CTYPE);
      initBlocks(CTYPE);
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
    propogate_exception(e, "In function UniTensor::setRawElem(double*):");
  }
}

void UniTensor::setElem(const Complex* _elem, bool _ongpu){
  try{
    if(typeID() == 1){
      std::ostringstream err;
      err<<"Can't set complex<double>* in a real tensor";
      throw std::runtime_error(exception_msg(err.str()));
    }
    else if(typeID() == 0){
      r_flag = RNULL;
      c_flag = CTYPE;
      uelemAlloc(CTYPE);
      initBlocks(CTYPE);
    }
    elemCopy(c_elem, _elem, m_elemNum * sizeof(Complex), ongpu, _ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(double*, bool=false):");
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

std::map<Qnum, Matrix> UniTensor::getBlocks(cflag _tp)const{
  std::map<Qnum, Matrix> mats;
  try{
    for(std::map<Qnum, Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.cm_elem, false, ongpu);
      mats.insert(std::pair<Qnum, Matrix>(it->first, mat));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlocks():");
  }
  return mats;
}

Matrix UniTensor::getBlock(cflag _tp, bool diag)const{
  try{
    Qnum q0(0);
    return getBlock(CTYPE, q0, diag);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(bool=false):");
    return Matrix();
  }
}

Matrix UniTensor::getBlock(cflag _tp, const Qnum& qnum, bool diag)const{
  try{
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
    propogate_exception(e, "In function UniTensor::getBlock(uni10::Qnum&):");
    return Matrix(0, 0);
  }
}

void UniTensor::set_zero(cflag _tp){
  try{
    elemBzero(c_elem, m_elemNum * sizeof(Complex), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero():");
  }
}

void UniTensor::set_zero(cflag _tp, const Qnum& qnum){
  try{
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
    propogate_exception(e, "In function UniTensor::set_zero(std::Qnum&):");
  }
}

void UniTensor::identity(cflag _tp){
  try{
    std::map<Qnum, Block>::iterator it;
    for ( it = blocks.begin() ; it != blocks.end(); it++ )
      setIdentity(it->second.cm_elem, it->second.Rnum, it->second.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity():");
  }
}

void UniTensor::identity(cflag _tp, const Qnum& qnum){
  try{
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
    propogate_exception(e, "In function UniTensor::identity(std::Qnum&):");
  }
}

void UniTensor::randomize(cflag _tp){
  try{
    elemRand(c_elem, m_elemNum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::randomize():");
  }
}

void UniTensor::orthoRand(cflag _tp){
  try{
    std::map<Qnum, Block>::iterator it;
    for ( it = blocks.begin() ; it != blocks.end(); it++ )
      orthoRandomize(it->second.cm_elem, it->second.Rnum, it->second.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand():");
  }
}

void UniTensor::orthoRand(cflag _tp, const Qnum& qnum){
  try{
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
    propogate_exception(e, "In function UniTensor::orthoRand(std::Qnum&):");
  }
}

UniTensor& UniTensor::transpose(cflag _tp){
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
    propogate_exception(e, "In function UniTensor::transpose():");
  }
  return *this;
}

UniTensor& UniTensor::permute(cflag _tp, int rowBondNum){
  try{
    std::vector<int> ori_labels = labels;
    this->permute(CTYPE, ori_labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int):");
  }
  return *this;
}

UniTensor& UniTensor::permute(cflag _tp, int* newLabels, int rowBondNum){
  try{
    std::vector<int> _labels(newLabels, newLabels + bonds.size());
    this->permute(CTYPE, _labels, rowBondNum);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::permute(int*, int):");
  }
  return *this;
}

UniTensor& UniTensor::permute(cflag _tp, const std::vector<int>& newLabels, int rowBondNum){
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

Complex UniTensor::at(cflag _tp, size_t idx)const{
  try{
    if(!(idx < m_elemNum)){
      std::ostringstream err;
      err<<"Index exceeds the number of elements("<<m_elemNum<<").";
      throw std::runtime_error(exception_msg(err.str()));
    }
    return getElemAt(idx, c_elem, ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(cflag, size_t):");
    return 0;
  }
}

UniTensor contract(cflag _tp, UniTensor& Ta, UniTensor& Tb, bool fast){
  try{
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
      for(int a = 0; a < interLabel.size(); a++)
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

UniTensor otimes(cflag _tp, const UniTensor & Ta, const UniTensor& Tb){
  try{
    UniTensor T1 = Ta;
    UniTensor T2 = Tb;
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
    return contract(CTYPE, T1, T2, true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function otimes(uni10::UniTensor&, uni10::UniTensor&):");
    return UniTensor();
  }
}

UniTensor& UniTensor::combineBond(cflag _tp, const std::vector<int>&cmbLabels){
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
    this->permute(CTYPE, rsp_labels, RBnum);
    UniTensor Tout(CTYPE, newBonds, reduced_labels);
    
    if(status & HAVEELEM)
      Tout.setElem(c_elem, ongpu);
    
    *this = Tout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::combineBond(std::vector<int>&):");
  }
  return *this;
}

void UniTensor::setRawElem(cflag _tp, const Block& blk){
  try{
    setRawElem(blk.getElem(CTYPE));
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setRawElem(uni10::Block&):");
  }
}

void UniTensor::putBlock(cflag _tp, const Block& mat){
  try{
    Qnum q0(0);
    putBlock(CTYPE, q0, mat);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Block&):");
  }
}

void UniTensor::putBlock(cflag _tp, const Qnum& qnum, const Block& mat){
  try{
    std::map<Qnum, Block>::iterator it;
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
    if(typeID() == 1){
      std::ostringstream err;
      err<<"Can't put COMPLEX Block in a REAL tensor."<<qnum<<std::endl;
      throw std::runtime_error(exception_msg(err.str()));
    }
    
    if(typeID() == 0){
      r_flag = RNULL;
      c_flag = CTYPE;
      uelemAlloc(CTYPE);
      initBlocks(CTYPE);
    }

    if(mat.cm_elem != it->second.cm_elem){
      if(mat.isDiag()){
        elemBzero(it->second.cm_elem, it->second.Rnum * it->second.Cnum * sizeof(Complex), ongpu);
        setDiag(it->second.cm_elem, mat.getElem(CTYPE), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
      }
      else
        elemCopy(it->second.cm_elem, mat.getElem(CTYPE), it->second.Rnum * it->second.Cnum * sizeof(Complex), ongpu, mat.isOngpu());
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Qnum&, uni10::Block&):");
  }
}

void UniTensor::addGate(cflag _tp, const std::vector<_Swap>& swaps){
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

Complex* UniTensor::getElem(cflag _tp){
  return c_elem;
}

Complex UniTensor::trace(cflag _tp)const{
  try{
    if(!(status & HAVEELEM)){
      std::ostringstream err;
      err<<"Cannot trace a tensor before setting its elements.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(status & HAVEBOND){
      size_t Rnum;
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
    propogate_exception(e, "In function UniTensor::trace():");
    return 0;
  }
}

UniTensor& UniTensor::partialTrace(cflag _tp, int la, int lb){
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
    propogate_exception(e, "In function UniTensor::partialTrace(int, int):");
  }
  return *this;
}

Matrix UniTensor::getRawElem(cflag _tp)const{
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
    propogate_exception(e, "In function UniTensor::getRawElem():");
    return Matrix();
  }
}

UniTensor& UniTensor::assign(cflag _tp, const std::vector<Bond>& _bond){
  try{
    UniTensor T(CTYPE, _bond);
    *this = T;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::assign(std::vector<Bond>&):");
  }
  return *this;
}

Complex UniTensor::at(cflag _tp, const std::vector<int>& idxs)const{
  try{
    std::vector<size_t> _idxs(idxs.size());
    for(int i = 0; i < idxs.size(); i++)
      _idxs[i] = idxs[i];
    return at(CTYPE, _idxs);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::at(std::vector<int>&):");
    return 0;
  }
}

Complex UniTensor::at(cflag _tp, const std::vector<size_t>& idxs)const{
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
    propogate_exception(e, "In function UniTensor::at(std::vector<size_t>&):");
    return 0;
  }
}
/******************ceeee*******************/
/******************Friend funcs*******************/

/*****************************************************/

/*
Real* UniTensor::getElem(){
  return elem;
}
*/

/**************** Private Functions ***********************/
/*********************  NO TYPE **************************/	    

void UniTensor::initUniT(){
  r_flag = RNULL; 
  c_flag = CNULL; 
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
    blocks[q0] = Block(1, 1);
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
}

void UniTensor::initUniT(int _typeID){
  if(_typeID == 0)
    initUniT();
  else if(_typeID == 1)
    initUniT(RTYPE);
  else if(_typeID == 2)
    initUniT(CTYPE);
}

size_t UniTensor::grouping(){
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
    Block blk(it->second, it2->second); // blk(Rnum, Cnum);
    off += blk.Rnum * blk.Cnum;
    blocks[it->first] = blk;
    Block* blkptr = &(blocks[it->first]);
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
void UniTensor::uelemFree(){
  if(elem != NULL)
    elemFree(elem, sizeof(Real) * m_elemNum, ongpu);
  if(c_elem != NULL)
    elemFree(c_elem, sizeof(Complex) * m_elemNum, ongpu);
}

/*********************  RREAL **********************/
size_t UniTensor::grouping(rflag _tp){
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
    Block blk(RTYPE, it->second, it2->second); // blk(Rnum, Cnum);
    off += blk.Rnum * blk.Cnum;
    blocks[it->first] = blk;
    Block* blkptr = &(blocks[it->first]);
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

void UniTensor::initUniT(rflag _tp){ //GPU
  r_flag = RTYPE; 
  c_flag = CNULL; 
  if(bonds.size()){
    m_elemNum = grouping(RTYPE);
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
    blocks[q0] = Block(RTYPE, 1, 1);
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
  uelemAlloc(RTYPE);
  initBlocks(RTYPE);  
  uelemBzero(RTYPE);
}

void UniTensor::initBlocks(rflag _tp){
  size_t offset = 0;
  for(std::map<Qnum, Block>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
    it->second.r_flag = RTYPE; 
    it->second.c_flag = CNULL; 
    it->second.m_elem = &(elem[offset]);
    it->second.ongpu = ongpu;
    offset += it->second.Rnum * it->second.Cnum;
  }
}

void UniTensor::uelemAlloc(rflag _tp){
  elem = (Real*)elemAlloc(sizeof(Real) * m_elemNum, ongpu); 
}

void UniTensor::uelemBzero(rflag _tp){
  elemBzero(elem, sizeof(Real) * m_elemNum, ongpu);
}

/*********************  rre **********************/
/*********************  CCOMPLEX **********************/

void UniTensor::initUniT(cflag _tp){ //GPU
  r_flag= RNULL;
  c_flag= CTYPE;
  if(bonds.size()){
    m_elemNum = grouping(CTYPE);
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

  uelemAlloc(CTYPE);
  initBlocks(CTYPE);  
  uelemBzero(CTYPE);
}

size_t UniTensor::grouping(cflag _tp){
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
    Block blk(CTYPE, it->second, it2->second); // blk(Rnum, Cnum);
    off += blk.Rnum * blk.Cnum;
    blocks[it->first] = blk;
    Block* blkptr = &(blocks[it->first]);
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


void UniTensor::initBlocks(cflag _tp){
  size_t offset = 0;
  for(std::map<Qnum, Block>::iterator it = blocks.begin() ; it != blocks.end(); it++ ){
    it->second.r_flag = RNULL; 
    it->second.c_flag = CTYPE; 
    it->second.cm_elem = &(c_elem[offset]);
    it->second.ongpu = ongpu;
    offset += it->second.Rnum * it->second.Cnum;
  }
}

void UniTensor::uelemAlloc(cflag _tp){
  c_elem = (Complex*)elemAlloc(sizeof(Complex) * m_elemNum, ongpu);
}

void UniTensor::uelemBzero(cflag _tp){
  elemBzero(c_elem, sizeof(Complex) * m_elemNum, ongpu);
}

/******************  cce  *****************************/

std::vector<UniTensor> UniTensor::_hosvd(rflag _tp, size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const{
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
  UniTensor T(*this);
  for(int t = 0; t < T.labels.size(); t++)
    T.labels[t] = t;
  std::vector<int>ori_labels = T.labels;
  std::vector<int>rsp_labels = T.labels;
  if(returnL)
    Ls.assign(modeNum, std::map<Qnum, Matrix>());
  std::vector<UniTensor> Us;
  UniTensor S(T);
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
    T.permute(RTYPE, rsp_labels, combNum);
    std::vector<Bond> bonds(T.bonds.begin(), T.bonds.begin() + combNum);
    bonds.push_back(combine(bonds).dummy_change(BD_OUT));
    Us.push_back(UniTensor(RTYPE, bonds));
    for(std::map<Qnum, Block>::iterator it = T.blocks.begin(); it != T.blocks.end(); it++){
      std::vector<Matrix> svd = it->second.svd(RTYPE);
      Us[m].putBlock(it->first, svd[0]);
      if(returnL)
        Ls[m][it->first] = svd[1];
    }
    for(int c = 0; c < combNum; c++)
      Us[m].labels[c] = fixedNum + m * combNum + c;
    Us[m].labels[combNum] = -m - 1;
    out_labels[fixedNum + m] = -m -1;
    UniTensor UT = Us[m];
    S *= UT.transpose(RTYPE);
  }
  S.permute(RTYPE, out_labels, fixedNum);
  Us.push_back(S);
  return Us;
}

std::vector<UniTensor> UniTensor::_hosvd(cflag _tp, size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const{
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
  UniTensor T(*this);
  for(int t = 0; t < T.labels.size(); t++)
    T.labels[t] = t;
  std::vector<int>ori_labels = T.labels;
  std::vector<int>rsp_labels = T.labels;
  if(returnL)
    Ls.assign(modeNum, std::map<Qnum, Matrix>());
  std::vector<UniTensor> Us;
  UniTensor S(T);
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
    T.permute(CTYPE, rsp_labels, combNum);
    std::vector<Bond> bonds(T.bonds.begin(), T.bonds.begin() + combNum);
    bonds.push_back(combine(bonds).dummy_change(BD_OUT));
    Us.push_back(UniTensor(CTYPE, bonds));
    for(std::map<Qnum, Block>::iterator it = T.blocks.begin(); it != T.blocks.end(); it++){
      std::vector<Matrix> svd = it->second.svd(CTYPE);
      Us[m].putBlock(it->first, svd[0]);
      if(returnL)
        Ls[m][it->first] = svd[1];
    }
    for(int c = 0; c < combNum; c++)
      Us[m].labels[c] = fixedNum + m * combNum + c;
    Us[m].labels[combNum] = -m - 1;
    out_labels[fixedNum + m] = -m -1;
    UniTensor UT = Us[m];
    S *= UT.transpose(CTYPE);
  }
  S.permute(CTYPE, out_labels, fixedNum);
  Us.push_back(S);
  return Us;
}

std::vector<UniTensor> UniTensor::_hosvd(size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const{
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
  UniTensor T(*this);
  for(int t = 0; t < T.labels.size(); t++)
    T.labels[t] = t;
  std::vector<int>ori_labels = T.labels;
  std::vector<int>rsp_labels = T.labels;
  if(returnL)
    Ls.assign(modeNum, std::map<Qnum, Matrix>());
  std::vector<UniTensor> Us;
  UniTensor S(T);
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
    Us.push_back(UniTensor(bonds));
    for(std::map<Qnum, Block>::iterator it = T.blocks.begin(); it != T.blocks.end(); it++){
      std::vector<Matrix> svd = it->second.svd();
      Us[m].putBlock(it->first, svd[0]);
      if(returnL)
        Ls[m][it->first] = svd[1];
    }
    for(int c = 0; c < combNum; c++)
      Us[m].labels[c] = fixedNum + m * combNum + c;
    Us[m].labels[combNum] = -m - 1;
    out_labels[fixedNum + m] = -m -1;
    UniTensor UT = Us[m];
    S *= UT.transpose();
  }
  S.permute(out_labels, fixedNum);
  Us.push_back(S);
  return Us;
}

void RtoC(UniTensor& UniT){
  try{
    if(UniT.typeID() == 0){
      std::ostringstream err;
      err<<"This matrix is EMPTY. Nothing to do.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(UniT.typeID() == 1){
      UniT.r_flag = RNULL;
      UniT.c_flag = CTYPE;
      UniT.c_elem = (Complex*)elemAlloc( UniT.m_elemNum * sizeof(Complex), UniT.ongpu);
      elemCast(UniT.c_elem, UniT.elem, UniT.m_elemNum, UniT.ongpu, UniT.ongpu);
      size_t offset = 0;
      for(std::map<Qnum, Block>::iterator it = UniT.blocks.begin() ; it != UniT.blocks.end(); it++ ){
        RtoC(it->second);
        it->second.cm_elem = &(UniT.c_elem[offset]);
        it->second.ongpu = UniT.ongpu;
        offset += it->second.Rnum * it->second.Cnum;
      }
      UniT.elem = NULL;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::RtoC( )");
  }
}
void RtoC(Block& mat){
  try{
    if(mat.typeID() == 0){
      std::ostringstream err;
      err<<"This matrix is EMPTY. Nothing to do.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(mat.typeID() == 1){
      mat.r_flag = RNULL;
      mat.c_flag = CTYPE;
      mat.cm_elem = (Complex*)elemAlloc(mat.elemNum() * sizeof(Complex), mat.ongpu);
      elemCast(mat.cm_elem, mat.m_elem, mat.elemNum(), mat.ongpu, mat.ongpu);
      if(mat.m_elem != NULL)
        elemFree(mat.m_elem, mat.elemNum() * sizeof(Real), mat.ongpu);
      mat.m_elem = NULL;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::RtoC():");
  }
}
}; /* namespace uni10 */

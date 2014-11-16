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
#include <uni10/tensor-network/UniTensor.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>

namespace uni10{

int64_t UniTensor::ELEMNUM = 0;
int UniTensor::COUNTER = 0;
size_t UniTensor::MAXELEMNUM = 0;
size_t UniTensor::MAXELEMTEN = 0;

UniTensor::UniTensor(): status(0){
  try{
    initUniT();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor():");
  }
}

UniTensor::UniTensor(double val): status(0){
  try{
    initUniT();
    if(ongpu)
      setElemAt(0, val, elem, ongpu);
    else
      elem[0] = val;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor UniTensor::UniTensor(double val):");
  }
}

UniTensor::UniTensor(const UniTensor& UniT):
  status(UniT.status), bonds(UniT.bonds), blocks(UniT.blocks), labels(UniT.labels),
  RBondNum(UniT.RBondNum), RQdim(UniT.RQdim), CQdim(UniT.CQdim), m_elemNum(UniT.m_elemNum), elem(NULL),
  QidxEnc(UniT.QidxEnc), RQidx2Off(UniT.RQidx2Off), CQidx2Off(UniT.CQidx2Off), RQidx2Dim(UniT.RQidx2Dim), CQidx2Dim(UniT.CQidx2Dim){
    try{
      RQidx2Blk.clear();
      if(UniT.status & HAVEBOND){
        for(std::map<int, Block*>::const_iterator it = UniT.RQidx2Blk.begin(); it != UniT.RQidx2Blk.end(); it++)
          RQidx2Blk[it->first] = &(blocks[(it->second)->qnum]);
      }
      elem = (DOUBLE*)elemAlloc(sizeof(DOUBLE) * m_elemNum, ongpu);
      std::map<Qnum,Block>::iterator it;
      for ( it = blocks.begin() ; it != blocks.end(); it++ )
        it->second.elem = &(elem[it->second.offset]);
      ELEMNUM += m_elemNum;
      COUNTER++;
      if(ELEMNUM > MAXELEMNUM)
        MAXELEMNUM = ELEMNUM;
      if(m_elemNum > MAXELEMTEN)
        MAXELEMTEN = m_elemNum;
      elemCopy(elem, UniT.elem, sizeof(DOUBLE) * UniT.m_elemNum, ongpu, UniT.ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In copy constructor UniTensor::UniTensor(uni10::UniTensor&):");
    }
  }

UniTensor& UniTensor::operator=(const UniTensor& UniT){
  try{
    bonds = UniT.bonds;
    blocks = UniT.blocks;
    labels = UniT.labels;
    RBondNum = UniT.RBondNum;
    RQdim = UniT.RQdim;
    CQdim = UniT.CQdim;
    QidxEnc = UniT.QidxEnc;
    RQidx2Off = UniT.RQidx2Off;
    CQidx2Off = UniT.CQidx2Off;
    RQidx2Dim = UniT.RQidx2Dim;
    CQidx2Dim = UniT.CQidx2Dim;
    RQidx2Blk.clear();
    if(UniT.status & HAVEBOND){
      for(std::map<int, Block*>::const_iterator it = UniT.RQidx2Blk.begin(); it != UniT.RQidx2Blk.end(); it++)
        RQidx2Blk[it->first] = &(blocks[(it->second)->qnum]);
    }
    ELEMNUM -= m_elemNum;	//free original memory
    if(elem != NULL)
      elemFree(elem, sizeof(DOUBLE) * m_elemNum, ongpu);
    status = UniT.status;
    m_elemNum = UniT.m_elemNum;
    elem = (DOUBLE*)elemAlloc(sizeof(DOUBLE) * m_elemNum, ongpu);
    std::map<Qnum,Block>::iterator it;
    for ( it = blocks.begin(); it != blocks.end(); it++ )
      it->second.elem = &(elem[it->second.offset]);
    ELEMNUM += m_elemNum;
    if(ELEMNUM > MAXELEMNUM)
      MAXELEMNUM = ELEMNUM;
    if(m_elemNum > MAXELEMTEN)
      MAXELEMTEN = m_elemNum;
    elemCopy(elem, UniT.elem, sizeof(DOUBLE) * UniT.m_elemNum, ongpu, UniT.ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::operator=(uni10::UniTensor&):");
  }
  return *this;
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

UniTensor::UniTensor(const std::string& fname): status(0){	//load Tensor from file
  try{
    int namemax = 32;
    if(fname.size() > namemax)
      name = fname.substr(fname.size() - namemax);
    else
      name = fname;
    FILE* fp = fopen(fname.c_str(), "r");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in reading file '" << fname <<"'.";
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
      double *tmp_elem = elem;
      size_t memsize = m_elemNum * sizeof(double);
      if(ongpu)
        tmp_elem = (double*)malloc(memsize);
      size_t num_el;
      fread(&num_el, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
      fread(tmp_elem, m_elemNum, sizeof(DOUBLE), fp);
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

void UniTensor::initUniT(){
	if(bonds.size()){
		grouping();
		if(!(blocks.size() > 0)){ //No block in Tensor, Error!
      std::ostringstream err;
      err<<"There is no symmetry block with the given bonds:\n";
      for(int b = 0; b < bonds.size(); b++)
        err<<"    "<<bonds[b];
      throw std::runtime_error(exception_msg(err.str()));
    }
		Block blk = blocks.rbegin()->second;
		m_elemNum = blk.offset + (blk.Rnum * blk.Cnum);
		labels.assign(bonds.size(), 0);
		for(int b = 0; b < bonds.size(); b++)
			labels[b] = b;
		status |= HAVEBOND;
	}
	else{
		Qnum q0(0);
		Block blk;
		blk.Rnum = 1;
		blk.Cnum = 1;
		blk.qnum = q0;
		blk.offset = 0;
		blocks[q0] = blk;
		RBondNum = 0;
		RQdim = 0;
		CQdim = 0;
		m_elemNum = 1;
		status |= HAVEELEM;
	}
	elem = NULL;
	elem = (DOUBLE*)elemAlloc(sizeof(DOUBLE) * m_elemNum, ongpu);
	std::map<Qnum,Block>::iterator it;
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		it->second.elem = &(elem[it->second.offset]);
	elemBzero(elem, sizeof(DOUBLE) * m_elemNum, ongpu);
	ELEMNUM += m_elemNum;
	COUNTER++;
	if(ELEMNUM > MAXELEMNUM)
		MAXELEMNUM = ELEMNUM;
	if(m_elemNum > MAXELEMTEN)
		MAXELEMTEN = m_elemNum;
}

UniTensor::~UniTensor(){
	elemFree(elem, sizeof(DOUBLE) * m_elemNum, ongpu);
	ELEMNUM -= m_elemNum;
	COUNTER--;
}

size_t UniTensor::elemNum()const{return m_elemNum;}
size_t UniTensor::inBondNum()const{return RBondNum;}
size_t UniTensor::bondNum()const{return bonds.size();}

std::vector<Qnum> UniTensor::blockQnum()const{
	std::vector<Qnum> keys;
	for(std::map<Qnum,Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++)
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
    for(std::map<Qnum,Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
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

size_t UniTensor::blockNum()const{
	return blocks.size();
}

void UniTensor::profile(){
  std::cout<<"\n===== Tensor profile =====\n";
	std::cout<<"Existing Tensors: " << COUNTER << std::endl;
	std::cout<<"Allocated Elements: " << ELEMNUM << std::endl;
	std::cout<<"Max Allocated Elements: " << MAXELEMNUM << std::endl;
	std::cout<<"Max Allocated Elements for a Tensor: " << MAXELEMTEN << std::endl;
  std::cout<<"============================\n\n";
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

void UniTensor::setLabel(const std::vector<int>& newLabels){
  try{
    std::set<int> labelS(&(newLabels[0]), &(newLabels[newLabels.size()]));
    if(!(bonds.size() == labelS.size())){
      throw std::runtime_error(exception_msg("The size of input vector(labels) does not match for the number of bonds."));
    }
    labels = newLabels;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setLabel(std::vector<int>&):");
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
      size_t memsize = m_elemNum * sizeof(DOUBLE);
      double* tmp_elem = elem;
      if(ongpu){
        tmp_elem = (double*)malloc(memsize);
        elemCopy(tmp_elem, elem, memsize, false, ongpu);
      }
      fwrite(tmp_elem, m_elemNum, sizeof(DOUBLE), fp);
      if(ongpu)
        free(tmp_elem);
    }
    fclose(fp);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::save(std::string&):");
  }
}
/*------------------- SET ELEMENTS -----------------*/

void UniTensor::randomize(){
	elemRand(elem, m_elemNum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::orthoRand(const Qnum& qnum){
  try{
  std::map<Qnum, Block>::iterator it = blocks.find(qnum);
  if(it == blocks.end()){
    std::ostringstream err;
    err<<"There is no block with the given quantum number "<<qnum;
    throw std::runtime_error(exception_msg(err.str()));
  }
  Block& block = it->second;
	orthoRandomize(block.elem, block.Rnum, block.Cnum, ongpu);
	status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::orthoRand(std::Qnum&):");
  }
}

void UniTensor::orthoRand(){
	std::map<Qnum,Block>::iterator it;
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		orthoRandomize(it->second.elem, it->second.Rnum, it->second.Cnum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::identity(const Qnum& qnum){
  try{
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    setIdentity(block.elem, block.Rnum, block.Cnum, ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::identity(std::Qnum&):");
  }
}

void UniTensor::identity(){
	std::map<Qnum,Block>::iterator it;
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		setIdentity(it->second.elem, it->second.Rnum, it->second.Cnum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::set_zero(const Qnum& qnum){
  try{
    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    Block& block = it->second;
    elemBzero(block.elem, block.Rnum * block.Cnum * sizeof(DOUBLE), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero(std::Qnum&):");
  }
}

void UniTensor::set_zero(){
  try{
    elemBzero(elem, m_elemNum * sizeof(DOUBLE), ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::set_zero():");
  }
}

double* UniTensor::getElem(){
  return elem;
}

void UniTensor::setElem(const double* _elem, bool _ongpu){
  try{
    elemCopy(elem, _elem, m_elemNum * sizeof(double), ongpu, _ongpu);
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(double*, bool=false):");
  }
}

void UniTensor::setElem(const std::vector<double>& _elem, bool _ongpu){
  try{
    setElem(&_elem[0], _ongpu);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::setElem(std::vector<double>&, bool=false):");
  }
}

void UniTensor::setName(const std::string& _name){
	name = _name;
}

std::string UniTensor::getName(){
	return name;
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

//=============================ACCESS MEMORY EXPLICITLY=====================================

Matrix UniTensor::printRaw(bool flag)const{
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
		Bond rBond = combine(ins);
		Bond cBond = combine(outs);
		std::vector<Qnum> rowQ = rBond.Qlist();
		std::vector<Qnum> colQ = cBond.Qlist();
		size_t rowNum = rBond.dim();
		size_t colNum = cBond.dim();
		std::vector<size_t> idxs(bondNum, 0);

    if(flag){
      std::cout<< "     ";
      for(int q = 0; q < colQ.size(); q++)
        std::cout<< "   " << std::setw(2) << colQ[q].U1() << "," << colQ[q].prt();
      std::cout<< std::endl << std::setw(5) << "" << std::setw(colQ.size() * 7 + 2) <<std::setfill('-')<<"";
      std::cout<<std::setfill(' ');
    }
    int cnt = 0;
    int r = 0;
    int bend;
    std::vector<double> rawElem;
    while(1){
      if(flag){
        if(cnt % colNum == 0){
          std::cout<<"\n    |\n" << std::setw(2) << rowQ[r].U1() << "," << rowQ[r].prt() << "|";
          r++;
        }
        std::cout<< std::setw(7) << std::fixed << std::setprecision(3) << at(idxs);
      }
      rawElem.push_back(at(idxs));
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
    if(flag)
      std::cout <<"\n    |\n";
    return Matrix(rowNum, colNum, &rawElem[0]);
	}
	else if(status & HAVEELEM){
		std::cout<<"\nScalar: " << elem[0]<<"\n\n";
		return Matrix(1, 1, elem);
	}
	else{
		std::cout<<"NO ELEMENT IN THE TENSOR!!!\n";
		return Matrix(0, 0);
	}
}

void UniTensor::printRawElem()const{
  try{
	  printRaw(true);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::printRawElem():");
  }
}
Matrix UniTensor::getRawElem()const{
  try{
	  return printRaw(false);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getRawElem():");
    return Matrix();
  }
}

std::ostream& operator<< (std::ostream& os, const UniTensor& UniT){
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
    std::map<Qnum, Matrix> blocks = UniT.getBlocks();
    std::map<Qnum, Matrix>::const_iterator it;
    bool printElem = true;
    for ( it = blocks.begin() ; it != blocks.end(); it++ ){
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

Matrix UniTensor::getBlock(const Qnum& qnum, bool diag)const{
  try{
    std::map<Qnum, Block>::const_iterator it = blocks.find(qnum);
    if(it == blocks.end()){
      std::ostringstream err;
      err<<"There is no block with the given quantum number "<<qnum;
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag){
      Matrix mat(it->second.Rnum, it->second.Cnum, true, ongpu);
      getDiag(it->second.elem, mat.getElem(), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
      return mat;
    }
    else{
      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.elem, false, ongpu);
      return mat;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlock(uni10::Qnum&):");
    return Matrix(0, 0);
  }
}

std::map<Qnum, Matrix> UniTensor::getBlocks()const{
	std::map<Qnum, Matrix> mats;
  try{
    for(std::map<Qnum,Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.elem, false, ongpu);
      mats.insert(std::pair<Qnum, Matrix>(it->first, mat));
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::getBlocks():");
  }
	return mats;
}

void UniTensor::putBlock(const Qnum& qnum, const Matrix& mat){
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
    if(mat.isDiag()){
      elemBzero(it->second.elem, it->second.Rnum * it->second.Cnum * sizeof(DOUBLE), ongpu);
      setDiag(it->second.elem, mat.getElem(), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
    }
    else{
      elemCopy(it->second.elem, mat.getElem(), it->second.Rnum * it->second.Cnum * sizeof(DOUBLE), ongpu, mat.isOngpu());
    }
    status |= HAVEELEM;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::putBlock(uni10::Qnum&, uni10::Matrix&):");
  }
}

UniTensor& UniTensor::combineBond(const std::vector<int>&cmbLabels){
  try{
	if((status & HAVEBOND) == 0){
    std::ostringstream err;
    err<<"There is no bond in the tensor to be combined.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	if(!(cmbLabels.size() > 1)){
    std::ostringstream err;
    err<<"There should be at least two labels of bond in the input vector.";
    throw std::runtime_error(exception_msg(err.str()));
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
	UniTensor Tout(newBonds, reduced_labels);
	elemCopy(Tout.elem, elem, sizeof(DOUBLE) * m_elemNum, Tout.ongpu, ongpu);
	Tout.status |= HAVEELEM;
	*this = Tout;
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::combineBond(std::vector<int>&):");
  }
  return *this;
}

}; /* namespace uni10 */

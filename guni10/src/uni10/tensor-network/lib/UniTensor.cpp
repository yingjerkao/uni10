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
	initUniT();
}

UniTensor::UniTensor(double val): status(0){
	initUniT();
	if(ongpu)
		setElemAt(0, val, elem, ongpu);
	else
		elem[0] = val;
}

UniTensor::UniTensor(const UniTensor& UniT):
	status(UniT.status), bonds(UniT.bonds), blocks(UniT.blocks), labels(UniT.labels),
    RBondNum(UniT.RBondNum), RQdim(UniT.RQdim), CQdim(UniT.CQdim), m_elemNum(UniT.m_elemNum), elem(NULL),
	QidxEnc(UniT.QidxEnc), RQidx2Off(UniT.RQidx2Off), CQidx2Off(UniT.CQidx2Off), RQidx2Dim(UniT.RQidx2Dim), CQidx2Dim(UniT.CQidx2Dim){
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

UniTensor& UniTensor::operator=(const UniTensor& UniT){
	//std::cout<<"ASSING CONSTRUCTING " << this << std::endl;
	//name = UniT.name;
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
	return *this;
}

UniTensor& UniTensor::assign(const std::vector<Bond>& _bond){
	UniTensor T(_bond);
	return (*this = T);
}

UniTensor::UniTensor(const std::vector<Bond>& _bonds, const std::string& _name): name(_name), status(0), bonds(_bonds){
	//cout<<"CONSTRUCTING " << this << std::endl;
	initUniT();
}

UniTensor::UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
	initUniT();
	addLabel(_labels);
}
UniTensor::UniTensor(const std::vector<Bond>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
	initUniT();
	addLabel(_labels);
}

UniTensor::UniTensor(const std::string& fname): status(0){	//load Tensor from file
	int namemax = 32;
	if(fname.size() > namemax)
		name = fname.substr(fname.size() - namemax);
	else
		name = fname;
	FILE* fp = fopen(fname.c_str(), "r");
	assert(fp != NULL);
	int st;
	fread(&st, 1, sizeof(int), fp);
	int bondNum;
	fread(&bondNum, 1, sizeof(bondNum), fp);  //OUT: bondNum(4 bytes)
	size_t qnum_sz;
	fread(&qnum_sz, 1, sizeof(size_t), fp);	//OUT: sizeof(Qnum)
	assert(qnum_sz == sizeof(Qnum));
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
	assert(num_l == bonds.size());
	labels.assign(num_l, 0);
	fread(&(labels[0]), num_l, sizeof(int), fp);
	if(st & HAVEELEM){
		double *tmp_elem = elem;
		size_t memsize = m_elemNum * sizeof(double);
		if(ongpu)
			tmp_elem = (double*)malloc(memsize);
		size_t num_el;
		fread(&num_el, 1, sizeof(m_elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
		assert(num_el == m_elemNum);
		fread(tmp_elem, m_elemNum, sizeof(DOUBLE), fp);
		if(ongpu){
			elemCopy(elem, tmp_elem, memsize, ongpu, false);
			free(tmp_elem);
		}
		status |= HAVEELEM;
	}
	fclose(fp);
}

void UniTensor::initUniT(){
	if(bonds.size()){
		grouping();
		assert(blocks.size() > 0); //No block in Tensor, Error!
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
	//cout<<"DESTRUCTING " << this << std::endl;
	elemFree(elem, sizeof(DOUBLE) * m_elemNum, ongpu);
	ELEMNUM -= m_elemNum;
	COUNTER--;
}

size_t UniTensor::elemNum()const{return m_elemNum;}
int UniTensor::inBondNum()const{return RBondNum;}
size_t UniTensor::bondNum()const{return bonds.size();}

std::vector<Qnum> UniTensor::blockQnum()const{
	std::vector<Qnum> keys;
	for(std::map<Qnum,Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++)
		keys.push_back(it->first);
	return keys;
}
Qnum UniTensor::blockQnum(int idx)const{
	assert(idx < blocks.size());
	for(std::map<Qnum,Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
		if(idx == 0)
			return it->first;
		idx--;
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

void UniTensor::addLabel(int* newLabels){
	std::vector<int> labels(newLabels, newLabels + bonds.size());
	addLabel(labels);
}

void UniTensor::addLabel(const std::vector<int>& newLabels){
	std::set<int> labelS(&(newLabels[0]), &(newLabels[newLabels.size()]));
	assert(bonds.size() == labelS.size());
	labels = newLabels;
}

std::vector<int> UniTensor::label()const{
	return labels;
}

int UniTensor::label(int idx)const{
	assert(idx < labels.size());
	return labels[idx];
}

std::vector<Bond> UniTensor::bond()const{
	return bonds;
}

Bond UniTensor::bond(int idx)const{
	assert(idx < bonds.size());
	return bonds[idx];
}

UniTensor& UniTensor::permute(int rowBondNum){
	this->permute(labels, rowBondNum);
	return *this;
}
UniTensor& UniTensor::permute(int* newLabels, int rowBondNum){
	assert(status & HAVEBOND);
	std::vector<int> _labels(newLabels, newLabels + bonds.size());
	this->permute(_labels, rowBondNum);
	return *this;
}

void UniTensor::save(const std::string& fname){
	assert((status & HAVEBOND));   //If not INIT, NO NEED to write out to file
	FILE* fp = fopen(fname.c_str(), "w");
	assert(fp != NULL);
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
/*------------------- SET ELEMENTS -----------------*/

void UniTensor::randomize(){
	elemRand(elem, m_elemNum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::orthoRand(const Qnum& qnum){
	Block& block = blocks[qnum];
	orthoRandomize(block.elem, block.Rnum, block.Cnum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::orthoRand(){
	std::map<Qnum,Block>::iterator it;
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		orthoRandomize(it->second.elem, it->second.Rnum, it->second.Cnum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::identity(const Qnum& qnum){
	Block& block = blocks[qnum];
	setIdentity(block.elem, block.Rnum, block.Cnum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::identity(){
	std::map<Qnum,Block>::iterator it;
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		setIdentity(it->second.elem, it->second.Rnum, it->second.Cnum, ongpu);
	status |= HAVEELEM;
}

void UniTensor::set_zero(const Qnum& qnum){
	Block& block = blocks[qnum];
	elemBzero(block.elem, block.Rnum * block.Cnum * sizeof(DOUBLE), ongpu);
	status |= HAVEELEM;
}

void UniTensor::set_zero(){
	elemBzero(elem, m_elemNum * sizeof(DOUBLE), ongpu);
	status |= HAVEELEM;
}

void UniTensor::setName(const std::string& _name){
	name = _name;
}

std::string UniTensor::getName(){
	return name;
}

std::vector<_Swap> UniTensor::exSwap(const UniTensor& Tb) const{
	assert(status & Tb.status & HAVEBOND);
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
	std::vector<_Swap> swaps;
	_Swap sp;
	for(int i = 0; i < intersect.size(); i++)
		for(int j = 0; j < left.size(); j++){
			sp.b1 = intersect[i];
			sp.b2 = left[j];
			swaps.push_back(sp);
		}
	return swaps;
}

bool UniTensor::similar(const UniTensor& Tb)const{
	if(bonds.size() != Tb.bonds.size())
		return false;
	for(int b = 0; b < bonds.size(); b++){
		if(bonds[b] == Tb.bonds[b]);
		else return false;
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
		//std::cout<<"colNum: " << colNum<<std::endl;
		//std::cout<<"rowNum: " << rowNum<<std::endl;
		std::vector<int> idxs(bondNum, 0);

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
	printRaw(true);
}
Matrix UniTensor::rawElem()const{
	return printRaw(false);
}

std::ostream& operator<< (std::ostream& os, const UniTensor& UniT){
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
	return os;
}

Matrix UniTensor::getBlock(const Qnum& qnum, bool diag)const{
	std::map<Qnum, Block>::const_iterator it = blocks.find(qnum);
	if(it == blocks.end())
		return Matrix(0, 0);
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

std::map<Qnum, Matrix> UniTensor::getBlocks()const{
	std::map<Qnum, Matrix> mats;
	for(std::map<Qnum,Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
		Matrix mat(it->second.Rnum, it->second.Cnum, it->second.elem, false, ongpu);
		mats.insert(std::pair<Qnum, Matrix>(it->first, mat));
	}
	return mats;
}

void UniTensor::putBlock(const Qnum& qnum, const Matrix& mat){
	assert(blocks.find(qnum) != blocks.end());
	Block& blk = blocks[qnum];
	assert(mat.row() == blk.Rnum && mat.col() == blk.Cnum);
	if(mat.isDiag()){
		elemBzero(blk.elem, blk.Rnum * blk.Cnum * sizeof(DOUBLE), ongpu);
		setDiag(blk.elem, mat.getElem(), blk.Rnum, blk.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
	}
	else{
		elemCopy(blk.elem, mat.getElem(), blk.Rnum * blk.Cnum * sizeof(DOUBLE), ongpu, mat.isOngpu());
	}
	status |= HAVEELEM;
}

UniTensor& UniTensor::combineBond(const std::vector<int>&cmbLabels){
	assert(status & HAVEBOND);
	assert(cmbLabels.size() > 1);
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
	assert(mark == cmbLabels.size());

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
	return (*this = Tout);
}
}; /* namespace uni10 */

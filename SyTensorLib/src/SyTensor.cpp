#include "SyTensor.h"
int64_t SyTensor_t::ELEMNUM = 0;
int SyTensor_t::COUNTER = 0;
int64_t SyTensor_t::MAXELEMNUM = 0;
int64_t SyTensor_t::MAXELEMTEN = 0;

SyTensor_t::SyTensor_t(): status(0), elem(NULL), RBondNum(0), RQdim(0), CQdim(0), elemNum(0){
	COUNTER++;
}

SyTensor_t::SyTensor_t(const SyTensor_t& SyT):
	status(SyT.status), bonds(SyT.bonds), blocks(SyT.blocks), labels(SyT.labels),
    RBondNum(SyT.RBondNum), RQdim(SyT.RQdim), CQdim(SyT.CQdim), elemNum(SyT.elemNum), elem(NULL),
	QidxEnc(SyT.QidxEnc), RQidx2Off(SyT.RQidx2Off), CQidx2Off(SyT.CQidx2Off), RQidx2Dim(SyT.RQidx2Dim), CQidx2Dim(SyT.CQidx2Dim){
	//std::cout<<"COPY CONSTRUCTING " << this << std::endl;
	
	//status &= ~HAVELABEL;	//Labels are NOT copied to another tensor.

	RQidx2Blk.clear();
	for(std::map<int, Block_t*>::const_iterator it = SyT.RQidx2Blk.begin(); it != SyT.RQidx2Blk.end(); it++)
		RQidx2Blk[it->first] = &(blocks[(it->second)->qnum]);
	if(SyT.status & INIT){
		elem = (DOUBLE*)myMalloc(elem, sizeof(DOUBLE) * elemNum, status);
		std::map<Qnum_t,Block_t>::iterator it; 
		for ( it = blocks.begin() ; it != blocks.end(); it++ )
			it->second.elem = &(elem[it->second.offset]);
		ELEMNUM += elemNum;
		if(ELEMNUM > MAXELEMNUM)
			MAXELEMNUM = ELEMNUM;
		if(elemNum > MAXELEMTEN)
			MAXELEMTEN = elemNum;
		myMemcpy(elem, SyT.elem, sizeof(DOUBLE) * SyT.elemNum, status, SyT.status);
	}
	COUNTER++;
}

SyTensor_t& SyTensor_t::operator=(const SyTensor_t& SyT){
	//std::cout<<"ASSING CONSTRUCTING " << this << std::endl;
	//name = SyT.name;
	bonds = SyT.bonds;
	blocks = SyT.blocks;
	labels = SyT.labels;
	RBondNum = SyT.RBondNum;
	RQdim = SyT.RQdim;
	CQdim = SyT.CQdim;
	QidxEnc = SyT.QidxEnc;
	RQidx2Off = SyT.RQidx2Off;
	CQidx2Off = SyT.CQidx2Off;
	RQidx2Dim = SyT.RQidx2Dim;
	CQidx2Dim = SyT.CQidx2Dim;
	RQidx2Blk.clear();
	for(std::map<int, Block_t*>::const_iterator it = SyT.RQidx2Blk.begin(); it != SyT.RQidx2Blk.end(); it++)
		RQidx2Blk[it->first] = &(blocks[(it->second)->qnum]);
	if(SyT.status & INIT){
		ELEMNUM -= elemNum;	//free original memory
		if(elem != NULL)
			myFree(elem, sizeof(DOUBLE) * elemNum, status);
		status = SyT.status;
		elemNum = SyT.elemNum;
		elem = (DOUBLE*)myMalloc(elem, sizeof(DOUBLE) * elemNum, status);
		std::map<Qnum_t,Block_t>::iterator it; 
		for ( it = blocks.begin(); it != blocks.end(); it++ )
			it->second.elem = &(elem[it->second.offset]);
		ELEMNUM += elemNum;
		if(ELEMNUM > MAXELEMNUM)
			MAXELEMNUM = ELEMNUM;
		if(elemNum > MAXELEMTEN)
			MAXELEMTEN = elemNum;
		myMemcpy(elem, SyT.elem, sizeof(DOUBLE) * SyT.elemNum, status, SyT.status);
	}
	return *this;
}

SyTensor_t::SyTensor_t(std::vector<Bond_t>& _bonds, const std::string& _name): name(_name), status(0), bonds(_bonds){
	//cout<<"CONSTRUCTING " << this << std::endl;
	assert(_bonds.size() > 0); //No bond in Tensor, Error!
	initSyT();
	COUNTER++;
}

SyTensor_t::SyTensor_t(std::vector<Bond_t>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
	assert(_bonds.size() > 0); //No bond in Tensor, Error!
	initSyT();
	addLabel(_labels);
	COUNTER++;
}
SyTensor_t::SyTensor_t(std::vector<Bond_t>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
	assert(_bonds.size() > 0); //No bond in Tensor, Error!
	initSyT();
	addLabel(_labels);
	COUNTER++;
}

SyTensor_t::SyTensor_t(const std::string& fname): status(0){	//load Tensor from file
	name = fname;
	FILE* fp = fopen(fname.c_str(), "r");
	assert(fp != NULL);
	int st;
	fread(&st, 1, sizeof(int), fp);
	int bondNum;
	fread(&bondNum, 1, sizeof(bondNum), fp);  //OUT: bondNum(4 bytes)
	size_t qnum_sz;
	fread(&qnum_sz, 1, sizeof(size_t), fp);	//OUT: sizeof(Qnum_t)
	assert(qnum_sz == sizeof(Qnum_t));
	for(int b = 0; b < bondNum; b++){
		int num_q;
		bondType tp;
		fread(&tp, 1, sizeof(bondType), fp);	//OUT: Number of Qnums in the bond(4 bytes)
		fread(&num_q, 1, sizeof(int), fp);	//OUT: Number of Qnums in the bond(4 bytes)
		Qnum_t q0;
		std::vector<Qnum_t> qnums(num_q, q0);
		fread(&(qnums[0]), num_q, qnum_sz, fp);
		std::vector<int> qdegs(num_q, 0);
		fread(&(qdegs[0]), num_q, sizeof(int), fp);
		std::vector<Qnum_t> tot_qnums;
		for(int q = 0; q < num_q; q++)
			for(int d = 0; d < qdegs[q]; d++)
				tot_qnums.push_back(qnums[q]);
		Bond_t bd(tp, tot_qnums);
		bonds.push_back(bd);
	}
	initSyT();
	if(st & HAVELABEL){
		int num_l;
		fread(&num_l, 1, sizeof(int), fp);	//OUT: Number of Labels in the Tensor(4 bytes)
		assert(num_l == bonds.size());
		labels.assign(num_l, 0);
		fread(&(labels[0]), num_l, sizeof(int), fp);
		status |= HAVELABEL;
	}
	if(st & HAVEELEM){
		int num_el;
		fread(&num_el, 1, sizeof(elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
		assert(num_el == elemNum);
		fread(elem, elemNum, sizeof(DOUBLE), fp);
		status |= HAVEELEM;
	}
	fclose(fp);
}

SyTensor_t::~SyTensor_t(){
	//cout<<"DESTRUCTING " << this << std::endl;
	if(status & INIT){
		myFree(elem, sizeof(DOUBLE) * elemNum, status);
		ELEMNUM -= elemNum;
	}
	COUNTER--;
}

int64_t SyTensor_t::getElemNum()const{return elemNum;}
int SyTensor_t::getRBondNum()const{return RBondNum;}
int SyTensor_t::getBondNum()const{return bonds.size();}

std::vector<Qnum_t> SyTensor_t::qnums(){
	std::vector<Qnum_t> keys;
	for(std::map<Qnum_t,Block_t>::iterator it = blocks.begin(); it != blocks.end(); it++)
		keys.push_back(it->first);
	return keys;
}

void SyTensor_t::check(){
	std::cout<<"Existing Tensors: " << COUNTER << std::endl; 
	std::cout<<"Allocated Elem: " << ELEMNUM << std::endl;
	std::cout<<"Max Allocated Elem: " << MAXELEMNUM << std::endl;
	std::cout<<"Max Allocated Elem for a Tensor: " << MAXELEMTEN << std::endl;
}

void SyTensor_t::addLabel(int* newLabels){
	assert(status & INIT);
	std::vector<int> labels(newLabels, newLabels + bonds.size());
	addLabel(labels);
}

void SyTensor_t::addLabel(std::vector<int>& newLabels){
	assert(status & INIT);
	std::set<int> labelS(&(newLabels[0]), &(newLabels[newLabels.size()]));
	assert(bonds.size() == labelS.size());
	labels = newLabels;
	status |= HAVELABEL;
}

std::vector<int> SyTensor_t::getLabel()const{
	assert(status & HAVELABEL);
	return labels;
}
std::vector<_Swap> _recSwap(int* _ord, int n){	//Given the reshape order out to in. 
	int ordF[n];
	for(int i = 0; i < n; i++)
		ordF[i] = i;
	return _recSwap(_ord, n, ordF);
}
std::vector<_Swap> _recSwap(int* _ord, int n, int* ordF){	//Given the reshape order out to in. 
	int* ord = (int*)malloc(sizeof(int) * n);
	memcpy(ord, _ord, sizeof(int) * n);
	std::vector<_Swap> swaps;
	_Swap sg; 
	int tmp;
	for(int i = 0; i < n - 1; i++)
		for(int j = 0; j < n - i - 1; j++)
			if(ord[j] > ord[j + 1]){
				sg.b1 = ordF[ord[j + 1]]; 
				sg.b2 = ordF[ord[j]];
				tmp = ord[j];
				ord[j] = ord[j + 1];
				ord[j + 1] = tmp;
				swaps.push_back(sg);
			}
	free(ord);
	return swaps;
}

void SyTensor_t::reshape(int* newLabels, int rowBondNum){
	assert(status & INIT);
	std::vector<int> _labels(newLabels, newLabels + bonds.size());
	this->reshape(_labels, rowBondNum);

}
void SyTensor_t::initSyT(){
	grouping();
	assert(blocks.size() > 0); //No block in Tensor, Error!
	Block_t blk = blocks.rbegin()->second;
	elemNum = blk.offset + (blk.Rnum * blk.Cnum);
	elem = NULL;
	elem = (DOUBLE*)myMalloc(elem, sizeof(DOUBLE) * elemNum, status);
	//elem = (DOUBLE*)malloc(sizeof(DOUBLE) * elemNum);
	std::map<Qnum_t,Block_t>::iterator it; 
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		it->second.elem = &(elem[it->second.offset]);

	ELEMNUM += elemNum;
	if(ELEMNUM > MAXELEMNUM)
		MAXELEMNUM = ELEMNUM;
	if(elemNum > MAXELEMTEN)
		MAXELEMTEN = elemNum;
	membzero(elem, sizeof(DOUBLE) * elemNum, status);
	//memset(elem, 0, sizeof(DOUBLE) * elemNum);
	status |= INIT;
}


void SyTensor_t::save(const std::string& fname){
	assert((status & INIT));   //If not INIT, NO NEED to write out to file
	FILE* fp = fopen(fname.c_str(), "w");
	assert(fp != NULL);
	fwrite(&status, 1, sizeof(status), fp);	//OUT: status(4 bytes)
	int bondNum = bonds.size();
	fwrite(&bondNum, 1, sizeof(bondNum), fp);  //OUT: bondNum(4 bytes)
	size_t qnum_sz = sizeof(bonds[0].Qnums[0]);
	fwrite(&qnum_sz, 1, sizeof(size_t), fp);	//OUT: sizeof(Qnum_t)
	for(int b = 0; b < bondNum; b++){
		int num_q = bonds[b].Qnums.size();
		fwrite(&(bonds[b].type), 1, sizeof(bondType), fp);	//OUT: Number of Qnums in the bond(4 bytes)
		fwrite(&num_q, 1, sizeof(int), fp);	//OUT: Number of Qnums in the bond(4 bytes)
		fwrite(&(bonds[b].Qnums[0]), num_q, qnum_sz, fp);
		fwrite(&(bonds[b].Qdegs[0]), num_q, sizeof(int), fp);
	}
	if(status & HAVELABEL){
		int num_l = labels.size();
		fwrite(&num_l, 1, sizeof(int), fp);	//OUT: Number of Labels in the Tensor(4 bytes)
		fwrite(&(labels[0]), num_l, sizeof(int), fp);
	}
	if(status & HAVEELEM){
		fwrite(&elemNum, 1, sizeof(elemNum), fp);	//OUT: Number of elements in the Tensor(4 bytes)
		fwrite(elem, elemNum, sizeof(DOUBLE), fp);
	}
	fclose(fp);
}
/*------------------- SET ELEMENTS -----------------*/

void SyTensor_t::randomize(){
	assert((status & INIT));   //If not INIT, CANNOT add elements
	randomNums(elem, elemNum, status);
	status |= HAVEELEM;
}

void SyTensor_t::orthoRand(const Qnum_t& qnum){
	Block_t& block = blocks[qnum];
	orthoRandomize(block.elem, block.Rnum, block.Cnum);
}

void SyTensor_t::orthoRand(){
	std::map<Qnum_t,Block_t>::iterator it; 
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		orthoRandomize(it->second.elem, it->second.Rnum, it->second.Cnum);
	status |= HAVEELEM;
}

void SyTensor_t::eye(const Qnum_t& qnum){
	Block_t& block = blocks[qnum];
	myEye(block.elem, block.Rnum, block.Cnum, status);
}

void SyTensor_t::eye(){
	std::map<Qnum_t,Block_t>::iterator it; 
	for ( it = blocks.begin() ; it != blocks.end(); it++ )
		myEye(it->second.elem, it->second.Rnum, it->second.Cnum, status);
	status |= HAVEELEM;
}

void SyTensor_t::bzero(const Qnum_t& qnum){
	Block_t& block = blocks[qnum];
	membzero(block.elem, block.Rnum * block.Cnum * sizeof(DOUBLE), status);
	//memset(block.elem, 0, block.Rnum * block.Cnum * sizeof(DOUBLE));
}

void SyTensor_t::bzero(){
	membzero(elem, elemNum * sizeof(DOUBLE), status);
	//memset(elem, 0, elemNum * sizeof(DOUBLE));
}

void SyTensor_t::setName(const std::string& _name){
	name = _name;
}

std::string SyTensor_t::getName(){
	return name;
}

std::vector<_Swap> SyTensor_t::exSwap(const SyTensor_t& Tb) const{
	assert(status & HAVELABEL & Tb.status);
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

bool SyTensor_t::similar(const SyTensor_t& Tb)const{
	if(bonds.size() != Tb.bonds.size())
		return false;
	for(int b = 0; b < bonds.size(); b++){
		if(bonds[b] == Tb.bonds[b]);
		else return false;
	}
	return true;
}

//=============================ACCESS MEMORY EXPLICITLY=====================================

Matrix_t SyTensor_t::printRawElem(bool flag){
	if(status & HAVEELEM){
		int bondNum = bonds.size();
		int colNum = 1;
		int rowNum = 1;
		for(int b = bondNum - 1; b >= 0; b--){
			if(bonds[b].type == BD_COL)
				colNum *= bonds[b].dim;
			else
				rowNum *= bonds[b].dim;
		}
		std::vector<Qnum_t> rowQ;
		std::vector<Qnum_t> colQ;
		int Rnum = RBondNum;
		int Cnum = bondNum - RBondNum;
		std::vector<int> idxs(Rnum, 0);
		std::vector<int> qidxs(Rnum, 0);
		int bend;
		while(1){
			Qnum_t qnum;
			for(int b = 0; b < Rnum; b++)
				qnum = qnum * bonds[b].Qnums[qidxs[b]];
			rowQ.push_back(qnum);
			for(bend = Rnum - 1; bend >= 0; bend--){
				idxs[bend]++;
				if(idxs[bend] < bonds[bend].offsets[qidxs[bend]] + bonds[bend].Qdegs[qidxs[bend]])
					break;
				else{
					qidxs[bend]++;
					if(qidxs[bend] < bonds[bend].Qnums.size())
						break;
					else{
						qidxs[bend] = 0;
						idxs[bend] = 0;
					}
				}
			}
			if(bend < 0)
				break;
		}
		idxs.assign(Cnum, 0);
		qidxs.assign(Cnum, 0);
		while(1){
			Qnum_t qnum;
			for(int b = 0; b < Cnum; b++)
				qnum = qnum * bonds[Rnum + b].Qnums[qidxs[b]];
			colQ.push_back(qnum);
			for(bend = Cnum - 1; bend >= 0; bend--){
				idxs[bend]++;
				if(idxs[bend] < bonds[Rnum + bend].offsets[qidxs[bend]] + bonds[Rnum + bend].Qdegs[qidxs[bend]])
					break;
				else{
					qidxs[bend]++;
					if(qidxs[bend] < bonds[Rnum + bend].Qnums.size())
						break;
					else{
						qidxs[bend] = 0;
						idxs[bend] = 0;
					}
				}
			}
			if(bend < 0)
				break;
		}
		if(flag){
			std::cout<< "     ";
			for(int q = 0; q < colQ.size(); q++)
				std::cout<< "   " << std::setw(2) << colQ[q].getU1() << "," << colQ[q].getPrt();
			std::cout<< std::endl << std::setw(5) << "" << std::setw(colQ.size() * 7 + 2) <<std::setfill('-')<<"";
			std::cout<<std::setfill(' ');
		}
		idxs.assign(bondNum, 0);
		int cnt = 0;
		int r = 0;
		std::vector<double> rawElem;
		while(1){
			if(flag){
				if(cnt % colNum == 0){
					std::cout<<"\n    |\n" << std::setw(2) << rowQ[r].getU1() << "," << rowQ[r].getPrt() << "|";
					r++;
				}
				std::cout<< std::setw(7) << std::fixed << std::setprecision(3) << at(idxs);
			}
			rawElem.push_back(at(idxs));
			for(bend = bondNum - 1; bend >= 0; bend--){
				idxs[bend]++;
				if(idxs[bend] < bonds[bend].dim)
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
		return Matrix_t(rowNum, colNum, &rawElem[0]);
	}
	else{
		printf("NO ELEMENT IN THE TENSOR!!!\n");
		return Matrix_t(0, 0);
	}
}


std::ostream& operator<< (std::ostream& os, SyTensor_t& SyT){
	assert(SyT.status & SyT.INIT);
	int row = 0;
	int col = 0;
	for(int i = 0; i < SyT.bonds.size(); i++)
		if(SyT.bonds[i].type == BD_ROW)
			row++;
		else
			col++;
	int layer = std::max(row, col);
	int nmlen = SyT.name.length() + 2;
	int star = 12 + (14 - nmlen) / 2;
	os<<std::endl;
	for(int s = 0; s < star; s++)
		os << "*";
	if(SyT.name.length() > 0)
		os << " " << SyT.name << " ";
	for(int s = 0; s < star; s++)
		os<<"*";	
	os << "\n             ____________\n";
	os << "            |            |\n";
	int llab = 0;
	int rlab = 0;
	char buf[128];
	for(int l = 0; l < layer; l++){
		if(l < row && l < col){
			if(SyT.status & SyT.HAVELABEL){
				llab = SyT.labels[l];
				rlab = SyT.labels[row + l];
			}
			else{
				llab = l;
				rlab = row + l;
			}
			sprintf(buf, "    %5d___|%-4d    %4d|___%-5d\n", llab, SyT.bonds[l].dim, SyT.bonds[row + l].dim, rlab);
			os<<buf;
		}
		else if(l < row){
			if(SyT.status & SyT.HAVELABEL)
				llab = SyT.labels[l];
			else
				llab = l;
			sprintf(buf, "    %5d___|%-4d    %4s|\n", llab, SyT.bonds[l].dim, "");
			os<<buf;
		}
		else if(l < col){
			if(SyT.status & SyT.HAVELABEL)
				rlab = SyT.labels[row + l];
			else
				rlab = row + l;
			sprintf(buf, "    %5s   |%4s    %4d|___%-5d\n", "", "", SyT.bonds[row + l].dim, rlab);
			os << buf;
		}   
		os << "            |            |   \n";
	}   
	os << "            |____________|\n";

	os << "\n================BONDS===============\n";
	for(int b = 0; b < SyT.bonds.size(); b++){
		os << SyT.bonds[b];
	}   
	os<<"\n===============BLOCKS===============\n";
	std::map<Qnum_t,Block_t>::iterator it; 
	int Rnum, Cnum;
	bool printElem = true;
	for ( it = SyT.blocks.begin() ; it != SyT.blocks.end(); it++ ){

		Rnum = it->second.Rnum;
		Cnum = it->second.Cnum;
		os << "--- " << it->second.qnum << ": " << Rnum << " x " << Cnum << " = " << Rnum * Cnum << " ---\n\n";
		if((SyT.status & SyT.HAVEELEM) && printElem){ 
			double* elem = SyT.elem;
			for(int r = 0; r < Rnum; r++){
				for(int c = 0; c < Cnum; c++)
					os<< std::setw(7) << std::fixed << std::setprecision(3) << elem[it->second.offset + r * Cnum + c];
				os << "\n\n";
			}
		}
	}   
	os << "Total elemNum: "<<SyT.elemNum<<std::endl;
	os << "***************** END ****************\n\n";
	return os;
}

Matrix_t SyTensor_t::getBlock(Qnum_t qnum, bool diag){
	assert(blocks.find(qnum) != blocks.end());
	Block_t blk = blocks[qnum];
	if(diag){
		Matrix_t mat(blk.Rnum, blk.Cnum, true);
		int elemNum = blk.Rnum < blk.Cnum ? blk.Rnum : blk.Cnum;
		for(int i = 0; i < elemNum; i++)
			mat.elem[i] = blk.elem[i * blk.Cnum + i];
		return mat;
	}
	else{
		Matrix_t mat(blk.Rnum, blk.Cnum, blk.elem);
		return mat;
	}
}

std::map<Qnum_t, Matrix_t> SyTensor_t::getBlocks(){
	std::map<Qnum_t, Matrix_t> mats;
	for(std::map<Qnum_t,Block_t>::iterator it = blocks.begin(); it != blocks.end(); it++){
		Matrix_t mat(it->second.Rnum, it->second.Cnum, it->second.elem);
		mats.insert(std::pair<Qnum_t, Matrix_t>(it->first, mat));
	}
	return mats;
}

void SyTensor_t::putBlock(const Qnum_t& qnum, Matrix_t& mat){
	assert(blocks.find(qnum) != blocks.end());
	Block_t& blk = blocks[qnum];
	assert(mat.row() == blk.Rnum && mat.col() == blk.Cnum);
	if(mat.isDiag()){
		memset(blk.elem, 0, blk.Rnum * blk.Cnum * sizeof(DOUBLE));
		int elemNum = blk.Rnum < blk.Cnum ? blk.Rnum : blk.Cnum;
		for(int i = 0; i < elemNum; i++)
			blk.elem[i * blk.Cnum + i] = mat.elem[i];
	}
	else{
		memcpy(blk.elem, mat.elem, blk.Rnum * blk.Cnum * sizeof(DOUBLE));
	}
}

void SyTensor_t::combineIndex(const std::vector<int>&cmbLabels){
	assert(status & HAVELABEL);
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
	std::vector<Bond_t> newBonds;
	int RBnum = 0;
	for(int l = 0; l < labels.size(); l++){
		if(marked[l] && l == picked[0]){
			for(int ll = 0; ll < cmbLabels.size(); ll++){
				rsp_labels[enc] = cmbLabels[ll];
				enc++;
			}
			std::vector<Bond_t> tmpBonds;
			for(int p = 0; p < picked.size(); p++)
				tmpBonds.push_back(bonds[picked[p]]);
			if(bonds[picked[0]].type == BD_ROW)
				RBnum += picked.size();
			newBonds.push_back(combine(tmpBonds));
			reduced_labels[enc_r] = labels[l];
			enc_r++;
		}
		else if(marked[l] == 0){
			rsp_labels[enc] = labels[l];
			reduced_labels[enc_r] = labels[l];
			if(bonds[l].type == BD_ROW)
				RBnum++;
			newBonds.push_back(bonds[l]);
			enc_r++;
			enc++;
		}
	}
	/*
	for(int i = 0; i < reduced_labels.size(); i++)
		std::cout<<reduced_labels[i]<<std::endl;
	for(int i = 0; i < rsp_labels.size(); i++)
		std::cout<<rsp_labels[i]<<std::endl;
	*/
	/*
	std::map<Qnum_t,Block_t>::iterator it; 
	for ( it = blocks.begin() ; it != blocks.end(); it++ ){
		std::cout<<it->first<<" offset = "<<it->second.offset<<"\n";
	}
	for (it = Tout.blocks.begin() ; it != Tout.blocks.end(); it++ ){
		std::cout<<it->first<<" offset = "<<it->second.offset<<"\n";
	}
	*/
	this->reshape(rsp_labels, RBnum);
	SyTensor_t Tout(newBonds, reduced_labels);
	myMemcpy(Tout.elem, elem, sizeof(DOUBLE) * elemNum, Tout.status, Tout.status);
	Tout.status |= HAVEELEM;
	*this = Tout;
}

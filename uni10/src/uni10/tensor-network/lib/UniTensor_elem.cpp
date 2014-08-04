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
#include <uni10/tools/uni10_tools.h>
//using namespace uni10::datatype;
namespace uni10{
void UniTensor::grouping(){
	assert(bonds.size() > 0);
	blocks.clear();
	int row_bondNum = 0;
	int col_bondNum = 0;
	RQdim = 1;
	CQdim = 1;
	for(int i = 0; i < bonds.size(); i++){
		if(bonds[i].type() == BD_IN){
			RQdim *= bonds[i].Qnums.size();
			row_bondNum++;
		}
		else{
			CQdim *= bonds[i].Qnums.size();
			col_bondNum++;
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
	Block blk;
	int qidx;
	int off = 0;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		blk.Rnum = it->second;
		blk.Cnum = it2->second;
		blk.qnum = it2->first;
		blk.offset = off;
		off += blk.Rnum * blk.Cnum;
		blocks[it->first] = blk;
		Block* blkptr= &(blocks[it->first]);
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
}

void UniTensor::addGate(std::vector<_Swap> swaps){
	assert(status & HAVEBOND);
	assert(status & HAVEELEM);
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
	size_t E_off;
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
		E_off = RQidx2Blk[RQoff]->offset + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
		sB_rDim = RQidx2Dim[RQoff];
		sB_cDim = CQidx2Dim[CQoff];

		int sign01 = 0;
		for(int i = 0; i < swaps.size(); i++)
			sign01 ^= (bonds[swaps[i].b1].Qnums[Q_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Q_idxs[swaps[i].b2]].prtF());
		sign = sign01 ? -1 : 1;

		for(sB_r = 0; sB_r < sB_rDim; sB_r++)
			for(sB_c = 0; sB_c < sB_cDim; sB_c++)
				elem[E_off + (sB_r * B_cDim) + sB_c] *= sign;
	}
}

UniTensor& UniTensor::permute(const std::vector<int>& newLabels, int rowBondNum){
	assert(status & HAVEBOND);
	//assert(status & HAVELABEL);
	assert(labels.size() == newLabels.size());
	int bondNum = bonds.size();
	int rsp_outin[bondNum];	//rsp_outin[2] = 1 means the index "2" of UniTout is the index "1" of UniTin, opposite to the order in TensorLib
	int cnt = 0;
	for(int i = 0; i < bondNum; i++)
		for(int j = 0; j < bondNum; j++)
			if(labels[i] == newLabels[j]){
				rsp_outin[j] = i;
				cnt++;
			}

	assert(cnt == newLabels.size());
	bool inorder = true;
	for(int i = 1; i < bondNum; i++)
		if(((rsp_outin[i] + bondNum - i) % bondNum) != rsp_outin[0]){
			inorder = false;
			break;
		}
	if(inorder  && rsp_outin[0] == 0 && RBondNum == rowBondNum)	//do nothing
		return *this;
	else{
		std::vector<Bond> outBonds;
		for(int b = 0; b < bonds.size(); b++)
			outBonds.push_back(bonds[rsp_outin[b]]);
		for(int b = 0; b < bonds.size(); b++){
			if(b < rowBondNum)
				outBonds[b].change(BD_IN);
			else
				outBonds[b].change(BD_OUT);
		}
		UniTensor UniTout(outBonds, name);
		if(status & HAVEELEM){
			int sign = 1;
			//For Fermionic system
			std::vector<_Swap> swaps;
			if(Qnum::isFermionic()){
				int inLabelF[bondNum];
				int outLabelF[bondNum];
				int ordF[bondNum];
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

				int rspF_outin[bondNum];
				for(int i = 0; i < bondNum; i++)
					for(int j = 0; j < bondNum; j++)
						if(inLabelF[i] == outLabelF[j])
							rspF_outin[j] = i;
				swaps = recSwap(rspF_outin, bondNum, ordF);
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
			size_t Ein_off, Eot_off;
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
				Ein_off = RQidx2Blk[Qin_RQoff]->offset + (RQidx2Off[Qin_RQoff] * Bin_cDim) + CQidx2Off[Qin_CQoff];
				Eot_off = UniTout.RQidx2Blk[Qot_RQoff]->offset + (UniTout.RQidx2Off[Qot_RQoff] * Bot_cDim) + UniTout.CQidx2Off[Qot_CQoff];
				sBin_rDim = RQidx2Dim[Qin_RQoff];
				sBin_cDim = CQidx2Dim[Qin_CQoff];
				sBot_cDim = UniTout.CQidx2Dim[Qot_CQoff];
				int cnt_ot = 0;
				sBin_idxs.assign(bondNum, 0);
				if(Qnum::isFermionic()){
					int sign01 = 0;
					for(int i = 0; i < swaps.size(); i++)
						sign01 ^= (bonds[swaps[i].b1].Qnums[Qin_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Qin_idxs[swaps[i].b2]].prtF());
					sign = sign01 ? -1 : 1;
				}
				for(sBin_r = 0; sBin_r < sBin_rDim; sBin_r++)
					for(sBin_c = 0; sBin_c < sBin_cDim; sBin_c++){
						sBot_r = cnt_ot / sBot_cDim;
						sBot_c = cnt_ot % sBot_cDim;
						UniTout.elem[Eot_off + (sBot_r * Bot_cDim) + sBot_c] = sign * elem[Ein_off + (sBin_r * Bin_cDim) + sBin_c];
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
			UniTout.status |= HAVEELEM;
		}
		*this = UniTout;
		this->addLabel(newLabels);
		return *this;
	}
}

void UniTensor::addRawElem(std::vector<DOUBLE> rawElem){
  addRawElem(&rawElem[0]);
}
void UniTensor::addRawElem(DOUBLE* rawElem){
	assert((status & HAVEBOND));   //If not INIT, CANNOT add elements
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
		E_off = RQidx2Blk[RQoff]->offset + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
		sB_rDim = RQidx2Dim[RQoff];
		sB_cDim = CQidx2Dim[CQoff];

		sB_idxs.assign(bondNum, 0);
		for(sB_r = 0; sB_r < sB_rDim; sB_r++)
			for(sB_c = 0; sB_c < sB_cDim; sB_c++){
				elem[E_off + (sB_r * B_cDim) + sB_c] = rawElem[R_off];
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
	status |= HAVEELEM;
}

<<<<<<< HEAD
=======
/*
>>>>>>> eb8ff2a4ac324948ab953f0306a14dc5c5d03a77
void UniTensor::elemSet(const Qnum& qnum, double* _elem){
	Block& block = blocks[qnum];
	memcpy(block.elem, _elem, block.Rnum * block.Cnum * sizeof(DOUBLE));
	status |= HAVEELEM;
}

void UniTensor::elemSet(double* _elem){
	memcpy(elem, _elem, elemNum() * sizeof(DOUBLE));
	status |= HAVEELEM;
<<<<<<< HEAD
}
=======
}*/
>>>>>>> eb8ff2a4ac324948ab953f0306a14dc5c5d03a77

double UniTensor::at(std::vector<int> idxs)const{
	assert(status & HAVEBOND);
	int bondNum = bonds.size();
	std::vector<int> Qidxs(bondNum, 0);
	for(int b = 0; b < bondNum; b++){
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
		int boff = blk->offset + (blkRoff * B_cDim) + blkCoff;
		int cnt = 0;
		std::vector<int> D_acc(bondNum, 1);
		for(int b = bondNum	- 1; b > 0; b--)
			D_acc[b - 1] = D_acc[b] * bonds[b].Qdegs[Qidxs[b]];
		for(int b = 0; b < bondNum; b++)
			cnt += (idxs[b] - bonds[b].offsets[Qidxs[b]]) * D_acc[b];
		return elem[boff + (cnt / sB_cDim) * B_cDim + cnt % sB_cDim];
	}
	else{
		return 0.0;;
	}
}
double UniTensor::operator[](size_t idx){
	assert(idx < m_elemNum);
	return elem[idx];
}

UniTensor& UniTensor::transpose(){
	assert(status & HAVEBOND);
	int bondNum = bonds.size();
	int rsp_outin[bondNum];	//rsp_outin[2] = 1 means the index "2" of UniTout is the index "1" of UniTin, opposite to the order in TensorLib
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
	//if(status & HAVELABEL)
	for(int b = 0; b < bonds.size(); b++){
		outBonds.push_back(bonds[rsp_outin[b]]);
		outLabels[b] = labels[rsp_outin[b]];
	}
	//else
	//	for(int b = 0; b < bonds.size(); b++)
	//		outBonds.push_back(bonds[rsp_outin[b]]);

	for(int b = 0; b < bondNum; b++){
		if(b < cbondNum)
			outBonds[b].m_type = BD_IN;
		else
			outBonds[b].m_type = BD_OUT;
	}
	UniTensor UniTout(outBonds, name);
	//if(status & HAVELABEL)
	UniTout.addLabel(outLabels);
	if(status & HAVEELEM){
		std::map<Qnum,Block>::iterator it_in;
		std::map<Qnum,Block>::iterator it_out;
		double* elem_in;
		double* elem_out;
		size_t Rnum, Cnum;
		for ( it_in = blocks.begin() ; it_in != blocks.end(); it_in++ ){
			it_out = UniTout.blocks.find((it_in->first));
			Rnum = it_in->second.Rnum;
			Cnum = it_in->second.Cnum;
			elem_in = it_in->second.elem;
			elem_out = it_out->second.elem;
			for(size_t i = 0; i < Rnum; i++)
				for(size_t j = 0; j < Cnum; j++)
					elem_out[j * Rnum + i] = elem_in[i * Cnum + j];
		}
		UniTout.status |= HAVEELEM;
	}
	return *this = UniTout;
}

bool UniTensor::elemCmp(const UniTensor& UniT)const{
	double diff;
	if(m_elemNum == UniT.m_elemNum){
		for(size_t i = 0; i < m_elemNum; i++){
			diff = fabs(elem[i] - UniT.elem[i]);
			if(diff > 1E-6)
				return false;
		}
	}
	else
		return false;
	return true;
}
double UniTensor::trace()const{
	assert(status & HAVEELEM);
	if(status & HAVEBOND){
		size_t Rnum;
		DOUBLE trVal = 0;
		for(std::map<Qnum, Block>::const_iterator it = blocks.begin() ; it != blocks.end(); it++ ){
			assert(it->second.Rnum == it->second.Cnum);
			Rnum = it->second.Rnum;
			for(size_t r = 0; r < Rnum; r++)
				trVal += it->second.elem[r * Rnum + r];
		}
		return trVal;
	}
	else
		return elem[0];
}

UniTensor& UniTensor::partialTrace(int la, int lb){
	assert(status & HAVEELEM);
	assert(status & HAVEBOND);
	assert(bonds.size() > 2 && la != lb);
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
	assert(enc == newLabels.size());

	UniTensor Tt(newBonds, newLabels);
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
	assert(tQdim == bonds[ib].Qnums.size());
	Qnum q0(0, PRT_EVEN);
	for(int q = 0; q < tQdim; q++){
		assert(bonds[ia].Qnums[q] * bonds[ib].Qnums[q] == q0);
		assert(bonds[ia].Qdegs[q] == bonds[ib].Qdegs[q]);
	}
	/*END*/
	int tBnum = Tt.bonds.size();
	//std::vector<int> Q_idxs(bondNum, 0);
	std::vector<int> Qt_Bdims(tBnum, 0);
	for(int b = 0; b < tBnum; b++)
		Qt_Bdims[b] = Tt.bonds[b].Qnums.size();

	int Qt_off;
	int Q_off;
	//int tmp;
	int Qt_RQoff, Qt_CQoff;
	int Q_RQoff, Q_CQoff;
	size_t sBt_rDim, sBt_cDim;	//sub-block of a Qidx of Tt
	size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
	size_t Bt_cDim;
	//int B_cDim;
<<<<<<< HEAD
	int64_t Et_off;
	std::vector<int64_t> E_offs(tQdim);
	std::vector<int> B_cDims(tQdim);
=======
	size_t Et_off;
	std::vector<size_t> E_offs(tQdim);
	std::vector<size_t> B_cDims(tQdim);
>>>>>>> eb8ff2a4ac324948ab953f0306a14dc5c5d03a77
	int tQdim2 = tQdim * tQdim;
	int Qenc = Q_acc[ia] + Q_acc[ib];
	for(std::map<int, size_t>::iterator it = Tt.QidxEnc.begin(); it != Tt.QidxEnc.end(); it++){
		Qt_off = it->first;

		Qt_RQoff = Qt_off / Tt.CQdim;
		Qt_CQoff = Qt_off % Tt.CQdim;
		Bt_cDim = Tt.RQidx2Blk[Qt_RQoff]->Cnum;
		Et_off = Tt.RQidx2Blk[Qt_RQoff]->offset + (Tt.RQidx2Off[Qt_RQoff] * Bt_cDim) + Tt.CQidx2Off[Qt_CQoff];
		sBt_rDim = Tt.RQidx2Dim[Qt_RQoff];
		sBt_cDim = Tt.CQidx2Dim[Qt_CQoff];

		for(int q = 0; q < tQdim; q++){
			Q_off = Qt_off * tQdim2 + q * Qenc;
			Q_RQoff = Q_off / CQdim;
			Q_CQoff = Q_off % CQdim;
			B_cDims[q] = RQidx2Blk[Q_RQoff]->Cnum;
			E_offs[q] = RQidx2Blk[Q_RQoff]->offset + (RQidx2Off[Q_RQoff] * B_cDims[q]) + CQidx2Off[Q_CQoff];
		}
		int tQdeg, sB_c_off;
		DOUBLE trVal;
		for(size_t sB_r = 0; sB_r < sBt_rDim; sB_r++)
			for(size_t sB_c = 0; sB_c < sBt_cDim; sB_c++){
				trVal = 0;
				for(int q = 0; q < tQdim; q++){
					tQdeg = bonds[ia].Qdegs[q];
					sB_c_off = sB_c * (tQdeg * tQdeg);
					for(int t = 0; t < tQdeg; t++){
						trVal += elem[E_offs[q] + (sB_r * B_cDims[q]) + sB_c_off + t * (tQdeg + 1)];
					}
				}
				Tt.elem[Et_off + sB_r * Bt_cDim + sB_c] = trVal;
			}
		Tt.status |= HAVEELEM;
	}
	return *this = Tt;
}
}; /* namespace uni10 */

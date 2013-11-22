#include "SyTensor.h"

void SyTensor_t::grouping(){
	blocks.clear();
	int row_bondNum = 0;
	int col_bondNum = 0;
	RQdim = 1;
	CQdim = 1;
	for(int i = 0; i < bonds.size(); i++){
		if(bonds[i].type == BD_ROW){
			RQdim *= bonds[i].Qnums.size();
			row_bondNum++;
		}
		else{
			CQdim *= bonds[i].Qnums.size();
			col_bondNum++;
		}
	}
	RBondNum = row_bondNum;
	map<Qnum_t,int> row_QnumMdim;
	vector<int> row_offs(row_bondNum, 0);
	map<Qnum_t,vector<int> > row_Qnum2Qidx;
	Qnum_t qnum;
	int dim;
	int boff = 0;
	RQidx2Off = (int*)calloc(RQdim, sizeof(int));
	RQidx2Dim = (int*)calloc(RQdim, sizeof(int));
	CQidx2Off = (int*)calloc(CQdim, sizeof(int));
	CQidx2Dim = (int*)calloc(CQdim, sizeof(int));
	Qidx = (bool*)calloc(RQdim * CQdim, sizeof(bool));
	RQidx2Blk = (Block_t**)calloc(RQdim, sizeof(Block_t*));
	if(row_bondNum){
		while(1){
			qnum.set();
			dim = 1;
			for(int b = 0; b < row_bondNum; b++){
				qnum = qnum * bonds[b].Qnums[row_offs[b]];
				dim *= bonds[b].Qdegs[row_offs[b]];
			}
			if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
				RQidx2Off[boff] = row_QnumMdim[qnum];
				RQidx2Dim[boff] = dim;
				row_QnumMdim[qnum] += dim;
			}
			else{
				RQidx2Off[boff] = 0;
				RQidx2Dim[boff] = dim;
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
		qnum.set();
		row_QnumMdim[qnum] = 1;
		row_Qnum2Qidx[qnum].push_back(0);
		RQidx2Dim[0] = 1;
	}
	map<Qnum_t,int>::iterator itt;
	map<Qnum_t,int> col_QnumMdim;
	vector<int> col_offs(col_bondNum, 0);
	map<Qnum_t,vector<int> > col_Qnum2Qidx;
	boff = 0;
	if(col_bondNum){
		while(1){
			qnum.set();
			dim = 1;
			for(int b = 0; b < col_bondNum; b++){
				qnum = qnum * bonds[b + row_bondNum].Qnums[col_offs[b]];
				dim *= bonds[b + row_bondNum].Qdegs[col_offs[b]];
			}
			if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
				if(col_QnumMdim.find(qnum) != col_QnumMdim.end()){
					CQidx2Off[boff] = col_QnumMdim[qnum];
					CQidx2Dim[boff] = dim;
					col_QnumMdim[qnum] += dim;
				}
				else{
					CQidx2Off[boff] = 0;
					CQidx2Dim[boff] = dim;
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
		qnum.set();
		if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
			col_QnumMdim[qnum] = 1;
			col_Qnum2Qidx[qnum].push_back(0);
			CQidx2Dim[0] = 1;
		}
	}

	map<Qnum_t,int>::iterator it;
	map<Qnum_t,int>::iterator it2;
	Block_t blk;
	int off = 0;
	QidxNum = 0;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		blk.Rnum = it->second;
		blk.Cnum = it2->second;
		blk.qnum = it2->first;
		blk.offset = off;
		off += blk.Rnum * blk.Cnum;
		blocks[it->first] = blk;
		Block_t* blkptr= &(blocks[it->first]);
		vector<int>& tmpRQidx = row_Qnum2Qidx[it->first];
		vector<int>& tmpCQidx = col_Qnum2Qidx[it->first];	
		for(int i = 0; i < tmpRQidx.size(); i++){
			RQidx2Blk[tmpRQidx[i]] = blkptr;
			for(int j = 0; j < tmpCQidx.size(); j++){
				Qidx[tmpRQidx[i] * CQdim + tmpCQidx[j]] = true;
				QidxNum++;
			}
		}
	}
	Eidx2Qidx = (int*)malloc(QidxNum * sizeof(int));
	Eidx = (int*)malloc(QidxNum * sizeof(int));
	int qidx;
	int enc = 0;
	int elemEnc = 0;
	for(int r = 0; r < RQdim; r++)
		for(int c = 0; c < CQdim; c++){
			qidx = r * CQdim + c;
			if(Qidx[qidx]){
				Eidx2Qidx[enc] = qidx;
				Eidx[enc] = elemEnc;
				elemEnc += RQidx2Dim[r] * CQidx2Dim[c];
				enc++;
			}
		}
}

void SyTensor_t::addGate(vector<_Swap> swaps){
	assert(status & HAVEELEM);
	assert(!(status & ONGPU));
	int sign = 1;
	int bondNum = bonds.size();
	vector<int> Q_idxs(bondNum, 0); 
	vector<int> Q_Bdims(bondNum, 0);
	for(int b = 0; b < bondNum; b++)
		Q_Bdims[b] = bonds[b].Qnums.size();
	int Q_off;
	int tmp;
	int RQoff, CQoff;
	int sB_r, sB_c;	//sub-block of a Qidx
	int sB_rDim, sB_cDim;	//sub-block of a Qidx
	int B_cDim;
	int64_t E_off;
	for(int enc = 0; enc < QidxNum; enc++){
		Q_off = Eidx2Qidx[enc];
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
			sign01 ^= (bonds[swaps[i].b1].Qnums[Q_idxs[swaps[i].b1]].getPrtF() & bonds[swaps[i].b2].Qnums[Q_idxs[swaps[i].b2]].getPrtF());
		sign = sign01 ? -1 : 1;
		
		for(sB_r = 0; sB_r < sB_rDim; sB_r++)
			for(sB_c = 0; sB_c < sB_cDim; sB_c++)
				elem[E_off + (sB_r * B_cDim) + sB_c] *= sign;
	}
}

void SyTensor_t::reshape(vector<int>& newLabels, int rowBondNum){
	assert(status & INIT);
	assert(status & HAVELABEL);
	assert(labels.size() == newLabels.size());
	int bondNum = bonds.size();
	int rsp_outin[bondNum];	//rsp_outin[2] = 1 means the index "2" of SyTout is the index "1" of SyTin, opposite to the order in TensorLib
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
	if(inorder  && rsp_outin[0] == 0 && RBondNum == rowBondNum);	//do nothing
	else{
		vector<Bond_t> outBonds;
		for(int b = 0; b < bonds.size(); b++)
			outBonds.push_back(bonds[rsp_outin[b]]);
		for(int b = 0; b < bonds.size(); b++){
			if(b < rowBondNum){
				if(outBonds[b].type == BD_COL)
					for(int q = 0; q < outBonds[b].Qnums.size(); q++)
						outBonds[b].Qnums[q] = -outBonds[b].Qnums[q];
				outBonds[b].type = BD_ROW;
			}
			else{
				if(outBonds[b].type == BD_ROW)
					for(int q = 0; q < outBonds[b].Qnums.size(); q++)
						outBonds[b].Qnums[q] = -outBonds[b].Qnums[q];
				outBonds[b].type = BD_COL;
			}
		}
		SyTensor_t SyTout(outBonds, name);
		if(status & HAVEELEM){
			assert(!(status & ONGPU));
			//For Fermionic system
#ifdef FERMIONIC
			int sign = 1;
			int inLabelF[bondNum];
			int outLabelF[bondNum];
			int ordF[bondNum];
			for(int b = 0; b < RBondNum; b++){
				inLabelF[b] = labels[b];
				ordF[b] = b;
			}
			for(int b = 0; b < SyTout.RBondNum; b++)
				outLabelF[b] = newLabels[b];
			for(int b = bondNum - 1; b >= RBondNum; b--){
				ordF[b] = bondNum - b + RBondNum - 1;
				inLabelF[ordF[b]] = labels[b];
			}
			for(int b = bondNum - 1; b >= SyTout.RBondNum; b--)
				outLabelF[bondNum - b + SyTout.RBondNum - 1] = newLabels[b];

			int rspF_outin[bondNum];
			for(int i = 0; i < bondNum; i++)
				for(int j = 0; j < bondNum; j++)
					if(inLabelF[i] == outLabelF[j])
						rspF_outin[j] = i;	
			vector<_Swap> swaps = _recSwap(rspF_outin, bondNum, ordF);
#endif
			//End Fermionic system
			vector<int> Qin_idxs(bondNum, 0); 
			vector<int> Qot_idxs(bondNum, 0); 
			vector<int> Qin_Bdims(bondNum, 0);
			for(int b = 0; b < bondNum; b++)
				Qin_Bdims[b] = bonds[b].Qnums.size();
			int Qin_off, Qot_off;
			int tmp;
			int Qin_RQoff, Qin_CQoff;
			int Qot_RQoff, Qot_CQoff;
			int sBin_r, sBin_c;	//sub-block of a Qidx
			int sBin_rDim, sBin_cDim;	//sub-block of a Qidx
			int sBot_cDim;	//sub-block of a Qidx
			int sBot_r, sBot_c;
			int Bin_cDim, Bot_cDim;
			int64_t Ein_off, Eot_off;
			vector<int> sBin_idxs(bondNum, 0);
			vector<int> sBin_sBdims(bondNum, 0);
			vector<int> Qot_acc(bondNum, 1); 
			vector<int> sBot_acc(bondNum, 1);
			for(int b = bondNum	- 1; b > 0; b--)
				Qot_acc[b - 1] = Qot_acc[b] * SyTout.bonds[b].Qnums.size();
			for(int enc = 0; enc < QidxNum; enc++){
				Qin_off = Eidx2Qidx[enc];
				tmp = Qin_off;
				for(int b = bondNum - 1; b >= 0; b--){
					Qin_idxs[b] = tmp % Qin_Bdims[b];
					sBin_sBdims[b] = bonds[b].Qdegs[Qin_idxs[b]];
					tmp /= Qin_Bdims[b];
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
				Qot_RQoff = Qot_off / SyTout.CQdim;
				Qot_CQoff = Qot_off % SyTout.CQdim;
				Bin_cDim = RQidx2Blk[Qin_RQoff]->Cnum;
				Bot_cDim = SyTout.RQidx2Blk[Qot_RQoff]->Cnum;
				Ein_off = RQidx2Blk[Qin_RQoff]->offset + (RQidx2Off[Qin_RQoff] * Bin_cDim) + CQidx2Off[Qin_CQoff];
				Eot_off = SyTout.RQidx2Blk[Qot_RQoff]->offset + (SyTout.RQidx2Off[Qot_RQoff] * Bot_cDim) + SyTout.CQidx2Off[Qot_CQoff];
				sBin_rDim = RQidx2Dim[Qin_RQoff];
				sBin_cDim = CQidx2Dim[Qin_CQoff];
				sBot_cDim = SyTout.CQidx2Dim[Qot_CQoff];
				int cnt_ot = 0;
				sBin_idxs.assign(bondNum, 0);
#ifdef FERMIONIC
				int sign01 = 0;
				for(int i = 0; i < swaps.size(); i++)
					sign01 ^= (bonds[swaps[i].b1].Qnums[Qin_idxs[swaps[i].b1]].getPrtF() & bonds[swaps[i].b2].Qnums[Qin_idxs[swaps[i].b2]].getPrtF());
				sign = sign01 ? -1 : 1;
#endif
				for(sBin_r = 0; sBin_r < sBin_rDim; sBin_r++)
					for(sBin_c = 0; sBin_c < sBin_cDim; sBin_c++){
						sBot_r = cnt_ot / sBot_cDim;
						sBot_c = cnt_ot % sBot_cDim;
#ifdef FERMIONIC
						SyTout.elem[Eot_off + (sBot_r * Bot_cDim) + sBot_c] = sign * elem[Ein_off + (sBin_r * Bin_cDim) + sBin_c];
#else
						SyTout.elem[Eot_off + (sBot_r * Bot_cDim) + sBot_c] = elem[Ein_off + (sBin_r * Bin_cDim) + sBin_c];
#endif
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
			SyTout.status |= HAVEELEM;
		}
		*this = SyTout;
		this->addLabel(newLabels);
	}
}

void SyTensor_t::addRawElem(DOUBLE* rawElem){
	assert((status & INIT));   //If not INIT, CANNOT add elements
	assert(!(status & ONGPU));
	int bondNum = bonds.size();
	vector<int> Q_idxs(bondNum, 0); 
	vector<int> Q_Bdims(bondNum, 0);
	vector<int> sB_idxs(bondNum, 0);
	vector<int> sB_sBdims(bondNum, 0);
	vector<int> rAcc(bondNum, 1); 
	for(int b = 0; b < bondNum; b++)
		Q_Bdims[b] = bonds[b].Qnums.size();
	for(int b = bondNum - 1; b > 0; b--)
		rAcc[b - 1] = rAcc[b] * bonds[b].dim;
	int Q_off;
	int tmp;
	int RQoff, CQoff;
	int sB_r, sB_c;	//sub-block of a Qidx
	int sB_rDim, sB_cDim;	//sub-block of a Qidx
	int B_cDim;
	int64_t E_off;
	int64_t R_off;
	for(int enc = 0; enc < QidxNum; enc++){
		Q_off = Eidx2Qidx[enc];
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


double SyTensor_t::at(vector<int> idxs)const{
	int bondNum = bonds.size();
	vector<int> Qidxs(bondNum, 0);
	for(int b = 0; b < bondNum; b++){
		for(int q = bonds[b].offsets.size() - 1; q >= 0; q--){
			if(idxs[b] < bonds[b].offsets[q])
				continue;
			Qidxs[b] = q;
			break;
		}
	}
	vector<int> Q_acc(bondNum, 1); 
	for(int b = bondNum	- 1; b > 0; b--)
		Q_acc[b - 1] = Q_acc[b] * bonds[b].Qnums.size();
	int Qoff = 0;
	for(int b = 0; b < bondNum; b++)
		Qoff += Q_acc[b] * Qidxs[b];
	assert(!(status & ONGPU));
	if(Qidx[Qoff]){
		int Q_RQoff = Qoff / CQdim;
		int Q_CQoff = Qoff % CQdim;
		int B_cDim = RQidx2Blk[Q_RQoff]->Cnum;
		int sB_cDim = CQidx2Dim[Q_CQoff];
		int boff = (RQidx2Blk[Q_RQoff])->offset + (RQidx2Off[Q_RQoff] * B_cDim) + CQidx2Off[Q_CQoff];
		int cnt = 0;
		vector<int> D_acc(bondNum, 1); 
		for(int b = bondNum	- 1; b > 0; b--)
			D_acc[b - 1] = D_acc[b] * bonds[b].Qdegs[Qidxs[b]];
		for(int b = 0; b < bondNum; b++)
			cnt += (idxs[b] - bonds[b].offsets[Qidxs[b]]) * D_acc[b];
		return elem[boff + (cnt / sB_cDim) * B_cDim + cnt % sB_cDim];
	}
	else{
		return 0.0;
	}
}

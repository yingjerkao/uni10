#include "SyTensor.h"

SyTensor_t::SyTensor_t(vector<Bond_t>& _bonds, const string& _name): name(_name){
	//cout<< "Contructing SyTensor " << this << endl;
	assert(_bonds.size() > 0); //No bond in Tensor, Error!
	bonds = _bonds;
	grouping();
	assert(blocks.size() > 0); //No block in Tensor, Error!
	Block_t blk = blocks.rbegin()->second;
	elemNum = blk.offset + (blk.Rnum * blk.Cnum);
	elem = (DOUBLE*)malloc(sizeof(DOUBLE) * elemNum);
	memset(elem, 0, sizeof(DOUBLE) * elemNum);
	status |= INIT;
}
void SyTensor_t::reshape(){
	assert(SyTin->status & INIT);
	assert(SyTin->status & HAVELABE);
	assert(SyTin->labels.size() == newLabels.size());
	int bondNum = SyTin->bonds.size();
	int rsp_outin[bondNum];	//rsp_outin[2] = 1 means the index "2" of SyTout is the index "1" of SyTin, opposite to the order in TensorLib
	int rsp_inout[bondNum];	//rsp_inout[2] = 1 means the index "2" of SyTin is the index "1" of SyTout, the same as the order in TensorLib
	int cnt = 0;
	for(int i = 0; i < bondNum; i++)
		for(int j = 0; j < bondNum; j++)
			if(SyTin->labels[i] == newLabels[j]){
				rsp_outin[j] = i;	
				rsp_inout[i] = j;	
				cnt++;
			}
	assert(cnt == newLabels.size());
	bool inorder = true;
	for(int i = 1; i < bondNum; i++)
		if(((rsp_outin[i] + bondNum - i) % bondNum) != rsp_outin[0]){
			inorder = false;
			break;
		}
	if(inorder && rsp_outin[0] == 0 && SyTin->RBondNum == rowBondNum){	//no need to reshape, clone
			clone(SyTin, SyTout);
			printf("CLONE!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	else if(inorder  && rsp_outin[0] != 0 && (bondNum - SyTin->RBondNum) == rowBondNum){
		transpose(SyTin, SyTout);
		printf("TRANSPOSE!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	else{
		SyTout->bonds.resize(newLabels.size());
		for(int b = 0; b < SyTin->bonds.size(); b++)
			SyTout->bonds[b] = SyTin->bonds[rsp_outin[b]];
		for(int b = 0; b < SyTin->bonds.size(); b++){
			if(b < rowBondNum){
				if(SyTout->bonds[b].type == BD_COL)
					for(int q = 0; q < SyTout->bonds[b].Qnums.size(); q++)
						SyTout->bonds[b].Qnums[q] = -SyTout->bonds[b].Qnums[q];
				SyTout->bonds[b].type = BD_ROW;
			}
			else{
				if(SyTout->bonds[b].type == BD_ROW)
					for(int q = 0; q < SyTout->bonds[b].Qnums.size(); q++)
						SyTout->bonds[b].Qnums[q] = -SyTout->bonds[b].Qnums[q];
				SyTout->bonds[b].type = BD_COL;
			}
		}
		initSyT(SyTout);
		addLabel(SyTout, newLabels);
		if(SyTin->status & HAVEELEM){
			int QcolNum = 1;
			int QcolNum_in = 1;
			for(int b = bondNum - 1; b >= 0; b--)
				if(SyTout->bonds[b].type == BD_COL)
					QcolNum *= SyTout->bonds[b].Qnums.size();
				else
					break;
			for(int b = bondNum - 1; b >= 0; b--)
				if(SyTin->bonds[b].type == BD_COL)
					QcolNum_in *= SyTin->bonds[b].Qnums.size();
				else
					break;
			
			vector<int> Qin_acc(bondNum, 1); 
			vector<int> Dupb(bondNum, 0);
			for(int b = bondNum	- 1; b > 0; b--)
				Qin_acc[b - 1] = Qin_acc[b] * SyTin->bonds[b].Qnums.size();
			int Qoff = 0;	//for SyTout
			int Qoff_in;
			int bend;
			int BcolNum, BcolNum_in;
			int DcolNum, DcolNum_in;
			int boff = 0, boff_in = 0;
			int cnt, cnt_in;
			vector<int> Qidxs(bondNum, 0); 
			vector<int> idxs;
			vector<int> Din_acc(bondNum, 1);	//Degeneracy acc
    
			while(1){
				if(SyTout->Qidx[Qoff]){
					Qoff_in = 0;
					DcolNum = 1;
					DcolNum_in = 1;
					for(int b = 0; b < bondNum; b++){
						Qoff_in += Qidxs[b] * Qin_acc[rsp_outin[b]];
						Dupb[b] = SyTout->bonds[b].Qdegs[Qidxs[b]];
						if(SyTout->bonds[b].type == BD_COL)
							DcolNum *= SyTout->bonds[b].Qdegs[Qidxs[b]];
						if(SyTin->bonds[b].type == BD_COL)
							DcolNum_in *= SyTin->bonds[b].Qdegs[Qidxs[rsp_inout[b]]];
					}
					assert(SyTin->Qidx[Qoff_in]);
					for(int b = bondNum	- 1; b > 0; b--)
						Din_acc[b - 1] = Din_acc[b] * SyTin->bonds[b].Qdegs[Qidxs[rsp_inout[b]]];
					BcolNum = (SyTout->RQidx2Blk[Qoff / QcolNum])->Cnum;
					boff = (SyTout->RQidx2Blk[Qoff / QcolNum])->offset + SyTout->RQidx2Off[Qoff / QcolNum] * BcolNum + SyTout->CQidx2Off[Qoff % QcolNum];
					BcolNum_in = (SyTin->RQidx2Blk[Qoff_in / QcolNum_in])->Cnum;
					boff_in = (SyTin->RQidx2Blk[Qoff_in / QcolNum_in])->offset + SyTin->RQidx2Off[Qoff_in / QcolNum_in] * BcolNum_in + SyTin->CQidx2Off[Qoff_in % QcolNum_in];
					cnt = 0;
					cnt_in = 0;
					idxs.assign(bondNum, 0);
					while(1){
						SyTout->elem[boff + (cnt / DcolNum) * BcolNum + cnt % DcolNum] = SyTin->elem[boff_in + (cnt_in / DcolNum_in) * BcolNum_in + cnt_in % DcolNum_in];
						cnt++;
						for(bend = bondNum - 1; bend >= 0; bend--){
							idxs[bend]++;
							if(idxs[bend] < Dupb[bend]){
								cnt_in += Din_acc[rsp_outin[bend]];
								break;
							}
							else{
								cnt_in -= Din_acc[rsp_outin[bend]] * (idxs[bend] - 1);
								idxs[bend] = 0;
							}
						}
						if(bend < 0)
							break;
					}
				}
				for(bend = bondNum - 1; bend >= 0; bend--){
					Qidxs[bend]++;
					if(Qidxs[bend] < SyTout->bonds[bend].Qnums.size())
						break;
					else
						Qidxs[bend] = 0;
				}
				Qoff++;
				if(bend < 0)
					break;
			}
			SyTout->status |= HAVEELEM;
		}
	}
	
}

void SyTensor_t::grouping(){
	blocks.clear();
	int row_bondNum = 0;
	int col_bondNum = 0;
	int row_Qdim = 1;
	int col_Qdim = 1;
	for(int i = 0; i < bonds.size(); i++){
		if(bonds[i].type == BD_ROW){
			row_Qdim *= bonds[i].Qnums.size();
			row_bondNum++;
		}
		else{
			col_Qdim *= bonds[i].Qnums.size();
			col_bondNum++;
		}
	}
	RBondNum = row_bondNum;
	map<Qnum_t,int> row_QnumMdim;
	vector<int> row_offs(row_bondNum, 0);
	map<Qnum_t,vector<int> > row_Qnum2Qidx;
	Qnum_t qnum(0, 0);
	int dim;
	int boff = 0;
	RQidx2Off.assign(row_Qdim, 0);
	CQidx2Off.assign(col_Qdim, 0);
	RQidx2Blk.assign(row_Qdim, NULL);
	Qidx.assign(row_Qdim * col_Qdim, false);
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
				row_QnumMdim[qnum] += dim;
			}
			else{
				RQidx2Off[boff] = 0;
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
	}
	map<Qnum_t,int>::iterator itt;
	for ( itt = row_QnumMdim.begin() ; itt != row_QnumMdim.end(); itt++ ){
		cout<< itt->first << "  ";
		cout<< itt->second;
	}

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
					col_QnumMdim[qnum] += dim;
				}
				else{
					CQidx2Off[boff] = 0;
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
		}
	}

	map<Qnum_t,int>::iterator it;
	map<Qnum_t,int>::iterator it2;
	Block_t blk(0, 0, 0);
	int off = 0;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		blk.Rnum = it->second;
		blk.Cnum = it2->second;
		blk.offset = off;
		off += blk.Rnum * blk.Cnum;
		blocks[it->first] = blk;
		Block_t* blkptr= &(blocks[it->first]);
		vector<int>& tmpRQidx = row_Qnum2Qidx[it->first];
		vector<int>& tmpCQidx = col_Qnum2Qidx[it->first];
		for(int i = 0; i < tmpRQidx.size(); i++){
			RQidx2Blk[tmpRQidx[i]] = blkptr;
			for(int j = 0; j < tmpCQidx.size(); j++)
				Qidx[tmpRQidx[i] * col_Qdim + tmpCQidx[j]] = true;
		}
	}
}


ostream& operator<< (ostream& os, SyTensor_t& SyT){
	assert(SyT.status & INIT);
	int row = 0;
	int col = 0;
	for(int i = 0; i < SyT.bonds.size(); i++)
		if(SyT.bonds[i].type == BD_ROW)
			row++;
		else
			col++;
	int layer = max(row, col);
	int nmlen = SyT.name.length() + 2;
	int star = 12 + (14 - nmlen) / 2;
	cout<<endl;
	for(int s = 0; s < star; s++)
		cout << "*";
	cout << " " << SyT.name << " ";
	for(int s = 0; s < star; s++)
		cout<<"*";	
	cout << "\n             ____________\n";
	cout << "            |            |\n";
	int llab = 0;
	int rlab = 0;
	char buf[128];
	for(int l = 0; l < layer; l++){
		if(l < row && l < col){
			if(SyT.status & HAVELABEL){
				llab = SyT.labels[l];
				rlab = SyT.labels[row + l];
			}
			else{
				llab = l;
				rlab = row + l;
			}
			sprintf(buf, "    %5d___|%-4d    %4d|___%-5d\n", llab, SyT.bonds[l].dim, SyT.bonds[row + l].dim, rlab);
			cout<<buf;
		}
		else if(l < row){
			if(SyT.status & HAVELABEL)
				llab = SyT.labels[l];
			else
				llab = l;
			sprintf(buf, "    %5d___|%-4d    %4s|\n", llab, SyT.bonds[l].dim, "");
			cout<<buf;
		}
		else if(l < col){
			if(SyT.status & HAVELABEL)
				rlab = SyT.labels[row + l];
			else
				rlab = row + l;
			sprintf(buf, "    %5s   |%4s    %4d|___%-5d\n", "", "", SyT.bonds[row + l].dim, rlab);
			cout << buf;
		}   
		cout << "            |            |   \n";
	}   
	cout << "            |____________|\n";

	cout << "\n================BONDS===============\n";
	for(int b = 0; b < SyT.bonds.size(); b++){
		cout << SyT.bonds[b];
	}   
	cout<<"\n===============BLOCKS===============\n";
	map<Qnum_t,Block_t>::iterator it; 
	int Rnum, Cnum;
	bool printElem = true;
	for ( it = SyT.blocks.begin() ; it != SyT.blocks.end(); it++ ){
		cout << "--- " << it->first << it->second;
		Rnum = it->second.Rnum;
		Cnum = it->second.Cnum;
		cout << " ---\n\n";

		if(SyT.status & HAVEELEM && printElem){
			for(int r = 0; r < Rnum; r++){
				for(int c = 0; c < Cnum; c++)
					cout<< setw(6) << setprecision(3) << SyT.elem[it->second.offset + r * Cnum + c];
				cout << "\n\n";
			}
		}
	}   
	cout << "Total elemNum: "<<SyT.elemNum<<endl;
	cout << "***************** END ****************\n\n";
	return os;
}

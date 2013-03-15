
void addRawElem(SyTensor_t* SyT, DOUBLE* rawElem){
	assert((SyT->status & INIT));   //If not INIT, CANNOT add elements
	int bondNum = SyT->bonds.size();
	vector<int> Qidxs(bondNum, 0); 
	vector<int> idxs(bondNum, 0); 
	vector<int> idxUpb(bondNum, 0); //upper bound of a index of some block
	vector<int> rAcc(bondNum, 1); 
	map<Qnum_t,Block_t>::iterator git;
	map<Qnum_t,int> grp_cnt;
	for (git = SyT->blocks.begin(); git != SyT->blocks.end(); git++)
		grp_cnt[git->first] = git->second.offset;
	for(int b = bondNum - 1; b > 0; b--)
		rAcc[b - 1] = rAcc[b] * SyT->bonds[b].dim;
	int bend = 0;
	int64_t roff = 0;
	int64_t boff = 0;
	int Qoff = 0;
	int QcolNum = 1;
	int BcolNum = 0;
	int DcolNum = 0;	//Degeneracy column number
	int cnt;
	for(int b = bondNum - 1; b >= 0; b--)
		if(SyT->bonds[b].type == BD_COL)
			QcolNum *= SyT->bonds[b].Qnums.size();
		else
			break;
	while(1){
		if(SyT->Qidx[Qoff]){
			roff = 0;
			DcolNum = 1;
			for(int b = 0; b < bondNum; b++){
				idxs[b] = SyT->bonds[b].offsets[Qidxs[b]];
				idxUpb[b] = idxs[b] + SyT->bonds[b].Qdegs[Qidxs[b]];
				roff += idxs[b] * rAcc[b];
				if(SyT->bonds[b].type == BD_COL)
					DcolNum *= SyT->bonds[b].Qdegs[Qidxs[b]];
			}

			BcolNum = (SyT->RQidx2Blk[Qoff / QcolNum])->Cnum;
			boff = (SyT->RQidx2Blk[Qoff / QcolNum])->offset + SyT->RQidx2Off[Qoff / QcolNum] * BcolNum + SyT->CQidx2Off[Qoff % QcolNum];

			cnt = 0;
			while(1){
				SyT->elem[boff + (cnt / DcolNum) * BcolNum + cnt % DcolNum] = rawElem[roff];
				for(bend = bondNum - 1; bend >= 0; bend--){
					idxs[bend]++;
					if(idxs[bend] < idxUpb[bend]){
						roff += rAcc[bend];
						break;
					}
					else{
						roff -= rAcc[bend] * (idxs[bend] - SyT->bonds[bend].offsets[Qidxs[bend]] - 1);
						idxs[bend] = SyT->bonds[bend].offsets[Qidxs[bend]];
					}
				}
				if(bend < 0)
					break;
				cnt++;
			}
		}   
		for(bend = bondNum - 1; bend >= 0; bend--){
			Qidxs[bend]++;
			if(Qidxs[bend] < SyT->bonds[bend].Qnums.size())
				break;
			else
				Qidxs[bend] = 0;
		}
		Qoff++;
		if(bend < 0)
			break;
	}
	SyT->status |= HAVEELEM;
}


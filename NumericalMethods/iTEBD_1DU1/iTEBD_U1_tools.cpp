void bondcat(UniTensor& T, UniTensor& L, int bidx);
void bondrm(UniTensor& T, UniTensor& L, int bidx);
Bond mkBond(bondType btype, map<Qnum, int>& trunc);
void setTruncation(UniTensor& theta, UniTensor& GA, UniTensor& GB, UniTensor& LA, size_t chi);
void sv_merge(vector<double>& svs, vector<size_t>& bidxs, size_t bidx, Matrix& sv_mat, size_t chi);

void bondcat(UniTensor& T, UniTensor& L, int bidx){
  vector<int> labels = T.label();
  size_t inBondNum = T.inBondNum();
  if(bidx == 0){
    /* -- L -- T --
     *         |
     */
    int labelL[] = {labels[0], 77};
    int labelT[] = {77, labels[1], labels[2]};
    L.setLabel(labelL);
    T.setLabel(labelT);
    T = L * T;
  }
  else if(bidx == 1){
    /* -- T -- L --
     *    |
     */
    int labelL[] = {77, labels[1]};
    int labelT[] = {labels[0], 77, labels[2]};
    L.setLabel(labelL);
    T.setLabel(labelT);
    T *= L;
  }
  else{
    throw runtime_error("\nIn bondcat(): Invalid Usage.");
  }
  T.permute(labels, inBondNum);
}

void bondrm(UniTensor& T, UniTensor& L, int bidx){
	UniTensor invL(L);
  double *elem = invL.getElem();
  for(int i = 0; i < invL.elemNum(); i++)
      elem[i] = elem[i] == 0 ? 0 : (1 / elem[i]);
	bondcat(T, invL, bidx);
}

Bond mkBond(bondType btype, map<Qnum, int>& trunc){
  for(map<Qnum, int>::iterator it = trunc.begin(); it != trunc.end(); it++){
    map<Qnum, int>::iterator nit = trunc.find(-(it->first));
    if(nit != trunc.end()){
      int dim = it->second > nit->second ? it->second : nit->second;
      it->second = dim;
      nit->second = dim;
    }
    else{
      trunc[-(it->first)] = it->second;
    }
  }
  vector<Qnum> qnums;
  for(map<Qnum, int>::iterator it = trunc.begin(); it != trunc.end(); it++)
    qnums.insert(qnums.end(), it->second, it->first);
  return Bond(btype, qnums);
}

void setTruncation(UniTensor& theta, UniTensor& GA, UniTensor& GB, UniTensor& LA, size_t chi){
  map<Qnum, vector<Matrix> >svds;
  vector<Qnum> blk_qnums = theta.blockQnum();
  for(vector<Qnum>::iterator q = blk_qnums.begin(); q != blk_qnums.end(); q++)
    svds[*q] = theta.getBlock(*q).svd();
  vector<double> svs;
  vector<size_t> bidxs;
  for(size_t bidx = 0; bidx < blk_qnums.size(); bidx++)
    sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi);
  vector<size_t>dims(blk_qnums.size(), 0);
  for(int i = 0; i < bidxs.size(); i++)
    dims[bidxs[i]]++;

  vector<Qnum> qnums(svs.size());
  int cnt = 0;
  for(size_t bidx = 0; bidx < blk_qnums.size(); bidx++){
    for(int i = 0; i < dims[bidx]; i++){
      qnums[cnt] = blk_qnums[bidx];
      cnt++;
    }
  }
  Bond bdi_mid = Bond(BD_IN, qnums);
  Bond bdo_mid = Bond(BD_OUT, qnums);

  vector<Bond> bondGa = GA.bond();
  vector<Bond> bondGb = GB.bond();
  vector<Bond> bondLa;
  bondGa[2] = bdo_mid;
  bondGb[0] = bdi_mid;
  bondLa.push_back(bdi_mid);
  bondLa.push_back(bdo_mid);
  vector<int> labelGa = GA.label();
  GA.assign(bondGa);
  GB.assign(bondGb);
  LA.assign(bondLa);
  GA.setLabel(labelGa);
  map<Qnum, int> degs = bdi_mid.degeneracy();
  Matrix sv_mat(bdi_mid.dim(), bdo_mid.dim(), svs, true);
  double norm = sv_mat.norm();
  for(map<Qnum, int>::iterator it = degs.begin(); it != degs.end(); it++){
    map<Qnum, vector<Matrix> >::iterator sit = svds.find(it->first);
    if(sit == svds.end())
      throw runtime_error("\nFatal Error in setTruncation()");
    GA.putBlock(it->first, sit->second[0].resize(sit->second[0].row(), it->second));
    GB.putBlock(it->first, sit->second[2].resize(it->second, sit->second[2].col()));
    LA.putBlock(it->first, sit->second[1].resize(it->second, it->second) * (1/norm));
  }
}

void sv_merge(vector<double>& svs, vector<size_t>& bidxs, size_t bidx, Matrix& sv_mat, size_t chi){
  if(svs.size()){
    int len = svs.size() + sv_mat.elemNum();
    len = len < chi ? len : chi;
    vector<double> ori_svs = svs;
    vector<size_t> ori_bidxs = bidxs;
    svs.assign(len, 0);
    bidxs.assign(len, 0);
    int cnt = 0;
    int cur1 = 0;
    int cur2 = 0;
    while(cnt < len){
      if(cur1 < ori_svs.size() && cur2 < sv_mat.elemNum()){
        if(ori_svs[cur1] >= sv_mat[cur2]){
          svs[cnt] = ori_svs[cur1];
          bidxs[cnt] = ori_bidxs[cur1];
          cur1++;
        }
        else{
          svs[cnt] = sv_mat[cur2];
          bidxs[cnt] = bidx;
          cur2++;
        }
      }
      else if(cur2 < sv_mat.elemNum()){
        for(; cnt < svs.size(); cnt++){
          svs[cnt] = sv_mat[cur2];
					bidxs[cnt] = bidx;
					cur2++;
				}
				break;
			}
			else{
				for(; cnt < svs.size(); cnt++){
					svs[cnt] = ori_svs[cur1];
					bidxs[cnt] = ori_bidxs[cur1];
					cur1++;
				}
				break;
			}
			cnt++;
    }
  }
  else{
    bidxs.assign(sv_mat.elemNum(), bidx);
    svs.assign(sv_mat.elemNum(), 0);
    for(int i = 0; i < sv_mat.elemNum(); i++)
      svs[i] = sv_mat[i];
  }
}

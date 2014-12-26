void bondcat(UniTensor& T, const Matrix& L, int bidx);
void bondrm(UniTensor& T, const Matrix& L, int bidx);
UniTensor triHamiltonian(const UniTensor& H0);
UniTensor truncateCore(const UniTensor& C, size_t chi);
double measureObs(bool upT, const UniTensor& Ob, vector<UniTensor> Us, const vector<UniTensor>& Cs, const vector<Matrix>& Ls, Network& state, Network& measure);

void bondcat(UniTensor& T, const Matrix& L, int bidx){
  int inBondNum = T.inBondNum();
	vector<int> labels = T.label();
  vector<int> per_labels = labels;
  int l = labels[bidx];
  per_labels.erase(per_labels.begin() + bidx);
	per_labels.insert(per_labels.begin(), l);
  T.permute(per_labels, 1);
  T.putBlock(L * T.getBlock());
  T.permute(labels, inBondNum);

}

void bondrm(UniTensor& T, const Matrix& L, int bidx){
	Matrix invL = L;
  for(int i = 0; i < L.elemNum(); i++)
    invL[i] = invL[i] == 0 ? 0 : (1 / invL[i]);
	bondcat(T, invL, bidx);
}

UniTensor triHamiltonian(const UniTensor& H0){
  vector<Bond> bondI;
  bondI.push_back(H0.bond(0)), bondI.push_back(H0.bond(2));
  UniTensor I(bondI);
  I.identity();

  UniTensor H3 = otimes(H0, I) + otimes(I, H0);
  int labels[] = {0, 1, 2, 3, 4, 5};
  H3.setLabel(labels);
  int per_labels[] = {2, 0, 1, 5, 3, 4};
  H3.permute(per_labels, 3);
  H3 += otimes(H0, I);
  H3.setLabel(labels);
  return H3;
}

UniTensor truncateCore(const UniTensor& C, size_t chi){
  UniTensor truC = C;
  int initLabel[] = {0, 1, 2};
  truC.setLabel(initLabel);
  int coreBondNum = truC.bondNum();
  Bond bdi_chi(BD_IN, chi);
  truC.permute(1);
  for(int i = 0; i < coreBondNum; i++){
    vector<Bond> bonds = truC.bond();
    bonds[0] = bdi_chi;
    UniTensor newC(bonds);
    newC.putBlock(truC.getBlock().resize(chi, bonds[1].dim() * bonds[2].dim()));
    truC = newC;
    int per_label[3] = {1, 2, 0};
    truC.permute(per_label, 1);
  }
  truC.setLabel(C.label());
  return truC;
}

void simpleUpdate(bool upT, vector<UniTensor>& Us, vector<UniTensor>& Cs, vector<Matrix>& Ls, const UniTensor& expH, Network& simplex){
  for(int u = 0; u < Us.size(); u++)
    if(upT)
      bondcat(Us[u], Ls[u + 3], 2);
    else
      bondcat(Us[u], Ls[u], 1);
  simplex.putTensor("UA", Us[0]);
  simplex.putTensor("UB", Us[1]);
  simplex.putTensor("UC", Us[2]);
  simplex.putTensor("expH", expH);
  if(upT)
    simplex.putTensor("Cup", Cs[0]);
  else
    simplex.putTensor("Cdn", Cs[1]);
  UniTensor theta = simplex.launch();


  vector<Matrix> svdLs;
  vector<UniTensor> svdUs = theta.hosvd(3, svdLs);
  /*======================= Truncation ======================*/
  int chi = Cs[0].bond(0).dim();
  int d = Us[0].bond(0).dim();
  //Truncate core
  UniTensor trunC = truncateCore(svdUs[3], chi);
  double norm = trunC.getBlock().norm();
  if(upT)
    Cs[0] = trunC * (1/norm);
  else
    Cs[1] = trunC * (1/norm);

  //Truncate Lamda
  for(int i = 0; i < svdLs.size(); i++){
    if(upT)
      Ls[i] = svdLs[i].resize(chi, chi) * (1/norm);
    else
      Ls[3 + i] = svdLs[i].resize(chi, chi) * (1/norm);
  }

  //Truncate Us
  if(upT){
    vector<UniTensor> trunUs(3, UniTensor(Us[0].bond()));
    for(int i = 0; i < trunUs.size(); i++)
      trunUs[i].putBlock(svdUs[i].getBlock().resize(d * chi, chi));
    for(int i = 0; i < Us.size(); i++){
      int per_label[] = {0, 2, 1};
      Us[i] = trunUs[i].permute(per_label, 2);
    }
  }
  else
    for(int i = 0; i < Us.size(); i++)
      Us[i].putBlock(svdUs[i].getBlock().resize(d * chi, chi));
  /*============End of Truncation ================*/
  for(int u = 0; u < Us.size(); u++)
    if(upT)
      bondrm(Us[u], Ls[u + 3], 2);
    else
      bondrm(Us[u], Ls[u], 1);
}

double measureObs(bool upT, const UniTensor& Ob, vector<UniTensor> Us, const vector<UniTensor>& Cs, const vector<Matrix>& Ls, Network& state, Network& measure){
  for(int u = 0; u < Us.size(); u++){
    if(upT)
      bondcat(Us[u], Ls[u + 3], 2);
    else
      bondcat(Us[u], Ls[u], 1);
  }

  if(upT)
    state.putTensor("Cup", Cs[0]);
  else
    state.putTensor("Cdn", Cs[1]);
  state.putTensor("UA", Us[0]);
  state.putTensor("UB", Us[1]);
  state.putTensor("UC", Us[2]);
  UniTensor S = state.launch();
  measure.putTensor("S", S);
  measure.putTensorT("ST", S);
  measure.putTensor("Ob", Ob);
  return measure.launch()[0] / S.getBlock().norm();
}



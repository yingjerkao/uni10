void bondcat(UniTensor& T, Matrix& L, int bidx);
void bondrm(UniTensor& T, Matrix& L, int bidx);

void bondcat(UniTensor& T, Matrix& L, int bidx){
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

void bondrm(UniTensor& T, Matrix& L, int bidx){
	Matrix invL(L.row(), L.col(), L.isDiag());
  for(int i = 0; i < L.elemNum(); i++)
    invL[i] = L[i] == 0 ? 0 : (1 / L[i]);
	bondcat(T, invL, bidx);
}

void bondcat(UniTensor& T, Matrix& L, int bidx);
void bondrm(UniTensor& T, Matrix& L, int bidx);
void updateU(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
void updateR(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
void updateD(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
void updateL(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA);
double measure2(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD);
double measure(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& Ob, Network& measure_net, Network& norm_net);
double bond_expectation(UniTensor& ALL, UniTensor& BLL, UniTensor& Ob, Network& measure_net, Network& norm_net);

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

void updateU(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	bondcat(BL, LR, 4);
	iTEBD.putTensor("A", AL);
	iTEBD.putTensor("B", BL);
	iTEBD.putTensor("expH", expH);
  UniTensor C = iTEBD.launch();
	UniTensor Theta = C;
	bondcat(Theta, LD, 2);
	bondcat(Theta, LL, 3);

	vector<Matrix> svds = Theta.getBlock().svd();
	int dim = LU.row();
	LU = svds[1].resize(dim, dim);
	double norm = LU.norm();
	LU *= (1 / norm);
	int per_order[] = {3, 0, 4, 1, 2};
  BL.assign(BL.bond());
  BL.permute(per_order, 1);
	BL.putBlock(svds[2].resize(dim, svds[2].col()));
	int order[] = {0, 1, 2, 3, 4};
	BL.permute(order, 1);

	updateA.putTensor("B", BL);
	updateA.putTensor("C", C);
	AL = updateA.launch();
	AL *= (1 / norm);
  AL.setLabel(order);
	bondrm(BL, LR, 4);
}

void updateR(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	bondcat(BL, LU, 3);
	iTEBD.putTensor("A", AL);
	iTEBD.putTensor("B", BL);
	iTEBD.putTensor("expH", expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta = C;
	bondcat(Theta, LD, 1);
	bondcat(Theta, LL, 2);

	vector<Matrix> svds = Theta.getBlock().svd();
  int dim = LR.row();
	LR = svds[1].resize(dim, dim);
	double norm = LR.norm();
	LR *= (1 / norm);
	int per_order[] = {4, 0, 1, 2, 3};
	BL.assign(BL.bond());
  BL.permute(per_order, 1);
  BL.putBlock(svds[2].resize(dim, svds[2].col()));
	int order[] = {0, 1, 2, 3, 4};
	BL.permute(order, 1);

	updateA.putTensor("B", BL);
	updateA.putTensor("C", C);
	AL = updateA.launch();
	AL *= (1 / norm);
  AL.setLabel(order);
	bondrm(BL, LU, 3);
}

void updateD(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	updateU(BL, AL, LD, LL, LU, LR, expH, iTEBD, updateA);
}

void updateL(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD, Network& updateA){
	updateR(BL, AL, LD, LL, LU, LR, expH, iTEBD, updateA);
}

double measure2(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& expH, Network& iTEBD){
	UniTensor A = AL;
	UniTensor B = BL;
	bondcat(B, LR, 4);
	iTEBD.putTensor("A", &A);
	iTEBD.putTensor("B", &B);
	iTEBD.putTensor("expH", &expH);
	UniTensor C = iTEBD.launch();
	UniTensor Theta = C;
	bondcat(Theta, LD, 2);
	bondcat(Theta, LL, 3);

  UniTensor Theta2 = Theta;
  UniTensor val = Theta * Theta2;
  double delta = 0.1;
  return -log(val[0])/delta/2;
}

double measure(UniTensor& AL, UniTensor& BL, Matrix& LU, Matrix& LR, Matrix& LD, Matrix& LL, UniTensor& Ob, Network& measure_net, Network& norm_net){
  UniTensor ALL = AL;
  UniTensor BLL = BL;
	bondcat(ALL, LD, 3);
	bondcat(ALL, LL, 4);
	bondcat(BLL, LU, 3);
	bondcat(BLL, LR, 4);

  double val = 0;

  // measure up bond
  UniTensor BLL1 = BLL;
  bondrm(BLL1, LU, 3);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);

  // measure right bond
  BLL1 = BLL;
  int rotateR_label[] = {0, 2, 3, 4, 1};
  ALL.permute(rotateR_label, 1);
  BLL1.permute(rotateR_label, 1);
  bondrm(BLL1, LR, 3);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);

  // measure down bond
  BLL1 = BLL;
  int rotateD_label[] = {0, 3, 4, 1, 2};
  ALL.permute(rotateD_label, 1);
  BLL1.permute(rotateD_label, 1);
  bondrm(BLL1, LD, 3);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);

  // measure down bond
  BLL1 = BLL;
  int rotateL_label[] = {0, 4, 1, 2, 3};
  ALL.permute(rotateL_label, 1);
  BLL1.permute(rotateL_label, 1);
  bondrm(BLL1, LL, 3);
  val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);
  return val / 4;
}

double bond_expectation(UniTensor& ALL, UniTensor& BLL, UniTensor& Ob, Network& measure_net, Network& norm_net){
	measure_net.putTensor("ALL", ALL);
	measure_net.putTensor("BLL", BLL);
	measure_net.putTensorT("ALLT", ALL);
	measure_net.putTensorT("BLLT", BLL);
	measure_net.putTensorT("Ob", Ob);
	UniTensor val = measure_net.launch();

	norm_net.putTensor("ALL", ALL);
	norm_net.putTensor("BLL", BLL);
	norm_net.putTensorT("ALLT", ALL);
	norm_net.putTensorT("BLLT", BLL);
	UniTensor norm = norm_net.launch();
  return val[0] / norm[0];
}


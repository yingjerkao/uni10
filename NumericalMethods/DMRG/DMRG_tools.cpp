enum Side{Left = -1, Center = 0, Right = 1};
UniTensor combineH(const UniTensor& H0, const UniTensor& HL, const UniTensor& HR);
UniTensor findGS(const UniTensor& SB, double& E0, Matrix& refState, int& iter);
int updateMPS(const UniTensor& GS, size_t chi, UniTensor& A, UniTensor& B, Matrix& lambda);
int updateMPS(const UniTensor& GS, size_t chi, UniTensor& A, UniTensor& B);
double sweep(int N, int chi, int range, int times, UniTensor& H0, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn);
double sweep(int N, int chi, int range, int times, vector<UniTensor>& H0s, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn);
double sweep(int N, int chi, int range, int times, UniTensor& H0, vector<UniTensor>& HLs, vector<UniTensor>& HRs, vector<UniTensor>& As, vector<UniTensor>& Bs, vector<Matrix>& Ls, Network& HLn, Network& HRn);
double sweep(int N, int chi, int range, int times, vector<UniTensor>& H0s, vector<UniTensor>& HLs, vector<UniTensor>& HRs, vector<UniTensor>& As, vector<UniTensor>& Bs, vector<Matrix>& Ls, Network& HLn, Network& HRn);
void bondcat(UniTensor& T, const Matrix& L, int bidx);
void bondrm(UniTensor& T, const Matrix& L, int bidx);
Matrix makeMPS(Side side, const UniTensor& A, const UniTensor& B, Matrix& L);
Matrix trialState(const UniTensor& A, const UniTensor& B, vector<Matrix>& Ls);
size_t hidx(int N, Side side, int l);

UniTensor combineH(const UniTensor& H0, const UniTensor& HL, const UniTensor& HR){
	vector<Bond> bond2;
  bond2.push_back(HL.bond(0));
  bond2.push_back(HL.bond(2));
	UniTensor IDL(bond2);
	IDL.identity();
  UniTensor IDR(bond2);
  bond2.clear();
  bond2.push_back(HR.bond(1));
  bond2.push_back(HR.bond(3));
  IDR.assign(bond2);
  IDR.identity();
	UniTensor IDd(HL.bond());
	IDd.identity();
	UniTensor IdD(HR.bond());
	IdD.identity();
	UniTensor SB = otimes(HL, IdD);
  SB += otimes(otimes(IDL, H0), IDR);
  SB += otimes(IDd, HR);
  return SB;
}

UniTensor findGS(const UniTensor& SB, double& E0, Matrix& refState, int& iter){
  Matrix blk = SB.getBlock();
  if(refState.col() != blk.col()){
    refState.resize(1, blk.col());
    refState.randomize();
  }
  iter= blk.lanczosEigh(E0, refState, 1000);

  vector<Bond> stateBonds;
  for(int b = 4; b < 8; b++)
    stateBonds.push_back(SB.bond(b));
	UniTensor GS(stateBonds, "GS");
  GS.setElem(refState.getElem());
  return GS.permute(2);
}

int updateMPS(const UniTensor& GS, size_t chi, UniTensor& A, UniTensor& B){
  Matrix lambda;
  return updateMPS(GS, chi, A, B, lambda);
}
int updateMPS(const UniTensor& GS, size_t chi, UniTensor& A, UniTensor& B, Matrix& lambda){
  int DL = GS.bond(0).dim() * GS.bond(1).dim();
  int DR = GS.bond(2).dim() * GS.bond(3).dim();
  int D = DL > DR ? DL : DR;
	D = D < chi ? D : chi;
	Bond bDi(BD_IN, D);
  Bond bDo(BD_OUT, D);

  vector<Bond> bondA;
  bondA.push_back(GS.bond(0));
  bondA.push_back(GS.bond(1));
  bondA.push_back(bDo);
  vector<Bond> bondB;
  bondB.push_back(bDi);
  bondB.push_back(GS.bond(2));
  bondB.push_back(GS.bond(3));
  vector<Matrix> svd = GS.getBlock().svd();

	A.assign(bondA);
  A.putBlock(svd[0].resize(bondA[0].dim() * bondA[1].dim(), bDo.dim()));
	B.assign(bondB);
  B.putBlock(svd[2].resize(bDi.dim(), bondB[1].dim() * bondB[2].dim()));
  lambda = svd[1].resize(chi, chi);
  return D;
}

void updateH(const UniTensor& HL, const UniTensor& HR, const UniTensor& A, const UniTensor& B, const UniTensor& H0L, const UniTensor H0R, Network& HLn, Network& HRn, UniTensor& newHL, UniTensor& newHR){
	vector<Bond> bond2;
	bond2.push_back(H0L.bond(0));
	bond2.push_back(H0L.bond(2));
	UniTensor Id(bond2, "Id");
	Id.identity();

  bond2.clear();
  bond2.push_back(HL.bond(0));
  bond2.push_back(HL.bond(2));
  UniTensor IDL(bond2);
  IDL.identity();
  bond2.clear();
  bond2.push_back(HR.bond(1));
  bond2.push_back(HR.bond(3));
  UniTensor IDR(bond2);
  IDR.identity();

  UniTensor HL2 = otimes(HL, Id);
  HL2 += otimes(IDL, H0L);
  UniTensor HR2 = otimes(Id, HR);
  HR2 += otimes(H0R, IDR);

  HLn.putTensor("Al", &A);
  HLn.putTensor("HL2", &HL2);
  HLn.putTensorT("AlT", &A);

  HRn.putTensor("Bl", &B);
  HRn.putTensor("HR2", &HR2);
  HRn.putTensorT("BlT", &B);
  newHL = HLn.launch();
  newHR = HRn.launch();
}

double sweep(int N, int chi, int range, int times, UniTensor& H0, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn){
  vector<UniTensor> H0s(1, H0);
  return sweep(N, chi, range, times, H0s, HLs, HRs, HLn, HRn);
}

double sweep(int N, int chi, int range, int times, vector<UniTensor>& H0s, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn){
  vector<UniTensor> As;
  vector<UniTensor> Bs;
  vector<Matrix> Ls;
  return sweep(N, chi, range, times, H0s, HLs, HRs, As, Bs, Ls, HLn, HRn);
}
double sweep(int N, int chi, int range, int times, UniTensor& H0, vector<UniTensor>& HLs, vector<UniTensor>& HRs, vector<UniTensor>& As, vector<UniTensor>& Bs, vector<Matrix>& Ls, Network& HLn, Network& HRn){
  vector<UniTensor> H0s(1, H0);
  return sweep(N, chi, range, times, H0s, HLs, HRs, As, Bs, Ls, HLn, HRn);
}

double sweep(int N, int chi, int range, int times, vector<UniTensor>& H0s, vector<UniTensor>& HLs, vector<UniTensor>& HRs, vector<UniTensor>& As, vector<UniTensor>& Bs, vector<Matrix>& Ls, Network& HLn, Network& HRn){
  assert(range < N);
  Matrix psi;
  int dir = Left;
  int cursor = 0;
  int cnt = 0;
  double E0;
  bool MPS_READY = false;
  if(As.size() == N && Bs.size() == N && Ls.size() == N){
    MPS_READY = true;
    for(int l = N - 2; l >= 0; l--) //Ls.size() is 2N-1
      Ls.push_back(Ls[l]);
  }
  while(cnt < times + 1){
    UniTensor* Hptr;
    for( ; abs(cursor) < range; cursor += dir){
      if(H0s.size() > 1)
        Hptr = &(H0s[hidx(N, Left, N + cursor - 1)]);
      else
        Hptr = &(H0s[0]);
      UniTensor SB = combineH(*Hptr, HLs[N + cursor - 2], HRs[N - cursor - 2]);
      if(MPS_READY){
        if(cursor < 0)
          psi = makeMPS(Left, As[N + cursor - 1], As[N + cursor], Ls[hidx(N, Left, N + cursor)]);
        else if(cursor > 0)
          psi = makeMPS(Right, Bs[N - cursor], Bs[N - cursor - 1], Ls[hidx(N, Left, N + cursor - 2)]);
        else
          psi = makeMPS(Center, As[N - 1], Bs[N - 1], Ls[hidx(N, Left, N + cursor - 1)]);
      }

      int iter;
      UniTensor GS = findGS(SB, E0, psi, iter);

      UniTensor A, B;
      Matrix L;
      int D = updateMPS(GS, chi, A, B, L);
      if(MPS_READY){
        if(cursor < 0){
          As[N + cursor - 1] = A;
          As[N + cursor] = B;
          bondcat(As[N + cursor], L, 0);
          As[N + cursor].permute(2);
          bondrm(As[N + cursor], Ls[hidx(N, Left, N + cursor)], 2);
          Ls[hidx(N, Left, N + cursor - 1)] = L;
        }
        else if(cursor > 0){
          Bs[N - cursor] = A;
          Bs[N - cursor - 1] = B;
          Bs[N - cursor].permute(1);
          bondcat(Bs[N - cursor], L, 2);
          bondrm(Bs[N - cursor], Ls[hidx(N, Left, N + cursor - 2)], 0);
          Ls[hidx(N, Left, N + cursor - 1)] = L;
        }
        else{
          As[N - 1] = A;
          Bs[N - 1] = B;
          Ls[hidx(N, Left, N + cursor - 1)] = L;
        }
      }
      cout<<"cursor = "<< cursor <<", D = " << chi << setprecision(10) << ", E = " << E0  << ", e = " << E0 / (2 * N) <<", iter = "<<iter<<endl;

      UniTensor newHL, newHR;
      updateH(HLs[N + cursor - 2], HRs[N - cursor - 2], A, B, *Hptr, *Hptr, HLn, HRn, newHL, newHR);
      if(N + cursor - 1 < HLs.size())
        HLs[N + cursor - 1] = newHL;
      else
        HLs.push_back(newHL);
      if(N - cursor - 1 < HRs.size())
        HRs[N - cursor - 1] = newHR;
      else{
        HRs.push_back(newHR);
      }
    }
    dir *= -1;
    cursor += 2*dir;
    cnt++;
    cout<<"sweep "<< cnt<<": E0 = "<< E0 / (2*N)<<endl;
  }
  return E0;
}

size_t hidx(int N, Side side, int l){
  if(side == Left)
    return l;
  else
    return (2*N-1) - l - 1;
}

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

Matrix makeMPS(Side side, const UniTensor& A, const UniTensor& B, Matrix& L){
  UniTensor tA = A;
  UniTensor tB = B;
  if(side == Left){
    if(tB.bond(2).dim() != L.row())
      return Matrix();
    bondcat(tB, L, 2);
  }
  else if(side == Right){
    if(tA.bond(0).dim() != L.row())
      return Matrix();
    bondcat(tA, L, 0);
  }
  else
    bondcat(tB, L, 0);
  int labelA[] = {-1, -2, 1};
  int labelB[] = {1, -3, -4};
  int labelS[] = {-1, -2, -3, -4};
  tA.setLabel(labelA);
  tB.setLabel(labelB);
  UniTensor S = contract(tA, tB, true);
  S.permute(labelS, 0);
  return S.getBlock();
}

Matrix trialState(const UniTensor& A, const UniTensor& B, vector<Matrix>& Ls){
  UniTensor tA = A;
  UniTensor tB = B;
  int labelA[] = {1, -3, -4};
  int labelB[] = {-1, -2, 1};
  int labelS[] = {-1, -2, -3, -4};
  tA.setLabel(labelA);
  tB.setLabel(labelB);
  bondcat(tA, Ls[1], 2);
  bondcat(tB, Ls[1], 0);
  tA.permute(labelA, 1);
  bondrm(tA, Ls[0], 0);
  UniTensor S = contract(tB, tA, true);
  S.permute(labelS, 0);
  return S.getBlock();
}

UniTensor makeL(const UniTensor& tA){
  UniTensor A = tA;
  int lA[] = {0, 1};
  int lAT[] = {2, 0};
  UniTensor AT = A;
  AT.transpose();
  A.setLabel(lA);
  AT.setLabel(lAT);
  return A * AT;
}

UniTensor makeR(const UniTensor& tB){
  UniTensor B = tB;
  int lB[] = {1, 0};
  int lBT[] = {0, 2};
  UniTensor BT = B;
  BT.transpose();
  B.setLabel(lB);
  BT.setLabel(lBT);
  return B * BT;
}

double mpsNorm(vector<UniTensor>& As, vector<UniTensor>& Bs, Matrix& L, Network& normL, Network& normR){
  UniTensor nL = makeL(As[0]);
  UniTensor nR = makeR(Bs[0]);
  for(int l = 1; l < As.size(); l++){
    normL.putTensorT("L", nL);
    normL.putTensor("A", As[l]);
    normL.putTensorT("AT", As[l]);
    nL = normL.launch();

    normR.putTensor("R", nR);
    normR.putTensor("B", Bs[l]);
    normR.putTensorT("BT", Bs[l]);
    nR = normR.launch();
  }
  bondcat(nR, L, 0);
  bondcat(nR, L, 1);
  return (nL * nR)[0];
}

double mpsExp2s(vector<UniTensor>& As, vector<UniTensor>& Bs, Matrix& L, int s, UniTensor& Ob, Network& normL, Network& normR, Network& expOb){
  int N = As.size();
  UniTensor nL = makeL(As[0]);
  UniTensor nR = makeR(Bs[0]);
  for(int l = 1; l < N - 1; l++){
    normL.putTensorT("L", nL);
    normL.putTensor("A", As[l]);
    normL.putTensorT("AT", As[l]);
    nL = normL.launch();

    normR.putTensor("R", nR);
    normR.putTensor("B", Bs[l]);
    normR.putTensorT("BT", Bs[l]);
    nR = normR.launch();
  }
  UniTensor BN = Bs[N - 1];
  bondcat(BN, L, 0);
  expOb.putTensor("L", nL);
  expOb.putTensor("R", nR);
  expOb.putTensor("A", As[N - 1]);
  expOb.putTensor("B", BN);
  expOb.putTensorT("AT", As[N - 1]);
  expOb.putTensorT("BT", BN);
  expOb.putTensor("Ob", Ob);
  return expOb.launch()[0];
}

double mpsPrdS(vector<UniTensor>& As, vector<UniTensor>& Bs, Matrix& L, vector<UniTensor>prdState, Network& expPrdL, Network& expPrdR){
  int N = As.size();
  UniTensor leftS = As[0];
  UniTensor rightS = Bs[0];
  for(int l = 0; l < N - 1; l++){
    expPrdL.putTensor("L", leftS);
    expPrdL.putTensor("A", As[l + 1]);
    expPrdL.putTensor("S", prdState[l]);
    leftS = expPrdL.launch();

    expPrdR.putTensor("R", rightS);
    expPrdR.putTensor("B", Bs[l + 1]);
    expPrdR.putTensor("S", prdState[2*N - l - 1]);
    rightS = expPrdR.launch();
  }

  int lA2[] = {0, 1};
  int lB2[] = {1, 0};
  bondcat(rightS, L, 0);
  leftS.setLabel(lA2);
  rightS.setLabel(lB2);
  leftS *= prdState[N - 1];
  rightS *= prdState[N];
  return (leftS * rightS)[0];
}

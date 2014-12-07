enum Side{Left = -1, Right= 1};
UniTensor combineH(const UniTensor& H0, const UniTensor& HL, const UniTensor& HR);
UniTensor findGS(const UniTensor& SB, double& E0, Matrix& refState, int& iter);
int updateMPS(const UniTensor& GS, size_t chi, UniTensor& A, UniTensor& B);
double sweep(int N, int chi, int range, int times, UniTensor& H0, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn);
double sweep(int N, int chi, int range, int times, vector<UniTensor>& H0s, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn);
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
  iter= blk.lanczosEigh(E0, refState, 200);

  vector<Bond> stateBonds;
  for(int b = 4; b < 8; b++)
    stateBonds.push_back(SB.bond(b));
	UniTensor GS(stateBonds, "GS");
  GS.setElem(refState.getElem());
  return GS.permute(2);
}

int updateMPS(const UniTensor& GS, size_t chi, UniTensor& A, UniTensor& B){
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
  assert(range < N);
  Matrix psi;
  int dir = Left; //direction
  int cursor = -1;
  int cnt = 0;
  double E0;
  while(cnt < times + 1){
    UniTensor* Hptr;
    for( ; abs(cursor) < range; cursor += dir){
      if(H0s.size() > 1)
        Hptr = &(H0s[hidx(N, Left, N + cursor - 1)]);
      else
        Hptr = &(H0s[0]);
      UniTensor SB = combineH(*Hptr, HLs[N + cursor - 2], HRs[N - cursor - 2]);
      int iter;
      UniTensor GS = findGS(SB, E0, psi, iter);

      UniTensor A, B;
      int D = updateMPS(GS, chi, A, B);
      //cout<<"cursor = "<< cursor <<", D = " << chi << setprecision(10) << ", E = " << E0  << ", e = " << E0 / (2 * N) <<", iter = "<<iter<<endl;

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
  }
  return E0;
}

size_t hidx(int N, Side side, int l){
  if(side == Left)
    return l;
  else
    return (2*N-1) - l - 1;

}

UniTensor combineH(const UniTensor& H0, const UniTensor& HL, const UniTensor& HR);
UniTensor findGS(const UniTensor& SB, double& E0, Matrix& refState, int& iter);
int updateMPS(const UniTensor& GS, size_t chi, UniTensor& A, UniTensor& B);
void sweep(int N, int chi, int range, int times, UniTensor& H0, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn);

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

void sweep(int N, int chi, int range, int times, UniTensor& H0, vector<UniTensor>& HLs, vector<UniTensor>& HRs, Network& HLn, Network& HRn){
  assert(range < N);
  Matrix psi;
  int dir = -1; //direction
  int idx = -1;
  int cnt = 0;
  while(cnt < times + 1){
    for( ; abs(idx) < range; idx += dir){
      UniTensor SB = combineH(H0, HLs[N + idx - 2], HRs[N - idx - 2]);
      double E0;
      int iter;
      UniTensor GS = findGS(SB, E0, psi, iter);

      UniTensor A, B;
      int D = updateMPS(GS, chi, A, B);
      cout<<"idx = "<< idx <<", D = " << chi << setprecision(10) << ", E = " << E0  << ", e = " << E0 / (2 * N) <<", iter = "<<iter<<endl;

      UniTensor newHL, newHR;
      updateH(HLs[N + idx - 2], HRs[N - idx - 2], A, B, H0, H0, HLn, HRn, newHL, newHR);
      if(N + idx - 1 < HLs.size())
        HLs[N + idx - 1] = newHL;
      else
        HLs.push_back(newHL);
      if(N - idx - 1 < HRs.size())
        HRs[N - idx - 1] = newHR;
      else{
        HRs.push_back(newHR);
      }
    }
    dir *= -1;
    idx += 2*dir;
    cnt++;
  }
}

//Bondcat.rm
void bondcat(UniTensor& T, const Matrix& L, int bidx);
void bondrm(UniTensor& T, const Matrix& L, int bidx);

//Heisenberg haniltonian in N sites.
UniTensor HeisenbergOneHalf(int N);
int flip(int i, int A, int B, map<int, vector<double> > &state);
UniTensor Exp_dT(UniTensor &T, double delta);

//IniTensor T, Lambda
UniTensor initTensor(int d1, int D1, int gr);
map<string, UniTensor> initLambda(int D1, int gr);

// (UniTensor A, put vector of subbonds of U inbond in 2nd argument, PATH/of/.netfile/you/want)
map<string, UniTensor> HOSVD(UniTensor &A, vector<Bond> &sbd, string PATH);
map<string, UniTensor> truncatePESS(map<string, UniTensor> &A, int CHI);
UniTensor constructT(map<string, UniTensor> &hosvd, string PATH);

//Return total dimensins of in/outbonds of tensor A
int inBondDim(UniTensor &A);
int outBondDim(UniTensor &A);

//Retrun Elem of Tensor
Matrix getTensorElem(UniTensor &T);

//SimpleUpdatePESS
void simpleUpdatePESS(bool upT, map<string, UniTensor> &US, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH,Network& simplex);
//void simpleUpdatePESSDir1(map<string, UniTensor> &Tu, map<string, UniTensor> &Td, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH);
//void simpleUpdatePESSDir2(map<string, UniTensor> &Td, map<string, UniTensor> &Tu, map<string, UniTensor> &Lb, map<string, UniTensor> &La, UniTensor &expH, string PATH);
//void simpleUpdatePESSDir3(map<string, UniTensor> &Td, map<string, UniTensor> &Tu, map<string, UniTensor> &Lb, map<string, UniTensor> &La, UniTensor &expH, string PATH);
//void simpleUpdatePESSDir4(map<string, UniTensor> &Tu, map<string, UniTensor> &Td, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH);
//void simpleUpdatePESSDir5(map<string, UniTensor> &Td, map<string, UniTensor> &Tu, map<string, UniTensor> &Lb, map<string, UniTensor> &La, UniTensor &expH, string PATH);
//void simpleUpdatePESSDir6(map<string, UniTensor> &Td, map<string, UniTensor> &Tu, map<string, UniTensor> &Lb, map<string, UniTensor> &La, UniTensor &expH, string PATH);
//void simpleUpdatePESS(UniTensor &Ta, map<string, UniTensor> &La, UniTensor &expH, string PATH);
map<string, UniTensor> invLambda(map<string, UniTensor>& L);

//Measure
UniTensor TLaunch(map<string, UniTensor> &T);

//Wite netfile
string netFileDir(string basedir, string filename);
void writeTbarNetFile(UniTensor &T, UniTensor &expH, string Ts, string expHs, string PATH);
void writeNewTNetFile(map<string, UniTensor> &hosvd, string PATH);
void writeCatNetFile(map<string, UniTensor> &Lambda, string PATH);
void writeMeasureNetFile(map<string,UniTensor> &Measure,string Tname,string Obsname,string PATH);

//================================================================//

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
  for(int i = 0; i < invL.elemNum(); i++)
    invL[i] = invL[i] == 0 ? 0 : (1 / invL[i]);
	bondcat(T, invL, bidx);
}

UniTensor TLaunch(map<string, UniTensor> &groupT){

  Network contractTup("M.net");
  contractTup.putTensor("U0", &groupT["U0"]);
  contractTup.putTensor("U1", &groupT["U1"]);
  contractTup.putTensor("U2", &groupT["U2"]);
  contractTup.putTensor("S0", &groupT["S0"]);

  UniTensor T = contractTup.launch();

  return T;
}

//================================================================//

// 1/2 Heisenberg haniltonian in N sites.
UniTensor HeisenbergOneHalf(int N){
  assert( N >= 2 );
  Matrix ham(pow(2, N), pow(2, N));
  ham.set_zero();
  map<int, vector<double> > state;
  for(int i = 0; i < N; i++)
    for(int j = 0; j < pow(2, i+1); j++)
      for(int k = 0; k < pow(2, N-1-i); k++){
        if( j % 2 == 0)
          state[j*pow(2, N-1-i)+k].push_back(0.5);
        else
          state[j*pow(2, N-1-i)+k].push_back(-0.5);
      }
  int j, b;
  double c = (N == 2) ? 0.5 : 1;        //If N were equal to 2, it would calculate twice.
  for(int a = 0; a < pow(2, N); a++){
    for(int i = 0; i < N; i++){
      j = (i+1) % N;
      if(state[a][i] == state[a][j]){
        ham.at(a, a) = ham.at(a, a) + 0.25 * c;
      }
      else{
        ham.at(a, a) = ham.at(a, a) - 0.25 * c;
        b = flip(a, i, j, state);
        ham.at(a , b) = 0.5;
        }
      }
    }
  //InitTensorH;
  Qnum q0(0);
  Bond bdi_d(BD_IN, d);
  Bond bdo_d(BD_OUT, d);
  vector<Bond> bondH(2*N);
  for(int i = 0; i < N; i++){
    bondH.at(i) = bdi_d;
    bondH.at(i+N) = bdo_d;
  }
  UniTensor H(bondH);
  H.putBlock(q0, ham);
  return H;
}

int flip(int a, int A, int B, map<int, vector<double> > &state){
  int P;
  vector<double> tmp = state[a];
  tmp[A] = (tmp[A] == 0.5) ? -0.5 : 0.5;
  tmp[B] = (tmp[B] == 0.5) ? -0.5 : 0.5;
  for(int i = 0; i < state.size(); i++){
    if(tmp == state[i])
      P = i;
  }
  return P;
}

UniTensor Exp_dT(UniTensor &T, double delta){
	assert(T.bondNum() == GR*2);
	int M = inBondDim(T);
	int N = outBondDim(T);
	assert(M == N);

	UniTensor U, UT, expT;
	UT = T;
	double *Tmp = (double*)malloc(sizeof(double)* M*N);
	for(int i = 0; i < M*N; i++)
		Tmp[i] = getTensorElem(T).getElem()[i] * -delta;
	UT.setRawElem(Tmp);

	//Eigenvalue decomposition
	//T = U*E*UT
	int ldA = N;
	int lwork = 10*N;
	double *work = (double*)malloc(sizeof(double) * lwork);
	double *Eig  = (double*)malloc(sizeof(double) * N);
	int info;
	dsyev_((char*)"V", (char*)"U", &N, Tmp, &ldA, Eig, work, &lwork, &info);
	assert(info == 0);
	UT.setRawElem(Tmp);
	U = UT.transpose();
	//exp(Eig)
	for(int i = 0; i < N; i++)
		Eig[i] = exp(Eig[i]);
	for(int i = 0; i < M; i++)
		for(int j = 0; j < N; j++)
			Tmp[i*N + j] *= Eig[i];
	UT.setRawElem(Tmp);
	//Reset U and UT label
	vector<int> Ulabel(T.bondNum());
	vector<int> UTlabel(T.bondNum());
	for(int i = 0; i < T.bondNum(); i++){
	  Ulabel.at(i) = i;
	  UTlabel.at(i) = i+GR;
	}
	U.setLabel(Ulabel);
	UT.setLabel(UTlabel);

	expT = U * UT;

  vector<int> labelexpT;
  labelexpT = T.label();
	expT.setLabel(labelexpT);

	return expT;

  Ulabel.clear();
  UTlabel.clear();
	free(Eig);  free(work);
}


//IniTensor T
UniTensor initTensor(int d1, int D1, int gr){

  Qnum q0(0);
  Bond bdi_d(BD_IN, d1);
  Bond bdi_D(BD_IN, D1);
  vector<Bond> bondT(gr*2);
  for(int i = 0; i < gr; i++){
    bondT.at(2*i) = bdi_d;
    bondT.at(2*i+1) = bdi_D;
  }

  vector<UniTensor> initT;
  UniTensor T(bondT);
  T.randomize();
  double norm = T.getBlock(q0).norm();
  Matrix T_elem = T.getBlock(q0) * (1/norm);
  T.putBlock(q0, T_elem);

  return T;
}

map<string, UniTensor> initLambda(int D1, int gr){

	map<string, UniTensor> Lambda;

	Qnum q0(0);
	Bond bdi_D(BD_IN, D1);
	Bond bdo_D(BD_OUT, D1);
	vector<Bond> bondL;
	bondL.push_back(bdi_D);
	bondL.push_back(bdo_D);

	UniTensor L(bondL);
	Matrix L_elem(D1, D1, true);
	L_elem.randomize();

	double normL = L_elem.norm();
	L_elem = L_elem * (1/ normL);
  L.putBlock(q0, L_elem);

  char Lname[4];
  for(int i = 0; i < gr; i++){
    sprintf(Lname, "L%d", i);
    Lambda[Lname] = L;
  }
  return Lambda;
}

map<string, UniTensor> HOSVD(UniTensor &A, vector<Bond> &sbd, string PATH){
    //In Kagome Lattice odd labels represent virtual bonds and even labels represent physical bonds
    //Number of subbonds in each inbond of U
    int sbn = sbd.size();
    //Greographical Representation
    int gr = A.bondNum() / sbn;
    assert(A.bondNum() % sbn == 0);

    //clone A, A1
    UniTensor A1 = A;
    //Praparing to do svd. The fisrt step is reshaping Tensor A1 to 1st model
    A1.permute(sbn);
    int iBD = inBondDim(A1);
    // iBD is the bond dimention of inBond numbers of U
    Bond bdo_Ui(BD_IN, iBD);
    Bond bdo_Uo(BD_OUT, iBD);

    //Initialize U's bonds
    vector<Bond> bond2(sbn+1);
    for(int i = 0; i < sbn; i++)
        bond2.at(i) = sbd[i];
    bond2.at(sbn) = bdo_Uo;

    vector<Bond> bondL;
    bondL.push_back(bdo_Ui);
    bondL.push_back(bdo_Uo);

    //Tmp container include  S, U and Lambda tensors; gr Unitary
    //tensors and Lambdas
    vector<UniTensor> hosvd(2*gr);
    Qnum q0(0);
    Matrix MA;

    //Initialize labels of A1 and core tensor
    vector<int> initLabel;
    vector<int> initLabelS;
    for(int i = 0; i < gr*sbn; i++)
        initLabel.push_back(i);
    A1.setLabel(initLabel);

    for(int i = 0; i < gr; i++)
        initLabelS.push_back(i);

    //Check whether fiting the restriction of HOSVD
    //For eample, in PESS, Each inbond of U, include a virturl bond and a physcial bond; Uinbond = {d, CHI};
    for(int i = 0; i < sbn; i++){
        for(int j = 0; j < gr; j++)
            assert(A1.bond()[sbn*j+i].dim() == sbd[i].dim());
    }

    //Reshap Tensor A1 to different models and do svd
    UniTensor TmpU;
    UniTensor TmpL;
    int tmp;
    vector<int> pLable;
    //Geographical
    for(int i = 0; i < gr; i++){
        //Create New Labels
        for(int j = 0; j < sbn*gr; j++){
            tmp = (sbn*i+j) % (gr*sbn);
            pLable.push_back(tmp);
        }
        A1.permute(pLable, sbn);
        MA = A1.getBlock(q0);
        TmpU.assign(bond2);
        TmpU.putBlock(q0, MA.svd()[0]);
        TmpL.assign(bondL);
//        double norm = MA.svd()[1].norm();
        TmpL.putBlock(q0, MA.svd()[1]);
        hosvd.at(i) = TmpU;
        hosvd.at(i+gr) = TmpL;
        pLable.clear();
    }
    //recycle A1
    A1.permute(initLabel, sbn);

    //Write .netfile
    FILE *fp;
    if((fp =fopen(PATH.c_str(),"r")) == NULL){
    //Check PATH is wether a path of a .net file.
    string type_net = ".net";
    int Size = PATH.size();
    char filenameExtension[8];
    size_t length =  PATH.copy(filenameExtension, 4, Size-4);
    filenameExtension[length] = '\0';
    assert(filenameExtension == type_net);
    //Open .net file
    FILE* corenet = fopen(PATH.c_str(), "w");
    char TL[8];
    char TbondL[32];
    //Write Tout Labels
    for(int i = 0; i <= (gr+1); i++){
        if(i == (gr+1)){
            strcpy(TbondL, "TOUT:");
            for(int j = 0; j < gr; j++){
                int a = j + gr*sbn;
                sprintf(TL, "%d", a);
                strcat(TbondL, TL);
                if(j == 0){
                    strcat(TbondL, ";");
                }
                if(j == (gr-1))
                    fprintf(corenet, "%s\n", TbondL);
                if(j != 0  && j < (gr-1)){
                    strcat(TbondL, ",");
                }
            }
        }
        if(i == gr){
            strcpy(TbondL, "A:");
            for(int j = 0; j < gr*sbn; j++){
                sprintf(TL, "%d", j);
                strcat(TbondL, TL);
                if(j == (A1.inBondNum()-1))
                    strcat(TbondL, ";");
                if(j == (gr*sbn-1))
                    fprintf(corenet, "%s\n", TbondL);
                if(j != (A1.inBondNum()-1)  && j < (gr*sbn-1))
                    strcat(TbondL, ",");
            }
        }
        if(i < gr){
            sprintf(TbondL, "U%d:", i);
            for(int j = 0; j <= sbn; j++){
                int a = i*sbn+j;
                int b = i+sbn*gr;
                if(j == sbn-1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL,";");
                }
                if(j < sbn -1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL, ",");
                }
                if(j == sbn){
                    sprintf(TL, "%d", b);
                    strcat(TbondL, TL);
                    fprintf(corenet,"%s\n", TbondL);
                }
            }
        }
    }
    fclose(corenet);
  }else
    fclose(fp);
    //Get coretensor
    Network CoreTensor(PATH);
    char TUname[8];
    char TSname[8];
    char LDname[8];
    for(int i = 0; i <= gr; i++){
        if(i == gr)
            CoreTensor.putTensor("A", &A1);
        if(i < gr){
            sprintf(TUname, "U%d", i);
            CoreTensor.putTensor(TUname, &hosvd[i]);
            }
    }
    //core tensor of mode 0 --> S0
    UniTensor core = CoreTensor.launch();
    core.setLabel(initLabelS);
/**************** Reshape core tensor S ****************/
/*  vector<int> SLabel;
    for(int i = 0; i < gr; i++){
        for(int j = 0; j < gr*sbn; j++){
           tmp = (i*sbn+j)%(gr*sbn);
           SLabel.push_back(tmp);
        }
        core.permute(SLabel, sbn);
        hosvd.at(i+gr) = core;
        SLabel.clear();
    }                                                   */
/********************************************************/
    map<string, UniTensor> hosvdm;
    hosvdm["S0"] = core;
    for(int i = 0; i < gr; i++){
        sprintf(TUname, "U%d", i);
        sprintf(LDname, "L%d", i);
        hosvdm[TUname] = hosvd[i];
        hosvdm[LDname] = hosvd[i+gr];
    }
    return hosvdm;
}

map<string, UniTensor> truncatePESS(map<string, UniTensor> &A, int CHI){

    assert(A.size() % 2 == 1);
    assert(CHI <= inBondDim(A["U0"]));

    for(map<string, UniTensor>::iterator it = A.find("U0"); it != A.end(); ++it)
      assert( outBondDim(it->second) == inBondDim(A.find("S0")->second));

    //clone A["S0"] for keeping the values
    //clone A["S)"]
    UniTensor TmpSS = A["S0"];
    map<string, UniTensor> hosvdt;

    Qnum q0(0);
    int Sbn = A["S0"].bondNum(); //bond number of core tensor
    int Ubn = A["U0"].bondNum();
    vector<int> SLabel, S1Label;
    vector<int> initLabel = A["S0"].label();

    //bond with dimension chi
    Bond bdi_x(BD_IN, CHI);
    Bond bdo_x(BD_OUT, CHI);

    vector<Bond> bondS;
    vector<Bond> bondU = A["U0"].bond();
    UniTensor TmpS;
    UniTensor TmpU;

    char Tname[8];
    int Tmp;
    Matrix MU(CHI,inBondDim(A["U0"]));
    Matrix TmpM;

    for(int i = 0; i < Sbn; i++){

        //Truncate U
        sprintf(Tname, "U%d", i);
        TmpM = MU;
        TmpM.setElem(A[Tname].getBlock(q0).transpose().getElem());
        bondU.at(Ubn-1) = bdo_x;
        TmpU.assign(bondU);
        TmpU.putBlock(q0, TmpM.transpose());
        hosvdt[Tname] = TmpU;

        //Truncate S
        //Reshaping to different modes
        for(int j = 0; j < Sbn; j++){
            Tmp = (i+j) % Sbn;
            SLabel.push_back(Tmp);
        }
        A["S0"].permute(SLabel, 1);
        TmpS = A["S0"];
        S1Label = A["S0"].label();
        bondS = A["S0"].bond();
        bondS.at(0) = bdi_x;
        A["S0"].assign(bondS);
        Matrix MS(CHI, outBondDim(TmpS));
        TmpM = MS;
        TmpM.setElem(TmpS.getBlock(q0).getElem());
        A["S0"].putBlock(q0, TmpM);
        A["S0"].setLabel(S1Label);
        SLabel.clear();
    }
    hosvdt["S0"] = A["S0"].permute(initLabel, 1);
    A["S0"] = TmpSS;
    return hosvdt;
}

UniTensor constructT(map<string, UniTensor> &hosvd, string PATH){

  //Write .netfile
  FILE* fp;
  if((fp = fopen(PATH.c_str(),"r")) == NULL)
    writeNewTNetFile(hosvd, PATH);
  else
    fclose(fp);

  Network contractNewT(PATH);

  for(map<string, UniTensor>::iterator it = hosvd.find("S0"); it != hosvd.end(); ++it)
    contractNewT.putTensor(it->first, &it->second);

  UniTensor newT = contractNewT.launch();

  return newT;
}

// Get Elem of the tensor
Matrix getTensorElem(UniTensor &T){

	Qnum q0(0);
	vector<Qnum> qnums;
	qnums.push_back(q0);

	Bond bdi_q(BD_IN, qnums);
	Bond bdo_q(BD_OUT, qnums);
	vector<Bond> bondq;
	bondq.push_back(bdi_q);
	bondq.push_back(bdo_q);

	UniTensor Tmp(bondq);
	Tmp = T;

	Matrix Elem;
	Elem = Tmp.getBlock(q0);

	return Elem;
}

//Return total dimensins of in/outbonds of tensor A
int inBondDim(UniTensor &A){

    int iBD = 1;
    for(int i = 0; i < A.inBondNum(); i++)
        iBD *= A.bond()[i].dim();

    return iBD;
}

int outBondDim(UniTensor &A){

    int oBD = 1;
    for(int i = 0; i < A.bondNum()-A.inBondNum(); i++)
        oBD *= A.bond()[A.inBondNum()+i].dim();

    return oBD;
}

void simpleUpdatePESS(bool upT, map<string, UniTensor> &US, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH, Network& simplex){
  // upT = true for up triangle, otherwise down triangle

  for(map<string, UniTensor>::iterator ita = La.begin(); ita != La.end(); ++ita){
    string U = ita -> first;
    U[0] = 'U';
    map<string, UniTensor>::iterator itb = Lb.find(ita->first);
    if(upT)
      bondcat(US[U], itb->second.getBlock(), 2);
    else
      bondcat(US[U], ita->second.getBlock(), 1);
  }

  //UniTensor T = TLaunch(groupTu);
  simplex.putTensor("U0", US["U0"]);
  simplex.putTensor("U1", US["U1"]);
  simplex.putTensor("U2", US["U2"]);
  simplex.putTensor("S0", US["S0"]);
  simplex.putTensor("expH", expH);


  UniTensor theta = simplex.launch();
  vector<Bond> bondU;
  bondU.push_back(theta.bond(0));
  bondU.push_back(theta.bond(1));
  string corePATH = netFileDir(PATH, "Core_S0.net");
  map<string, UniTensor> hosvd = HOSVD(theta, bondU, corePATH);

  // Truncation
  if(upT){
    int per_label[] = {0, 2, 1};
    US["U0"].permute(per_label, 2);
    US["U1"].permute(per_label, 2);
    US["U2"].permute(per_label, 2);
  }

  US["U0"].putBlock(hosvd["U0"].getBlock().resize(d*D, D));
  US["U1"].putBlock(hosvd["U1"].getBlock().resize(d*D, D));
  US["U2"].putBlock(hosvd["U2"].getBlock().resize(d*D, D));

  if(upT)
    for(map<string, UniTensor>::iterator it = La.begin(); it != La.end(); ++it){
      map<string, UniTensor>::iterator itL = hosvd.find(it->first);
      Matrix mat = itL->second.getBlock().resize(D, D);
      double norm = mat.norm();
//      cout<<norm<<endl;
      it->second.putBlock(mat * (1.0/norm));
    }
  else
    for(map<string, UniTensor>::iterator it = Lb.begin(); it != Lb.end(); ++it){
      map<string, UniTensor>::iterator itL = hosvd.find(it->first);
      Matrix mat = itL->second.getBlock().resize(D, D);
//      cout<<itL->second;
      double norm = mat.norm();
      it->second.putBlock(mat * (1.0/norm));
//      cout<<norm<<endl;
    }
  UniTensor S0 = hosvd["S0"];
  int coreBondNum = S0.bondNum();
  Bond bDi = Bond(BD_IN, D);
  Bond bDo = Bond(BD_OUT, D);
  int initLabel[] = {0,1,2};
  for(int i = 0; i < coreBondNum; i++){
    vector<Bond> bonds = S0.bond();
    bonds[0] = bDi;
    UniTensor newS0(bonds);
    newS0.putBlock(S0.getBlock().resize(D, bonds[1].dim() * bonds[2].dim()));
    S0 = newS0;
    int per_label[3] = {1, 2, 0};
    S0.permute(per_label, 1);
  }
  S0.setLabel(initLabel);
  Matrix sv = S0.getBlock();
  S0.putBlock(sv * (1 / sv.norm()));
  US["S0"] = S0;

  if(upT){
    int per_label[] = {0, 1, 2};
    US["U0"].permute(per_label, 2);
    US["U1"].permute(per_label, 2);
    US["U2"].permute(per_label, 2);
  }

  // End of Truncation

  for(map<string, UniTensor>::iterator ita = La.begin(); ita != La.end(); ++ita){
    string U = ita -> first;
    U[0] = 'U';
    map<string, UniTensor>::iterator itb = Lb.find(ita->first);
    if(upT){
      bondcat(US[U], itb->second.getBlock(), 2);
    }
    else{
      bondcat(US[U], ita->second.getBlock(), 1);
    }
  }
}
/*
void simpleUpdatePESSDir1(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH){

  simpleUpdatePESS(groupTd, Lb, La, expH, PATH);

  UniTensor U2 = groupTd["U2"];
  map<string, UniTensor> groupTutmp = groupTu;
  groupTutmp["U2"] = U2;
  groupTu = groupTutmp;

  simpleUpdatePESS(groupTu, La, Lb, expH, PATH);

}

void simpleUpdatePESSDir2(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH){

  simpleUpdatePESS(groupTu, La, Lb, expH, PATH);

  UniTensor U2 = groupTu["U0"];
  map<string, UniTensor> groupTdtmp = groupTd;
  groupTdtmp["U0"] = U2;
  groupTd = groupTdtmp;

  simpleUpdatePESS(groupTd, Lb, La, expH, PATH);

}

void simpleUpdatePESSDir3(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH){

  simpleUpdatePESS(groupTd, Lb, La, expH, PATH);

  UniTensor U2 = groupTd["U1"];
  map<string, UniTensor> groupTutmp = groupTu;
  groupTutmp["U1"] = U2;
  groupTu = groupTutmp;

  simpleUpdatePESS(groupTu, La, Lb, expH, PATH);

}

void simpleUpdatePESSDir4(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH){

  simpleUpdatePESS(groupTu, La, Lb, expH, PATH);

  UniTensor U2 = groupTu["U2"];
  map<string, UniTensor> groupTdtmp = groupTd;
  groupTdtmp["U2"] = U2;
  groupTd = groupTdtmp;

  simpleUpdatePESS(groupTd, Lb, La, expH, PATH);

}

void simpleUpdatePESSDir5(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH){

  simpleUpdatePESS(groupTd, Lb, La, expH, PATH);

  UniTensor U2 = groupTd["U0"];
  map<string, UniTensor> groupTutmp = groupTu;
  groupTutmp["U0"] = U2;
  groupTu = groupTutmp;

  simpleUpdatePESS(groupTu, La, Lb, expH, PATH);

}

void simpleUpdatePESSDir6(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, map<string, UniTensor> &La, map<string, UniTensor> &Lb, UniTensor &expH, string PATH){

  simpleUpdatePESS(groupTu, La, Lb, expH, PATH);

  UniTensor U2 = groupTu["U1"];
  map<string, UniTensor> groupTdtmp = groupTd;
  groupTdtmp["U2"] = U2;
  groupTd = groupTdtmp;

  simpleUpdatePESS(groupTd, Lb, La, expH, PATH);

}
*/

map<string, UniTensor> invLambda(map<string, UniTensor>& L){

	Qnum q0(0);
	map<string, UniTensor> invL = L;

	for(map<string, UniTensor>::iterator it=L.begin(); it!=L.end(); ++it){
		assert(inBondDim(it->second) == outBondDim(it->second));
		Matrix mat(inBondDim(it->second),outBondDim(it->second), true);
		for(int i = 0; i < inBondDim(it->second); i++){
			if(it->second.getBlock(q0, true)[i] == 0)
				mat[i] = 0;
			else{
				mat[i] = 1 / it->second.getBlock(q0,true)[i];
			}
		}
		invL[it->first].putBlock(q0,mat);
	}
  return invL;
}

string netFileDir(string basedir, string filename){

  //Check basePATH
  char *cPATH = (char*)basedir.c_str();
  char lastNotation[4];
  int lengthb = basedir.size();
  strcpy(lastNotation, cPATH +(lengthb-1));
  string slash = "/";
  assert(lastNotation == slash);

  //chack Filename Extension
  string type_net = ".net";
  int Size = filename.size();
  char filenameExtension[8];
  size_t length =  filename.copy(filenameExtension, 4, Size-4);
  filenameExtension[length] = '\0';
  assert(filenameExtension == type_net);
  //Open .net file

  char dir[64];
  strcpy(dir, basedir.c_str());
  strcat(dir, filename.c_str());

  return dir;
}


void writeTbarNetFile(UniTensor &T, UniTensor &expH, string Ts, string expHs, string PATH){
    //Write .netfile
    int gr = T.bondNum() / 2;
    //Check PATH is wether a path of a .net file.
    string type_net = ".net";
    int Size = PATH.size();
    char filenameExtension[8];
    size_t length =  PATH.copy(filenameExtension, 4, Size-4);
    filenameExtension[length] = '\0';
    assert(filenameExtension == type_net);
    //Open .net file
    FILE* Tbarnet = fopen(PATH.c_str(), "w");
    char TL[8];
    char TbondL[32];
    //Write Tout Labels
    for(int i = 0; i < 3; i++){
      if(i == 2){
        strcpy(TbondL, "TOUT:");
        for(int j = 0; j < 2*gr; j++){
          sprintf(TL, "%d", j+1);
          strcat(TbondL, TL);
          if(j == (2*gr-1)){
            fprintf(Tbarnet, "%s\n", TbondL);
            break;
          }
          strcat(TbondL, ",");
        }
      }
      if(i == 1){
        strcpy(TbondL, Ts.c_str());
        strcat(TbondL, ":");
        for(int j = 0; j < 2*gr; j++){
          int a = (j%2 == 0) ? -(j+1) : (j+1);
          sprintf(TL, "%d", a);
          strcat(TbondL, TL);
          if(j == T.bondNum()-1)
            fprintf(Tbarnet, "%s\n", TbondL);
          if(j == T.inBondNum()-1)
            strcat(TbondL, ";");
          else
            strcat(TbondL, ",");
        }
      }
      if(i == 0){
        strcpy(TbondL, expHs.c_str());
        strcat(TbondL, ":");
        for(int j = 0; j < expH.inBondNum();j++){
          int a = 2*j + 1;
          sprintf(TL, "%d", a);
          strcat(TbondL, TL);
          if(j == expH.inBondNum()-1)
            break;
          strcat(TbondL, ",");
        }
        strcat(TbondL, ";");
        for(int j = 0; j < expH.bondNum()-expH.inBondNum(); j++){
          int a = 2*j + 1;
          sprintf(TL, "-%d", a);
          strcat(TbondL, TL);
          if(j == expH.bondNum()-expH.inBondNum()-1){
            fprintf(Tbarnet, "%s\n", TbondL);
            break;
          }
        strcat(TbondL, ",");
      }
    }
  }
  fclose(Tbarnet);
}

void writeNewTNetFile(map<string, UniTensor> &hosvd, string PATH){

    int gr = hosvd["S0"].bondNum();
    int sbn = hosvd["U0"].bondNum()-1;

    FILE* newTnet = fopen(PATH.c_str(), "w");

    char TL[8];
    char TbondL[32];
    //Write Tout Labels
    for(int i = 0; i <= (gr+1); i++){
        if(i == gr){
            strcpy(TbondL, "S0:");
            for(int j = 0; j < gr; j++){
                int a = j + gr*sbn;
                sprintf(TL, "%d", a);
                strcat(TbondL, TL);
                if(j == 0){
                    strcat(TbondL, ";");
                }
                if(j == (gr-1))
                    fprintf(newTnet, "%s\n", TbondL);
                if(j != 0  && j < (gr-1)){
                    strcat(TbondL, ",");
                }
            }
        }
        if(i == (gr+1)){
            strcpy(TbondL, "TOUT:");
            for(int j = 0; j < gr*sbn; j++){
                sprintf(TL, "%d", j);
                strcat(TbondL, TL);
                if(j == (gr*sbn-1))
                    fprintf(newTnet, "%s\n", TbondL);
                if(j < (gr*sbn-1))
                    strcat(TbondL, ",");
            }
        }
        if(i < gr){
            sprintf(TbondL, "U%d:", i);
            for(int j = 0; j <= sbn; j++){
                int a = i*sbn+j;
                int b = i+sbn*gr;
                if(j == sbn-1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL,";");
                }
                if(j < sbn -1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL, ",");
                }
                if(j == sbn){
                    sprintf(TL, "%d", b);
                    strcat(TbondL, TL);
                    fprintf(newTnet,"%s\n", TbondL);
                }
            }
        }
    }
    fclose(newTnet);
}

void writeMeasureNetFile(map<string,UniTensor> &Measure, string Tname, string Obsname, string PATH){

  FILE *netfile = fopen(PATH.c_str(), "w");
  char TL[8];
  char TbondL[32];
  for(int i = 0; i < 4; i++){
    if(i == 0){
      strcpy(TbondL, Tname.c_str());
      strcat(TbondL, ":");
      for(int j = 0; j < Measure[Obsname].bondNum() / 2; j++){
        if(j != (Measure[Obsname].bondNum() / 2) - 1 ){
          sprintf(TL, "%d," ,j+1);
          strcat(TbondL, TL);
          sprintf(TL, "-%d,",j+1);
          strcat(TbondL, TL);
        }else{
          sprintf(TL, "%d," ,j+1);
          strcat(TbondL, TL);
          sprintf(TL, "-%d", j+1);
          strcat(TbondL, TL);
          fprintf(netfile, "%s\n", TbondL);
        }
      }
    }
    if(i == 1){
      strcpy(TbondL, Tname.c_str());
      strcat(TbondL, "T:");
      for(int j = 0; j < Measure[Obsname].bondNum() / 2; j++){
        if(j != (Measure[Obsname].bondNum() / 2) - 1 ){
          sprintf(TL, "%lu," ,j+1+Measure[Obsname].bondNum()/2);
          strcat(TbondL, TL);
          sprintf(TL, "-%d,",j+1);
          strcat(TbondL, TL);
        }else{
          sprintf(TL, "%lu," ,j+1+Measure[Obsname].bondNum()/2);
          strcat(TbondL, TL);
          sprintf(TL, "-%d", j+1);
          strcat(TbondL, TL);
          fprintf(netfile, "%s\n", TbondL);
        }
      }
    }
    if(i == 2){
      strcpy(TbondL, Obsname.c_str());
      strcat(TbondL, ":");
      for(int j = 0; j < Measure[Obsname].bondNum(); j++){
        if( j == Measure[Obsname].inBondNum()-1){
          sprintf(TL, "%d;" ,j+1);
          strcat(TbondL, TL);
        }if( j == Measure[Obsname].bondNum()-1){
          sprintf(TL, "%d" ,j+1);
          strcat(TbondL, TL);
          fprintf(netfile, "%s\n", TbondL);
        }if (j != Measure[Obsname].bondNum()-1 && j != Measure[Obsname].inBondNum()-1){
          sprintf(TL, "%d," ,j+1);
          strcat(TbondL, TL);
        }
      }
    }if(i == 3){
      strcpy(TbondL, "TOUT:");
      fprintf(netfile, "%s\n", TbondL);
    }
  }
  fclose(netfile);
}

void writeCatNetFile(map<string, UniTensor> &Lambda, string PATH){

  FILE *netFile = fopen(PATH.c_str(), "w");
  char TL[16];
  char TbondL[32];

  int b = 1;
  for(map<string, UniTensor>::iterator it = Lambda.begin(); it != Lambda.end(); ++it){
    if(it->second.bondNum() == 2){
      strcpy(TbondL, it->first.c_str());
      strcat(TbondL, ":");
      sprintf(TL, "%d;-%d", b, b);
      strcat(TbondL, TL);
      fprintf(netFile, "%s\n", TbondL);
      b += 2;
    }else{
      for(int i = 0; i < 2; i++){
        if(i == 0){
          assert(it->second.inBondNum() == it->second.bondNum());
          strcpy(TbondL, it->first.c_str());
          strcat(TbondL, ":");
          for(int j = 0; j < it->second.bondNum(); j++){
            if(j % 2 == 0){
              sprintf(TL, "%d,", j);
              strcat(TbondL, TL);
            }
            if(j % 2 == 1 && j != it->second.bondNum()-1){
              sprintf(TL, "-%d,", j);
              strcat(TbondL, TL);
            }
            if(j == it->second.bondNum()-1){
              sprintf(TL, "-%d", j);
              strcat(TbondL, TL);
              fprintf(netFile, "%s\n", TbondL);
            }
          }
        }if(i == 1){
          strcpy(TbondL, "TOUT:");
          for(int j = 0; j < it->second.bondNum(); j++){
            if(j == it->second.bondNum()-1){
              sprintf(TL, "%d", j);
              strcat(TbondL, TL);
              fprintf(netFile, "%s\n", TbondL);
            }else{
              sprintf(TL, "%d,", j);
              strcat(TbondL, TL);
            }
          }
        }
      }
    }
  }
  fclose(netFile);
}

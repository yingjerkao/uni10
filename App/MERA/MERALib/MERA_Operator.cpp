//Ascend the Observables "Ob" to the higher layer observables in MERA. Here we consider translation invariant;
SyTensor_t Ascend(SyTensor_t& Ob, SyTensor_t& W, SyTensor_t& U, Network_t& asdL, Network_t& asdC, Network_t& asdR);
//Descend the Density matrix "Rho" to the lower layer density matrix in MERA. Here we consider translation invariant;
//void Descend(Tensor* Rho, Tensor* W, Tensor* U, Tensor* Tout);
//Here we assume that W of a layer are the same(that is W1=W2). Diagram list is a list of all the diagrams in MERA.
//void UpdateW(Tensor* W, Tensor* Rho, Tensor* Ob, Tensor* U);
//Here we assume that W of a layer are the same(that is W1=W2). Diagram list is a list of all the diagrams in MERA.
//void UpdateU(Tensor* U, Tensor* Rho, Tensor* Ob, Tensor* W);
//Diagonalize Observable and keep only the state with minimal energy to generate a pure desity matrix.
//void Filter(Tensor* Ob, Tensor* newRho, double* Ground_state_energy);

//For Ternary MERA


SyTensor_t Ascend(SyTensor_t& Ob, SyTensor_t& W, SyTensor_t& U, Network_t& asdL, Network_t& asdC, Network_t& asdR){	
	SyTensor_t W2 = W;
	SyTensor_t UT = U;
	UT.transpose();
	SyTensor_t WT = W;
	WT.transpose();
	SyTensor_t W2T = W;
	W2T.transpose();

	vector<SyTensor_t*> tens;
	tens.push_back(&W);
	tens.push_back(&W2);
	tens.push_back(&U);
	tens.push_back(&Ob);
	tens.push_back(&UT);
	tens.push_back(&WT);
	tens.push_back(&W2T);
	
	for(int i = 0; i < tens.size(); i++){
		asdL.replaceWith(i, tens[i], 1);
		asdC.replaceWith(i, tens[i], 1);
		asdR.replaceWith(i, tens[i], 1);
	}

	SyTensor_t newH = asdL.launch("Ob");
	newH += asdC.launch();
	newH += asdR.launch();
	newH *= 1.0 / 3;
	return newH;
}

SyTensor_t Descend(SyTensor_t& Rho, SyTensor_t& W, SyTensor_t& U, Network_t& desL, Network_t& desC, Network_t& desR){
	SyTensor_t W2 = W;
	SyTensor_t UT = U;
	UT.transpose();
	SyTensor_t WT = W;
	WT.transpose();
	SyTensor_t W2T = W;
	W2T.transpose();

	vector<SyTensor_t*> tens;
	tens.push_back(&UT);
	tens.push_back(&WT);
	tens.push_back(&W2T);
	tens.push_back(&Rho);
	tens.push_back(&W);
	tens.push_back(&W2);
	tens.push_back(&U);
	
	for(int i = 0; i < tens.size(); i++){
		desL.replaceWith(i, tens[i], 1);
		desC.replaceWith(i, tens[i], 1);
		desR.replaceWith(i, tens[i], 1);
	}

	SyTensor_t newR = desL.launch("Rho");
	newR += desC.launch();
	newR += desR.launch();
	newR *= 1.0 / 3;
	return newR;
}

SyTensor_t UpdateW(SyTensor_t& W, SyTensor_t& Rho, SyTensor_t& Ob, SyTensor_t& U, Network_t& W1EnvL, Network_t& W1EnvC, Network_t& W1EnvR, Network_t& W2EnvL, Network_t& W2EnvC, Network_t& W2EnvR){
	SyTensor_t W1 = W;
	SyTensor_t W2 = W;
	SyTensor_t UT = U;
	UT.transpose();
	SyTensor_t W1T = W;
	W1T.transpose();
	SyTensor_t W2T = W;
	W2T.transpose();

	vector<SyTensor_t*> tens1;
	tens1.push_back(&W2);
	tens1.push_back(&U);
	tens1.push_back(&Ob);
    tens1.push_back(&UT);
    tens1.push_back(&W1T);
	tens1.push_back(&W2T);
	tens1.push_back(&Rho);

	vector<SyTensor_t*> tens2;
	tens2.push_back(&W1);
	tens2.push_back(&U);
	tens2.push_back(&Ob);
    tens2.push_back(&UT);
    tens2.push_back(&W1T);
	tens2.push_back(&W2T);
	tens2.push_back(&Rho);

	for(int i = 0; i < tens1.size(); i++){
		W1EnvL.replaceWith(i, tens1[i]);
		W1EnvC.replaceWith(i, tens1[i]);
        W1EnvR.replaceWith(i, tens1[i]);
	}

	for(int i = 0; i < tens2.size(); i++){
		W2EnvL.replaceWith(i, tens2[i]);
		W2EnvC.replaceWith(i, tens2[i]);
        W2EnvR.replaceWith(i, tens2[i]);
	}

	SyTensor_t env = W1EnvL.launch();
    env += W1EnvC.launch();
	env += W1EnvR.launch();
	env += W2EnvL.launch();
	env += W2EnvC.launch();
	env += W2EnvR.launch();
	SyTensor_t newW = W;

	/*
    Matrix_t Env = printRawElem(env);
	cout<<Env;
	int Rnum = Env.row();
	int Cnum = Env.col();
	for(int i = 0; i < Cnum; i++){
		double tmp = Env.elem[3 * Cnum + i];
		Env.elem[3 * Cnum + i] = Env.elem[4 * Cnum + i];
		Env.elem[4 * Cnum + i] = tmp;
	}
	cout<<env;
	cout<<Env;

	vector<Matrix_t> ret_E = Env.svd();

	cout << ret_E[1];
	//cout << ret_E[0];
	//cout << env;
	*/
    /*
	for(int i = 0; i < qnums.size(); i++){
		Matrix_t mat = env.getBlock(qnums[i]);
		vector<Matrix_t> rets = mat.svd();
		rets[0].transpose();
		rets[2].transpose();
		Matrix_t blk = rets[2] * rets[0];
		blk *= -1;
		newW.putBlock(qnums[i], blk);
	}
	*/
	map<Qnum_t, Matrix_t> blocks = env.getBlocks();
	for(map<Qnum_t,Matrix_t>::iterator it = blocks.begin(); it != blocks.end(); it++){
		//cout<<it->second;
		vector<Matrix_t> rets = it->second.svd();
		//cout << it->second;
		//Matrix_t mat_n = rets[0]*rets[1]*rets[2];
		//cout << mat_n;
		rets[0].transpose();
		rets[2].transpose();
		Matrix_t blk = rets[2] * rets[0];
		blk *= -1;
		newW.putBlock(it->first, blk);
	}
	
	/*
	ret_E[0].transpose();
	ret_E[2].transpose();
	Matrix_t NW = ret_E[2]*ret_E[0];
	NW*=(-1);
	*/
	//printRawElem(V);
	//printRawElem(u);
	//return newW;
	return newW;
}

SyTensor_t UpdateU(SyTensor_t& U, SyTensor_t& Rho, SyTensor_t& Ob, SyTensor_t& W, Network_t& UEnvL, Network_t& UEnvC, Network_t& UEnvR){
	SyTensor_t W1 = W;
	SyTensor_t W2 = W;
	SyTensor_t UT = U;
	UT.transpose();
	SyTensor_t W1T = W;
	W1T.transpose();
	SyTensor_t W2T = W;
	W2T.transpose();

	vector<SyTensor_t*> tens;
	tens.push_back(&Ob);
	tens.push_back(&UT);
	tens.push_back(&W1T);
	tens.push_back(&W2T);
    tens.push_back(&Rho);
    tens.push_back(&W1);
	tens.push_back(&W2);

	for(int i = 0; i < tens.size(); i++){
		UEnvL.replaceWith(i, tens[i]);
		UEnvC.replaceWith(i, tens[i]);
        UEnvR.replaceWith(i, tens[i]);
	}

	SyTensor_t env = UEnvL.launch();
    env += UEnvC.launch();
	env += UEnvR.launch();
	SyTensor_t newU = U;

	map<Qnum_t, Matrix_t> blocks = env.getBlocks();
	for(map<Qnum_t,Matrix_t>::iterator it = blocks.begin(); it != blocks.end(); it++){
		vector<Matrix_t> rets = it->second.svd();
		rets[0].transpose();
		rets[2].transpose();
		Matrix_t blk = rets[2] * rets[0];
		blk *= -1;
		newU.putBlock(it->first, blk);
	}
	return newU;
}

SyTensor_t Filter(SyTensor_t& Ob, Qnum_t& Qnum){
	//cout << Qnum << endl;
	int RBondNum = Ob.getRBondNum();	// row bond number;
	int bondNum = Ob.getBondNum();	//bond number
	vector<int> labels(bondNum, 0);
	vector<int> rsp_labels(bondNum, 0);	//reshape order
	for(int b = 0; b < bondNum; b++){
		labels[b] = b;
		if(b < RBondNum)
			rsp_labels[b] = RBondNum - 1 - b;
		else
			rsp_labels[b] = bondNum - 1 - b + RBondNum;
	}
	SyTensor_t ObR = Ob;
	ObR.addLabel(labels);
	ObR.reshape(rsp_labels, RBondNum);
	SyTensor_t Her = (Ob + ObR)*(1.0/2);

	Matrix_t mat = Her.getBlock(Qnum);
	vector<Matrix_t> outs = mat.diagonalize();
	Matrix_t EigV = outs[1];
	vector<double> GS;
	int N = EigV.row();
	for(int i = 0; i < N; i++)
		GS.push_back(EigV.elem[i]);
	Matrix_t Dens = mat;
	double TrRho = 0.0;
	for(int i = 0; i < N; i++)
		TrRho += GS[i]*GS[i];
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++)
			Dens.elem[i*N+j] = (GS[i]*GS[j]) / TrRho;
	}
	/*
	for(int i = 0; i < GS.size(); i++)
		cout << i << "\t" << GS[i] << endl;
	*/
	SyTensor_t newRho = Her;
	newRho.bzero();
	newRho.putBlock(Qnum, Dens);
	return newRho;
}

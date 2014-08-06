//Ascend the Observables "Ob" to the higher layer observables in MERA. Here we consider translation invariant;
UniTensor Ascend(UniTensor& Ob, UniTensor& W, UniTensor& U, Network& asdL, Network& asdC, Network& asdR);
//Descend the Density matrix "Rho" to the lower layer density matrix in MERA. Here we consider translation invariant;
//void Descend(Tensor* Rho, Tensor* W, Tensor* U, Tensor* Tout);
//Here we assume that W of a layer are the same(that is W1=W2). Diagram list is a list of all the diagrams in MERA.
//void UpdateW(Tensor* W, Tensor* Rho, Tensor* Ob, Tensor* U);
//Here we assume that W of a layer are the same(that is W1=W2). Diagram list is a list of all the diagrams in MERA.
//void UpdateU(Tensor* U, Tensor* Rho, Tensor* Ob, Tensor* W);
//Diagonalize Observable and keep only the state with minimal energy to generate a pure desity matrix.
//void Filter(Tensor* Ob, Tensor* newRho, double* Ground_state_energy);

//For Ternary MERA


UniTensor Ascend(UniTensor& Ob, UniTensor& W, UniTensor& U, Network& asdL, Network& asdC, Network& asdR){	
	UniTensor W2 = W;
	UniTensor UT = U;
	UT.transpose();
	UniTensor WT = W;
	WT.transpose();
	UniTensor W2T = W;
	W2T.transpose();

	vector<UniTensor*> tens;
	tens.push_back(&W);
	tens.push_back(&W2);
	tens.push_back(&U);
	tens.push_back(&Ob);
	tens.push_back(&UT);
	tens.push_back(&WT);
	tens.push_back(&W2T);
	
	for(int i = 0; i < tens.size(); i++){
		asdL.putTensor(i, tens[i], 1);
		asdC.putTensor(i, tens[i], 1);
		asdR.putTensor(i, tens[i], 1);
	}
	
	UniTensor newH = asdL.launch("Ob");
	newH += asdC.launch();
	newH += asdR.launch();
	newH *= 1.0 / 3;
	return newH;
}

UniTensor Descend(UniTensor& Rho, UniTensor& W, UniTensor& U, Network& desL, Network& desC, Network& desR){
	UniTensor W2 = W;
	UniTensor UT = U;
	UT.transpose();
	UniTensor WT = W;
	WT.transpose();
	UniTensor W2T = W;
	W2T.transpose();

	vector<UniTensor*> tens;
	tens.push_back(&UT);
	tens.push_back(&WT);
	tens.push_back(&W2T);
	tens.push_back(&Rho);
	tens.push_back(&W);
	tens.push_back(&W2);
	tens.push_back(&U);
	
	for(int i = 0; i < tens.size(); i++){
		desL.putTensor(i, tens[i], 1);
		desC.putTensor(i, tens[i], 1);
		desR.putTensor(i, tens[i], 1);
	}

	UniTensor newR = desL.launch("Rho");
	newR += desC.launch();
	newR += desR.launch();
	newR *= 1.0 / 3;
	return newR;
}

UniTensor UpdateW(UniTensor& W, UniTensor& Rho, UniTensor& Ob, UniTensor& U, Network& W1EnvL, Network& W1EnvC, Network& W1EnvR, Network& W2EnvL, Network& W2EnvC, Network& W2EnvR){
	UniTensor W1 = W;
	UniTensor W2 = W;
	UniTensor UT = U;
	UT.transpose();
	UniTensor W1T = W;
	W1T.transpose();
	UniTensor W2T = W;
	W2T.transpose();

	vector<UniTensor*> tens1;
	tens1.push_back(&W2);
	tens1.push_back(&U);
	tens1.push_back(&Ob);
    tens1.push_back(&UT);
    tens1.push_back(&W1T);
	tens1.push_back(&W2T);
	tens1.push_back(&Rho);

	vector<UniTensor*> tens2;
	tens2.push_back(&W1);
	tens2.push_back(&U);
	tens2.push_back(&Ob);
    tens2.push_back(&UT);
    tens2.push_back(&W1T);
	tens2.push_back(&W2T);
	tens2.push_back(&Rho);

	for(int i = 0; i < tens1.size(); i++){
		W1EnvL.putTensor(i, tens1[i]);
		W1EnvC.putTensor(i, tens1[i]);
        W1EnvR.putTensor(i, tens1[i]);
	}

	for(int i = 0; i < tens2.size(); i++){
		W2EnvL.putTensor(i, tens2[i]);
		W2EnvC.putTensor(i, tens2[i]);
        W2EnvR.putTensor(i, tens2[i]);
	}

	UniTensor env = W1EnvL.launch();
    env += W1EnvC.launch();
	env += W1EnvR.launch();
	env += W2EnvL.launch();
	env += W2EnvC.launch();
	env += W2EnvR.launch();
	UniTensor newW = W;

	/*
    Matrix Env = printRawElem(env);
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

	vector<Matrix> ret_E = Env.svd();

	cout << ret_E[1];
	//cout << ret_E[0];
	//cout << env;
	*/
    /*
	for(int i = 0; i < qnums.size(); i++){
		Matrix mat = env.getBlock(qnums[i]);
		vector<Matrix> rets = mat.svd();
		rets[0].transpose();
		rets[2].transpose();
		Matrix blk = rets[2] * rets[0];
		blk *= -1;
		newW.putBlock(qnums[i], blk);
	}
	*/
	map<Qnum, Matrix> blocks = env.getBlocks();
	for(map<Qnum,Matrix>::iterator it = blocks.begin(); it != blocks.end(); it++){
		//cout<<it->second;
		vector<Matrix> rets = it->second.svd();
		//cout << it->second;
		//Matrix mat_n = rets[0]*rets[1]*rets[2];
		//cout << mat_n;
		rets[0].transpose();
		rets[2].transpose();
		Matrix blk = rets[2] * rets[0];
		blk *= -1;
		newW.putBlock(it->first, blk);
	}
	
	/*
	ret_E[0].transpose();
	ret_E[2].transpose();
	Matrix NW = ret_E[2]*ret_E[0];
	NW*=(-1);
	*/
	//printRawElem(V);
	//printRawElem(u);
	//return newW;
	return newW;
}

UniTensor UpdateU(UniTensor& U, UniTensor& Rho, UniTensor& Ob, UniTensor& W, Network& UEnvL, Network& UEnvC, Network& UEnvR){
	UniTensor W1 = W;
	UniTensor W2 = W;
	UniTensor UT = U;
	UT.transpose();
	UniTensor W1T = W;
	W1T.transpose();
	UniTensor W2T = W;
	W2T.transpose();

	vector<UniTensor*> tens;
	tens.push_back(&Ob);
	tens.push_back(&UT);
	tens.push_back(&W1T);
	tens.push_back(&W2T);
    tens.push_back(&Rho);
    tens.push_back(&W1);
	tens.push_back(&W2);

	for(int i = 0; i < tens.size(); i++){
		UEnvL.putTensor(i, tens[i]);
		UEnvC.putTensor(i, tens[i]);
        UEnvR.putTensor(i, tens[i]);
	}

	UniTensor env = UEnvL.launch();
    env += UEnvC.launch();
	env += UEnvR.launch();
	UniTensor newU = U;

	map<Qnum, Matrix> blocks = env.getBlocks();
	for(map<Qnum,Matrix>::iterator it = blocks.begin(); it != blocks.end(); it++){
		vector<Matrix> rets = it->second.svd();
		rets[0].transpose();
		rets[2].transpose();
		Matrix blk = rets[2] * rets[0];
		blk *= -1;
		newU.putBlock(it->first, blk);
	}
	return newU;
}

UniTensor Filter(UniTensor& Ob, Qnum& Qnum, int& Ei){
	//cout << Qnum << endl;
	int RBondNum = Ob.inBondNum();	// row bond number;
	int bondNum = Ob.inBondNum();	//bond number
	vector<int> labels(bondNum, 0);
	vector<int> per_labels(bondNum, 0);	//permute order
	for(int b = 0; b < bondNum; b++){
		labels[b] = b;
		if(b < RBondNum)
			per_labels[b] = RBondNum - 1 - b;
		else
			per_labels[b] = bondNum - 1 - b + RBondNum;
	}
	UniTensor ObR = Ob;
	ObR.addLabel(labels);
	ObR.permute(per_labels, RBondNum);
	UniTensor Her = (Ob + ObR)*(1.0/2);

	Matrix mat = Her.getBlock(Qnum);
	vector<Matrix> outs = mat.diagonalize();
	Matrix Eigv = outs[0];
	Matrix EigV = outs[1];
	Ei = Eigv[0];
	vector<double> GS;
	int N = EigV.row();
	for(int i = 0; i < N; i++)
		GS.push_back(EigV[i]);
	Matrix Dens = mat;
	double TrRho = 0.0;
	for(int i = 0; i < N; i++)
		TrRho += GS[i]*GS[i];
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++)
			Dens[i*N+j] = (GS[i]*GS[j]) / TrRho;
	}
	/*
	for(int i = 0; i < GS.size(); i++)
		cout << i << "\t" << GS[i] << endl;
	*/
	UniTensor newRho = Her;
	newRho.set_zero();
	newRho.putBlock(Qnum, Dens);
	return newRho;
}

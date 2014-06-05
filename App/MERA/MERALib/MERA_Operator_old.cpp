//Ascend the Observables "Ob" to the higher layer observables in MERA. Here we consider translation invariant;
void Ascend(Tensor* Ob, Tensor* W, Tensor* U, Tensor* Tout);
//Descend the Density matrix "Rho" to the lower layer density matrix in MERA. Here we consider translation invariant;
void Descend(Tensor* Rho, Tensor* W, Tensor* U, Tensor* Tout);
//Here we assume that W of a layer are the same(that is W1=W2). Diagram list is a list of all the diagrams in MERA.
void UpdateW(Tensor* W, Tensor* Rho, Tensor* Ob, Tensor* U);
//Here we assume that W of a layer are the same(that is W1=W2). Diagram list is a list of all the diagrams in MERA.
void UpdateU(Tensor* U, Tensor* Rho, Tensor* Ob, Tensor* W);
//Diagonalize Observable and keep only the state with minimal energy to generate a pure desity matrix.
void Filter(Tensor* Ob, Tensor* newRho, double* Ground_state_energy);

//For Ternary MERA
#define W_IN 1
#define W_OUT 3
#define U_IN 2
#define U_OUT 2
#define OBS_IN 2
#define OBS_OUT 2

void Ascend(Tensor* Ob, Tensor* W, Tensor* U, Tensor* Tout){
	//------This function depend on diagrams of MERA------
	Tensor* W2 = (Tensor*)calloc(1, sizeof(Tensor));
	clone(W, W2);
	Tensor* UT = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W1T = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W2T = (Tensor*)calloc(1, sizeof(Tensor));
	int WTorder[W_IN + W_OUT];
	for(int i = 0; i < W_IN + W_OUT; i++)
		WTorder[i] = (i + W_OUT) % (W_IN + W_OUT);
	int UTorder[U_IN + U_OUT];
	for(int i = 0; i < U_IN + U_OUT; i++)
		UTorder[i] = (i + U_OUT) % (U_IN + U_OUT);
	reshapeClone(U, UTorder, UT);
	reshapeClone(W, WTorder, W1T);
	clone(W1T, W2T);
	int n = Tout->elemNum, inc=1;
	double a = 1.0/3;
	Tensor* Atemp;
	Atemp = (Tensor*)calloc(1, sizeof(Tensor));
	memset(Tout->elem, 0, n * sizeof(double));
	AscendL(W, W1T, U, Ob, UT, W2, W2T, Atemp);
	daxpy(&n, &a, Atemp->elem, &inc, Tout->elem, &inc);
	AscendC(W, W1T, U, Ob, UT, W2, W2T, Atemp);
	daxpy(&n, &a, Atemp->elem, &inc, Tout->elem, &inc);
	AscendR(W, W1T, U, Ob, UT, W2, W2T, Atemp);
	daxpy(&n, &a, Atemp->elem, &inc, Tout->elem, &inc);
	Tout->status |= HAVEELEM;
	recycle(Atemp);
	recycle(W2);
	recycle(UT);
	recycle(W1T);
	recycle(W2T);
}

void Descend(Tensor* Rho, Tensor* W, Tensor* U, Tensor* Tout){
	//------This function depend on diagrams of MERA------
	Tensor* W2 = (Tensor*)calloc(1, sizeof(Tensor));
	clone(W, W2);
	Tensor* UT = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W1T = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W2T = (Tensor*)calloc(1, sizeof(Tensor));
	int WTorder[W_IN + W_OUT];
	for(int i = 0; i < W_IN + W_OUT; i++)
		WTorder[i] = (i + W_OUT) % (W_IN + W_OUT);
	int UTorder[U_IN + U_OUT];
	for(int i = 0; i < U_IN + U_OUT; i++)
		UTorder[i] = (i + U_OUT) % (U_IN + U_OUT);
	reshapeClone(U, UTorder, UT);
	reshapeClone(W, WTorder, W1T);
	clone(W1T, W2T);
	int n = Tout->elemNum, inc = 1;
	double a = 1.0/3;
	Tensor* Dtemp;
  	Dtemp = (Tensor*)calloc(1, sizeof(Tensor));
	memset(Tout->elem, 0, n * sizeof(double));
	DescendL(W1T, W, UT, Rho, U, W2T, W2, Dtemp);
	daxpy(&n, &a, Dtemp->elem, &inc, Tout->elem, &inc);
	DescendC(W1T, W, UT, Rho, U, W2T, W2, Dtemp);
	daxpy(&n, &a, Dtemp->elem, &inc, Tout->elem, &inc);
	DescendR(W1T, W, UT, Rho, U, W2T, W2, Dtemp);
	daxpy(&n, &a, Dtemp->elem, &inc, Tout->elem, &inc);	
	Tout->status |= HAVEELEM;
	recycle(Dtemp);
	recycle(W2);
	recycle(UT);
	recycle(W1T);
	recycle(W2T);
}

void UpdateW(Tensor* W, Tensor* Rho, Tensor* Ob, Tensor* U){
	//------This function depend on diagrams of MERA------
	Tensor* W1 = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W2 = (Tensor*)calloc(1, sizeof(Tensor));
	clone(W, W1);
	clone(W, W2);
	Tensor* UT = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W1T = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W2T = (Tensor*)calloc(1, sizeof(Tensor));
	int WTorder[W_IN + W_OUT];
	for(int i = 0; i < W_IN + W_OUT; i++)
		WTorder[i] = (i + W_OUT) % (W_IN + W_OUT);
	int UTorder[U_IN + U_OUT];
	for(int i = 0; i < U_IN + U_OUT; i++)
		UTorder[i] = (i + U_OUT) % (U_IN + U_OUT);
	reshapeClone(U, UTorder, UT);
	reshapeClone(W, WTorder, W1T);
	clone(W1T, W2T);
	int n = W->elemNum, inc = 1;
	double a = 1;
	double *gamma = (double*)malloc(n*sizeof(double));
	Tensor* Wtemp = (Tensor*)calloc(1, sizeof(Tensor));
	
	W1EnvL(Rho, W1T, U, Ob, UT, W2, W2T, Wtemp);
	dcopy(&n, Wtemp->elem, &inc, gamma, &inc);
	W1EnvC(Rho, W1T, U, Ob, UT, W2, W2T, Wtemp);
	daxpy(&n, &a, Wtemp->elem, &inc, gamma, &inc);
	W1EnvR(Rho, W1T, U, Ob, UT, W2, W2T, Wtemp);
	daxpy(&n, &a, Wtemp->elem, &inc, gamma, &inc);
	
	W2EnvL(Rho, W1, W1T, U, Ob, UT, W2T, Wtemp);
	daxpy(&n, &a, Wtemp->elem, &inc, gamma, &inc);
	W2EnvC(Rho, W1, W1T, U, Ob, UT, W2T, Wtemp);
	daxpy(&n, &a, Wtemp->elem, &inc, gamma, &inc);
	W2EnvR(Rho, W1, W1T, U, Ob, UT, W2T, Wtemp);
	daxpy(&n, &a, Wtemp->elem, &inc, gamma, &inc);
	
	int M = 1, N = 1;	//# of rows of updated W
	for(int i = 0; i<W_IN; i++)
		M *= W->bondDim[i];
	for(int i = W_IN; i<W_IN+W_OUT; i++)
		N *= W->bondDim[i];	//# of columns of updated W
	//gamma is a tensor of kind WT, thus gamma is a NxM matrix
	//Now we want to do SVD on gammaT(a MxN matrix), since M<=N
	//It is easy for us to treat gamma as a MxN matrix, since we want to do SVD in fortran arry format.
	//What we need to do is nothing but call gamma is a MxN matrix.(without transposition)
	assert(M <= N);
	int min = M;	//min = min(M,N)
	int ldA = M, ldu = M, ldvT = min;
	double *S = (double*)malloc(min*sizeof(double));
	double *vT = (double*)malloc(ldvT*N*sizeof(double));
	double *u = (double*)malloc(ldu*M*sizeof(double));
	int lwork = 12*N;
	double *work = (double*)malloc(lwork*sizeof(double));
	int info;
	//gammaT = u*S*vT
	dgesvd((char*)"A", (char*)"S", &M, &N, gamma, &ldA, S, u, &ldu, vT, &ldvT, work, &lwork, &info);
	free(gamma);
	assert(info == 0);
	double alpha = -1;
	double beta = 0;
	//WT = v*uT, but WT in fortran alignment is W in C.
	dgemm((char*)"T", (char*)"T", &N, &M, &M, &alpha, vT, &ldvT, u, &ldu, &beta, W->elem, &N);
	free(S), free(vT), free(u), free(work);
	recycle(Wtemp);
	recycle(UT);
	recycle(W1);
	recycle(W2);
	recycle(W1T);
	recycle(W2T);
}

void UpdateU(Tensor* U, Tensor* Rho, Tensor* Ob, Tensor* W){
	//------This function depend on diagrams of MERA------
	Tensor* W2 = (Tensor*)calloc(1, sizeof(Tensor));
	clone(W, W2);
	Tensor* UT = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W1T = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W2T = (Tensor*)calloc(1, sizeof(Tensor));
	int WTorder[W_IN + W_OUT];
	for(int i = 0; i < W_IN + W_OUT; i++)
		WTorder[i] = (i + W_OUT) % (W_IN + W_OUT);
	int UTorder[U_IN + U_OUT];
	for(int i = 0; i < U_IN + U_OUT; i++)
		UTorder[i] = (i + U_OUT) % (U_IN + U_OUT);
	reshapeClone(U, UTorder, UT);
	reshapeClone(W, WTorder, W1T);
	clone(W1T, W2T);
	int n = U->elemNum, inc = 1;
	double a = 1;
	double *gamma = (double*)malloc(n*sizeof(double));
	Tensor* Utemp = (Tensor*)calloc(1, sizeof(Tensor));
	UEnvL(Rho, W, W1T, Ob, UT, W2, W2T, Utemp);
	dcopy(&n, Utemp->elem, &inc, gamma, &inc);
	UEnvC(Rho, W, W1T, Ob, UT, W2, W2T, Utemp);
	daxpy(&n, &a, Utemp->elem, &inc, gamma, &inc);
	UEnvR(Rho, W, W1T, Ob, UT, W2, W2T, Utemp);
	daxpy(&n, &a, Utemp->elem, &inc, gamma, &inc);
	int M = 1, N = 1;	//# of rows of updated W
	for(int i = 0; i<U_IN; i++)
		M *= U->bondDim[i];
	for(int i = U_IN; i<U_IN+U_OUT; i++)
		N *= U->bondDim[i];	//# of columns of updated W
	//gamma is a tensor of kind WT, thus gamma is a NxM matrix
	//Now we want to do SVD on gammaT(a MxN matrix), since M<=N
	//It is easy for us to treat gamma as a MxN matrix, since we want to do SVD in fortran arry alignment.
	//What we need to do is nothing but call gamma is a MxN matrix.(without transposition)
	assert(M == N);
	int min = M; //min = min(M,N)
	int ldA = M, ldu = M, ldvT = min;
	double *S = (double*)malloc(min*sizeof(double));
	double *vT = (double*)malloc(ldvT*N*sizeof(double));
	double *u = (double*)malloc(ldu*M*sizeof(double));
	int lwork = 10*N;
	double *work = (double*)malloc(lwork*sizeof(double));
	int info;
	//gammaT = u*S*vT
	dgesvd((char*)"A", (char*)"S", &M, &N, gamma, &ldA, S, u, &ldu, vT, &ldvT, work, &lwork, &info);
	free(gamma);
	assert(info ==0);
	double alpha = -1;
	double beta = 0;
	//UT = v*uT, but UT in fortran alignment is U in C.
	dgemm((char*)"T", (char*)"T", &N, &M, &M, &alpha, vT, &ldvT, u, &ldu, &beta, U->elem, &N);
	free(S), free(vT), free(u), free(work);
	recycle(Utemp);
	recycle(W2);
	recycle(UT);
	recycle(W1T);
	recycle(W2T);
}

void Filter(Tensor* Ob, Tensor* Rho, double* GE){
	//------This function depend on diagrams of MERA------
	//------This function depend on Systems         ------
	int N = Ob->bondDim[0]*Ob->bondDim[1];
	assert(N == (Ob->bondDim[2]*Ob->bondDim[3]));
	double *A = (double*)malloc(sizeof(double)*N*N);
	double *elem = (double*)malloc(sizeof(double)*N*N);
	int i, j;
	//Construct a Hermition Operator A from Ob
	assert(OBS_IN == 2);
	int order[] = {1, 0, 3, 2};
	reshapeElem(Ob, order, elem);
	for(i = 0; i<N; i++)
		for(j = 0; j<N; j++)
			A[j*N+i] = (elem[i*N+j] + Ob->elem[i*N+j])/2;	//Symmetrize
	/*
	for(i = 0; i < N; i++)
		for(j = i; j < N; j++)
			assert(Ob->elem[i * N + j] == Ob->elem[j * N + i]);
	*/
	free(elem);
	int ldA = N;
	int lwork = 10*N;
	double* work= (double*)malloc(sizeof(double)*lwork);
	double* Eig = (double*)malloc(sizeof(double*)*N);	//Eigenvalue array, in ascending order
	int info;
	dsyev((char*)"V", (char*)"U", &N, A, &ldA, Eig, work, &lwork, &info);
	assert(info == 0);
	*GE = Eig[0];
	//ground state is the first eigenvector of A(The first column of A, A is in Fortran alignment)
	double *GS = A;	//Ground state, 1*N vector
	
	//Calculate the matrix of the pure ground state density matrix.
	//That is: Rho = |GS><GS|, in matrix representation <m|GS><GS|n>
	double TrRho = 0;
	for(i = 0; i<N; i++)
		TrRho += GS[i]*GS[i];
	for(i = 0; i<N; i++)
		for(j = 0; j<N; j++)
			Rho->elem[i*N+j] = GS[i]*GS[j]/TrRho;	//for TrRho = 1;
	Rho->status |= HAVEELEM;
	
	free(A);
	free(work), free(Eig);
}

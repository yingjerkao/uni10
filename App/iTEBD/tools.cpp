void printDMatrix(double *mt, int M, int N){	//OK
	printf("\n");
	for(int i = 0; i < M; i++){
		for(int j = 0; j < N; j++)
			printf("%8.3f", mt[i * N + j]);
		printf("\n");
	}
	printf("\n");
}

void diagonalize(double* Kij, int N, double* Eig, double* EigVec){	//OK
	memcpy(EigVec, Kij, N * N * sizeof(double));
	int ldA = N;
	int lwork = 4*N;
	double* work= (double*)malloc(sizeof(double)*lwork);
	int info;
	dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, &info);
	assert(info == 0);
	free(work);
}

void takeExp(double* Mij, int N, double alpha, double* expMij){
	double *M_eig = (double*)malloc(N * sizeof(double));
	double *M_eigV = (double*)malloc(N * N * sizeof(double));
	double *M_eigVT = (double*)malloc(N * N * sizeof(double));
	diagonalize(Mij, N, M_eig, M_eigVT);
	for(int i = 0; i < N; i++)
		M_eig[i] = exp(alpha * M_eig[i]);
	//transpose and multiplied by eigen values
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			M_eigV[i * N + j] = M_eigVT[j * N + i] * M_eig[j];
	myDgemm(M_eigV, M_eigVT, N, N, N, expMij);
	free(M_eig);
	free(M_eigV);
	free(M_eigVT);
}

void myDgesvd(double* Mij_ori, int M, int N, double* U, double* S, double* vT){ //not tested yet
	//Mij = U * S * VT
	double* Mij = (double*)malloc(M * N * sizeof(double));
	memcpy(Mij, Mij_ori, M * N * sizeof(double));
	int min = M < N ? M : N;	//min = min(M,N)
	int ldA = N, ldu = N, ldvT = min;
	int lwork = 12*N;
	double *work = (double*)malloc(lwork*sizeof(double));
	int info;
	//gammaT = u*S*vT
	dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, &info);
	assert(info == 0);
	free(work);
	free(Mij);
}


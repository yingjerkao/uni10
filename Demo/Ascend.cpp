#include <stdio.h>
#include <stdlib.h>

#include "network.h"

int main(int argc, char* argv[]){
	Tensor* W1 = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W1T = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W2 = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* W2T = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* U = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* UT = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* H = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* Tout = (Tensor*)calloc(1, sizeof(Tensor));
	int dim;
	if(argc < 2){
		printf("Give a chi!!\n");
		exit(0);
	}
	sscanf(argv[1], "%d", &dim);
	int bondNum = 4;
	int bondDim[4];
	for(int i = 0; i < bondNum; i++)
		bondDim[i] = dim;
	initTensor(H, bondNum, bondDim);
	initTensor(U, bondNum, bondDim);
	int W_bondDim[] = {dim, dim, dim, dim};
	initTensor(W1, bondNum, W_bondDim);
	initTensor(W2, bondNum, W_bondDim);
	OrthoRandomize(H, 2);
	//print_tensor(H, 2);
	OrthoRandomize(U, 2);
	OrthoRandomize(W1, 1);
	OrthoRandomize(W2, 1);
	int UT_rsp[] = {2, 3, 0, 1};
	int WT_rsp[] = {3, 0, 1, 2};
	reshapeClone(U, UT_rsp, UT);
	reshapeClone(W1, WT_rsp, W1T);
	reshapeClone(W2, WT_rsp, W2T);
	initNetwork((char*)"diasNet");
	time_t start, end;
	start = clock();		
	AscendL(W1, W1T, U, H, UT, W2, W2T, Tout);
	end = clock();
	printf("time = %f(sec)\n", float(end - start) / CLOCKS_PER_SEC);
	char outfn[] = "cTout";
	outTensorf(Tout, outfn);
	//print_tensor(Tout, 2);
}

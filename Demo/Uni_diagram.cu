#include <stdio.h>
#include <stdlib.h>
#define GPU 0
#include "network.h"
int main(){
	Tensor* H = (Tensor*)calloc(1, sizeof(Tensor));
	int bondNum = 4;
	int bondDim[] = {2, 2, 2, 2};
	initTensor(H, bondNum, bondDim);
	double H_elem[] = {1,  0,  0, 0,
					   0, -1,  0, 0,
					   0,  0, -1, 0,
					   0,  0,  0, 1};
	addElem(H, H_elem, HtoD);

	Tensor* U = (Tensor*)calloc(1, sizeof(Tensor));
	//int bondNum = 4;
	//int bondDim[] = {2, 2, 2, 2};
	initTensor(U, bondNum, bondDim);
	double U_elem[] = {0.5,  0.5,  0.5,  0.5,
					   0.5, -0.5,  0.5, -0.5,
					   0.5,  0.5, -0.5, -0.5,
					   0.5, -0.5, -0.5,  0.5};
	addElem(U, U_elem, HtoD);
	
	Tensor* UD = (Tensor*)calloc(1, sizeof(Tensor));
	int UD_rsp[] = {2, 3, 0, 1};
	reshapeClone(U, UD_rsp, UD);
	
	initNetwork((char*)"diasNet");
	Tensor* Tout = (Tensor*)calloc(1, sizeof(Tensor));
	uni_trans(U, UD, H, Tout);

	/*
	int UD_label[] = {1, 2, 5, 6};
	int H_label[] = {5, 6, 7, 8};
	int U_label[] = {7, 8, 3, 4};
	addLabel(UD, UD_label);
	addLabel(H, H_label);
	addLabel(U, U_label);
	Tensor* Ttmp = (Tensor*)calloc(1, sizeof(Tensor));
	contraction(H, U, Ttmp);
	Tensor* Tout = (Tensor*)calloc(1, sizeof(Tensor));
	contraction(UD, Ttmp, Tout);
	*/

	print_tensor(Tout, 2);
}

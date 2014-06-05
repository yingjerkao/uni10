#include <time.h>
#include <unistd.h> //the name of the header file that provides access to the POSIX operating system API.
#include <sys/stat.h>
#include <sys/types.h>
#include "network.h"

#define d 2 //the Hilbert space dimension of a particle
#define D 2 //the dimension of PEPS
#define CHI  20//the dimension of environment bonds (the dimension kept when contraction environment PEPS)

#include "Operator.cpp"
#include "Measurement.cpp"

#define DATA_NUM 	501 	//number of different transverse field values
#define MAX_ITER	1000 	//maximum number of iteration
#define DELTA   	(1.0E-10)	//precision of ground state, for contraction of environment E
#define Q_ONE		1		//Times of contracting LU, RD Environment before doing line contraction and measure E

int main(){
	//declare output file name.
	char* basedir = (char*)"../Data/";
	char  dirbuf[32];
	sprintf(dirbuf, "D%dCHI%d", D, CHI);
	char  dir[64];
	strcpy(dir, basedir);
	strcat(dir, dirbuf);
 	char *ret_fn = (char*)"result";
	char *out_fn = (char*)"out";
	char ret_path[128], out_path[128];
    strcpy(ret_path, dir); strcat(ret_path, ret_fn);
    strcpy(out_path, dir); strcat(out_path, out_fn);
	mkdir(dir,0755);

	//Giving lambda
	double lamb[DATA_NUM];
	double lambda;
	for (int dp = 0; dp < DATA_NUM; dp++)
		lamb[dp] = 0.5*dp;
	
	Tensor* A = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* B = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* a = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* b = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* H = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* E = (Tensor*)calloc(1, sizeof(Tensor));
	Tensor* C[4];
	Tensor* Ta[4];
	Tensor* Tb[4];
	double* LambdaU = (double*)malloc(D * sizeof(double));
	double* LambdaR = (double*)malloc(D * sizeof(double));
	double* LambdaD = (double*)malloc(D * sizeof(double));
    double* LambdaL = (double*)malloc(D * sizeof(double));
	initTensors(C, Ta, Tb, H, D*D);
	
	clock_t start, step. end;
	start = clock();

	FILE* ret = fopen(ret_path,"w");
	FILE* out = fopen(out_path,"w");

	for(int dp = 0; dp < DATA_NUM; dp++){
	lambda =lambda[dp];
	double g = lambda/4;
	fprintf(out,"\nlambda = %4.2\n", lambda);
	printf("lambda = %4.2f\n", lambda);
	double H_elem[] = { -1, g, g, 0,
						 g, 1, 0, g,
						 g, 0, 1, g,
						 0, g, g, -1 };
	addElem(H,H_elem);

	//get tensors from 2DiTEBD (simple update)
	
	inTensors(A, B, H, LambdaU, LambdaR, LambdaD, LambdaL, Lambda);
	AbsorbLambdas(A, B, LambdaU, LambdaR, LambdaD, LambdaL);
 	ContractPhys(A,a);
	ContractPhys(B,b);

	if(a->bondDim[0] != Ta[0]->bondDim[2])
		initTensors(C, Ta, Tb, H, a->bondDim[0]);
	//contract the environment
	double delta;
	double Ei, Ef;
	bool flag = true;
	int iter = 0;
	Ei = 0;
	do{
		iter++
	}	
	}
}
/*int main(){
	char *basedir = (char*)"../Data/CPU_TransIsing2/"
}*/


/*#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;


int main(){
	Qnum q0(0);
	vector<Qnum> qnums;
	qnums.push_back(q0); 
	qnums.push_back(q0);
	Bond bdi(BD_IN, qnums);
	Bond bdo(BD_OUT, qnums);
	vector<Bond> bonds;
	bonds.push_back(bdi);
	bonds.push_back(bdi);
	bonds.push_back(bdo);
	bonds.push_back(bdo);
	
	UniTensor H0(bonds);
	int l0145[] = {0, 1, 4, 5};
	H0.addLabel(l0145);

	double g = 0.5;
	double H_elem[] = {    1,      g,      g,     0,
						   g,     -1,      0,     g,
						   g,      0,     -1,     g,
						   0,      g,      g,     1};

	H0.addRawElem(H_elem);
	cout<<H0;


	UniTensor I(bonds);
	I.eye();
	int l2367[] = {2, 3, 6, 7};
	I.addLabel(l2367);
	cout<<I;

	UniTensor H1 = H0 * I;	
	int l01234567[] = {0, 1, 2, 3, 4, 5, 6, 7};
	H1.permute(l01234567, 4);
	cout<<H1;

	UniTensor H2 = H1;
	int l23016745[] = {2, 3, 0, 1, 6, 7, 4, 5};
	H2.permute(l23016745, 4);
	cout<<H2;

	UniTensor H3 = H1;
	int l02134657[] = {0, 2, 1, 3, 4, 6, 5, 7};
	H3.permute(l02134657, 4);
	cout<<H3;

	UniTensor H4 = H1;
	int l20316475[] = {2, 0, 3, 1, 6, 4, 7, 5};
	H4.permute(l20316475, 4);
	cout<<H4;

	UniTensor Hf = H1 + H2 + H3 + H4;
	cout<<Hf;

}*/


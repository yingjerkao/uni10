#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "network.h"
#include "Operator.cpp"
//System dependent variables
#define OBS_BOND_NUM 4		//the bond number of observable, here observable is H.
#define SITE_DIM	 2		//the Hilbert space dimension of a particle
#define SITE_NUM	 54
#define LAYER_NUM	 (int)ceil(log(SITE_NUM)/log(W_OUT))
#define CHI 		 12		//maximum bond dimension
//Iteration variables
#define Q_ONE		 10		//times for update U and W in one iteration
#define DELTA	  (1.0E-8)	//precision of ground state energy, used in doing iteration
#define MAX_ITER	 1000	//maximum number of iteration

//Initialize all Tensors used in MERA and also give the random elements of W and U in all layers.
//The elements of W and U of all layers are randomly given unitary operators in the beginning.
void initTensors(Tensor** Rho, Tensor** Obs, Tensor** W, Tensor** U);
int main(){
	MEM = 0;
	double	lambda = 0;
	int 	tau;				//layer index
	//Here We assume taht the system is translational invariant
	Tensor* Rho[LAYER_NUM];		//H[tau] is the tau'th layer of H
	Tensor* H[LAYER_NUM];
	Tensor* W[LAYER_NUM-1];
	Tensor* U[LAYER_NUM-1];
	
	initTensors(Rho, H, W, U);
	initNetwork((char*)"../MERALib/diasNet");
	clock_t start, step, end;
	start = clock();	
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};
	addElem(H[0], H_elem);
	//Setting Hamiltonians in every layer by ascending operator
	for(tau = 0; tau<LAYER_NUM-1; tau++){
		Ascend(H[tau], W[tau], U[tau], H[tau+1]);
	}
	//Calculate the highest pure density matrix by leaving the state with least energy in H[LAYER_NUM]
	//and filter out those high energy states.
	double Ei, Ef;
	printf("status = %d\n", H[LAYER_NUM-1]->status);
	Filter(H[LAYER_NUM-1], Rho[LAYER_NUM-1], &Ei);
	printf("Ei = %f\n", Ei);
	//Setting density matrice in every layer by descending operator
	//start from highest pure density just calculated above.
	printf("Descending...\n");
	for(tau = LAYER_NUM - 2; tau >= 0; tau--)
		Descend(Rho[tau+1], W[tau], U[tau], Rho[tau]);

	//Do iteration to find the ground state and ground state energy
	int iter = 0;	//counter for iterations
	double delta;
	bool flag;
	do{
		//Update W and U in each layer
		for(tau = 0; tau<LAYER_NUM - 1; tau++){
			for(int i = 0; i<Q_ONE; i++){
				UpdateW(W[tau], Rho[tau+1], H[tau], U[tau]);
				UpdateU(U[tau], Rho[tau+1], H[tau], W[tau]);
			}
			//use updated U and W of layer tau to calculate new H in layer tau+1
			Ascend(H[tau], W[tau], U[tau], H[tau+1]);
		}
		Filter(H[LAYER_NUM-1], Rho[LAYER_NUM-1], &Ef);
		printf("Ef = %f\n", Ef);
		//Renew density matrix
		for(tau = LAYER_NUM - 2; tau >= 0; tau--)
			Descend(Rho[tau+1], W[tau], U[tau], Rho[tau]);
		delta = fabs(Ef - Ei);
		Ei = Ef;
		iter++;
		step = clock();
		printf("Iteration %d\telapsed time = %-8.2f secs\tE = %12.10f\n", iter, (float)(step - start)/CLOCKS_PER_SEC, Ef);
	}while(!(flag = delta<=DELTA) && iter<MAX_ITER);
	end = clock();
	printf("Total elapsed time: %-8.2f secs.(CPU time)\n", (float)(end - start)/CLOCKS_PER_SEC);
}

void initTensors(Tensor** Rho, Tensor** Obs, Tensor** W, Tensor** U){
	int tau;
	//Set dimension of bond of each layer.
	int LayDim[LAYER_NUM];
	for(int bound = (int)ceil(log(CHI)/log(SITE_DIM)/W_OUT), tau = 0; tau <LAYER_NUM; tau++){
		if(tau == 0)
			LayDim[tau] = SITE_DIM;
		else if(tau < bound)
			LayDim[tau] = (int)pow(SITE_DIM, W_OUT*tau);
		else
			LayDim[tau] = CHI;
	}

	//Initialize and add elements to each layer's  U, W
	int Wbd[W_IN+W_OUT];	//bond dimension of W
	int Ubd[U_IN+U_OUT];
	for(tau = 0; tau < LAYER_NUM - 1; tau++){
		for(int i = 0; i < W_IN; i++)
			Wbd[i] = LayDim[tau+1];
		for(int i = W_IN; i < W_IN + W_OUT; i++)
			Wbd[i] = LayDim[tau];
		W[tau] = (Tensor*)calloc(1, sizeof(Tensor));
		initTensor(W[tau], W_IN + W_OUT, Wbd);
		OrthoRandomize(W[tau], W_IN);
		
		for(int i = 0; i < U_IN + U_OUT; i++)
			Ubd[i] = LayDim[tau];
		U[tau] = (Tensor*)calloc(1, sizeof(Tensor));
		initTensor(U[tau], U_IN + U_OUT, Ubd);
		OrthoRandomize(U[tau], U_IN);
	}

	//Initialize each layer's Rho and Obs
	int Obsbd[OBS_BOND_NUM];		//bond dimensions of Obs
	int Rhobd[OBS_BOND_NUM];
	for(tau = 0; tau<LAYER_NUM; tau++){
		for(int i = 0; i<OBS_BOND_NUM; i++q
			Rhobd[i] = LayDim[tau];		
		Rho[tau] = (Tensor*)calloc(1, sizeof(Tensor));
		initTensor(Rho[tau], OBS_BOND_NUM, Rhobd);
		
		for(int i = 0; i<OBS_BOND_NUM; i++)
			Obsbd[i] = LayDim[tau];	
		Obs[tau] = (Tensor*)calloc(1, sizeof(Tensor));
		initTensor(Obs[tau], OBS_BOND_NUM, Obsbd);
	}	
}

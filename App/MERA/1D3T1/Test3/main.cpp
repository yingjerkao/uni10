#include <iostream>
#include <time.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include "MERA_Operator2.cpp"
#define FERMION 1
#define OBS_BOND_NUM 4		//the bond number of observable, here observable is H.
#define SITE_DIM	 2		//the Hilbert space dimension of a particle
#define SITE_NUM	 162
#define LAYER_NUM	 (int)ceil(log(SITE_NUM)/log(W_OUT))
#define CHI 		 22		//maximum bond dimension
#define DATA_NUM 	 101		//number of data points; number of different transverse field values.
#define HxR	 		 1000			//the ration of external field in z direction and in x direction
//Iteration variables
#define Q_ONE		 10		//times for update U and W in one iteration
#define DELTA	  (1.0E-8)	//precision of ground state energy, used in doing iteration
#define MAX_ITER	 1000	//maximum number of iteration
#define ID 1

#define W_IN 1
#define W_OUT 3
#define U_IN 2
#define U_OUT 2
#define OBS_IN 2
#define OBS_OUT 2

void initTensors(vector<UniTensor**> Rhos, vector<UniTensor**> Obs, vector<UniTensor**> Ws, vector<UniTensor**> Us);

int main(){
	Qnum q10(1, PRT_EVEN);
	Qnum q_10(-1, PRT_EVEN);
	Qnum q30(3, PRT_EVEN);
#ifdef FERMION
	Qnum q_11(PRTF_ODD, -1, PRT_EVEN);
	Qnum q11(PRTF_ODD, 1, PRT_EVEN);
	Qnum q_31(PRTF_ODD, -3, PRT_EVEN);
#else
	Qnum q_11(-1, PRT_ODD);
	Qnum q11(1, PRT_ODD);
	Qnum q_31(-3, PRT_ODD);
#endif

	char *basedir = (char*)"TransIsing/";
	char dirbuf[32];
	sprintf(dirbuf, "N%dChi%d-%d/", SITE_NUM, CHI, ID);
	char dir[64];
	strcpy(dir, basedir);
	strcat(dir, dirbuf);
	char *ret_fn = (char*)"result";
	char *out_fn = (char*)"out";
	char ret_path[128], out_path[128];
	strcpy(ret_path, dir); strcat(ret_path, ret_fn);
	strcpy(out_path, dir); strcat(out_path, out_fn);
	mkdir(dir, 0755);
	//The system is 1D Ising with Transverse Magnetic Field.
	double lamb[DATA_NUM + 5];	//magnitude of transverse field
	double lambda;
	int tau;	//layer index
	vector<UniTensor**> Rhos;
	vector<UniTensor**> Hs;
	vector<UniTensor**> Ws;
	vector<UniTensor**> Us;
	int dp;

	initTensors(Rhos, Hs, Ws, Us);
	Network asdL("../Diagrams/AscendL");
	Network asdC("../Diagrams/AscendC");
	Network asdR("../Diagrams/AscendR");
	Network desL("../Diagrams/DescendL");
	Network desC("../Diagrams/DescendC");
	Network desR("../Diagrams/DescendR");
	Network W1EnvL("../Diagrams/W1EnvL");
	Network W1EnvC("../Diagrams/W1EnvC");
	Network W1EnvR("../Diagrams/W1EnvR");
	Network W2EnvL("../Diagrams/W2EnvL");
	Network W2EnvC("../Diagrams/W2EnvC");
	Network W2EnvR("../Diagrams/W2EnvR");
	Network UEnvL("../Diagrams/UEnvL");
	Network UEnvC("../Diagrams/UEnvC");
	Network UEnvR("../Diagrams/UEnvR");
	clock_t start, step, end;
	start = clock();

	//Givivg lamda
	for(int dp = 0; dp < 5; dp++)
		lamb[dp] = 0.1*dp;
	for(int dp = 5; dp<DATA_NUM + 5; dp++)
		lamb[dp] = 0.45+0.01*dp;
	FILE* ret = fopen(ret_path, "w");
	FILE* out = fopen(out_path, "w");
	fprintf(out, "Q_ONE = %d\nDELTA = %e\n\n", Q_ONE, DELTA);
	for(int dp = 0; dp<DATA_NUM; dp++){
		lambda = lamb[dp];
		double H_elem[] = {lambda,  0,  0,   -1,
							 0   ,  0, -1,    0,
							 0   , -1,  0,    0,
							-1   ,  0,  0, -lambda};
		
		Hs[0].addRawElem(H_elem);
		//Setting Hamiltonians in every layer by ascending operator
		for(tau = 0; tau<LAYER_NUM-1; tau++)
			Hs[tau+1] = Ascend(&Hs[tau], &Ws[tau], &Us[tau], &asdL, &asdC, &asdR);
		cout << "End of Ascend\n";

		//Calculate the highest pure density matrix by leaving the state with least energy in Hs[LAYER_NUM]
		//and filter out those high energy states.
		double Ei, Ef;
		Rhos[LAYER_NUM-1] = Filter(&Hs[LAYER_NUM-1], &q10, &Ei);

		//Setting density matrice in every layer by descending operator
		//start from highest pure density just calculated above.
		for(tau = LAYER_NUM - 2; tau >= 0; tau--)
			Rhos[tau] = Descend(&Rhos[tau+1], &Ws[tau], &Us[tau], &desL, &desC, &desR);
		
		//Do iteration to find the ground state and ground state energy
		int iter = 0;	//counter for iterations
		double delta;
		bool flag;
		fprintf(out, "\n\nlambda = %-4.2f\n", lambda);
		do{
			//Update W and U in each layer
			for(tau = 0; tau<LAYER_NUM - 1; tau++){
				for(int i = 0; i<Q_ONE; i++){
					Ws[tau] = UpdateW(&Ws[tau], &Rhos[tau+1], &Hs[tau], &Us[tau], &W1EnvL, &W1EnvC, &W1EnvR, &W2EnvL, &W2EnvC, &W2EnvR);
					Us[tau] = UpdateU(&Us[tau], &Rhos[tau+1], &Hs[tau], &Ws[tau], &UEnvL, &UEnvC, &UEnvR);
				}
				//use updated U and W of layer tau to calculate new H in layer tau+1
				Hs[tau+1] = Ascend(&Hs[tau], &Ws[tau], &Us[tau], &asdL, &asdC, &asdR);
			}
			Rhos[LAYER_NUM-1] = Filter(&Hs[LAYER_NUM-1], &q10, &Ef);
			//Renew density matrix
			for(tau = LAYER_NUM - 2; tau >= 0; tau--)
				Rhos[tau] = Descend(&Rhos[tau+1], &Ws[tau], &Us[tau], &desL, &desC, &desR);
			delta = fabs(Ef - Ei);
			Ei = Ef;
			iter++;
			step = clock();
			fprintf(out, "Iteration %d\telapsed time = %-8.2f secs\tE = %12.10f\n", iter, (float)(step - start)/CLOCKS_PER_SEC, Ef);
			fflush(out);
		}while(!(flag = delta<=DELTA) && iter<MAX_ITER);
		
		if(flag)
			fprintf(ret, "lambda = %-3.2f	GE = %-16.14f	Total %d times of iteration.\n\n", lambda, Ef, iter);
		else
			fprintf(ret, "\nlambda = %-3.2f	It fails to converge!!!!!!!!!!!!!!!!!!\n\n", lambda);
		fflush(ret);

		//Output Tensors
		char Rho_path[50], H_path[50], W_path[50], U_path[50];
		char *Rho_fn = (char*)"Rho";
		char *H_fn = (char*)"H";
		char *W_fn = (char*)"W";
		char *U_fn = (char*)"U";
		char tensor_dir[40];
		char buf[6];
		strcpy(tensor_dir, dir);
		strcat(tensor_dir, (char*)"lambda_");
		sprintf(buf, "%4.2f", lambda);
		strcat(tensor_dir, buf);
		strcat(tensor_dir, (char*)"/");
		mkdir(tensor_dir, 0755);
		for(int i = 0; i < LAYER_NUM; i++){
			char taubuf[2];
			sprintf(taubuf, "%d", i);
			strcpy(Rho_path, tensor_dir); 
			strcat(Rho_path, Rho_fn);
			strcat(Rho_path, taubuf);
			strcpy(H_path, tensor_dir); 
			strcat(H_path, H_fn);
			strcat(H_path, taubuf);
			Rhos[i].save(Rho_path);
			Hs[i].save(H_path);
		}
		for(int i = 0; i < LAYER_NUM - 1; i++){
			char taubuf[2];
			sprintf(taubuf, "%d", i);
			strcpy(W_path, tensor_dir); 
			strcat(W_path, W_fn);
			strcat(W_path, taubuf);
			strcpy(U_path, tensor_dir);
			strcat(U_path, U_fn);
			strcat(U_path, taubuf);
			Ws[i].save(W_path);
			Us[i].save(U_path);
		}
	}
	Us.clear();
	Hs.clear();
	Ws.clear();
	Rhos.clear();
	end = clock();
	fprintf(out, "Total elapsed time: %-8.2f secs.(CPU time)", (float)(end - start)/CLOCKS_PER_SEC);
	fclose(out);
	fclose(ret);
}



void initTensors(vector<UniTensor**> Rhos, vector<UniTensor**> Obs, vector<UniTensor**> Ws, vector<UniTensor**> Us){
	Qnum q10(1, PRT_EVEN);
	Qnum q_10(-1, PRT_EVEN);
	Qnum q30(3, PRT_EVEN);
#ifdef FERMION
	Qnum q_11(PRTF_ODD, -1, PRT_EVEN);
	Qnum q11(PRTF_ODD, 1, PRT_EVEN);
	Qnum q_31(PRTF_ODD, -3, PRT_EVEN);
#else
	Qnum q_11(-1, PRT_ODD);
	Qnum q11(1, PRT_ODD);
	Qnum q_31(-3, PRT_ODD);
#endif

	vector<Bond> bonds;
	vector<Qnum> qnums;
	qnums.push_back(q10);qnums.push_back(q_11);
	Bond bdi(BD_IN, qnums);
	Bond bdo(BD_OUT, qnums);
	Bond bdi1;
	Bond bdo1;
	int dim = 0;
	map<Qnum, int> degs;
	
	bonds.push_back(bdi);
	bonds.push_back(bdi);
	bonds.push_back(bdo);
	bonds.push_back(bdo);

	UniTensor U;
	UniTensor H0;
	UniTensor W1;
	UniTensor Rho;
	
	for(int tau = 0; tau < LAYER_NUM - 1; tau++){
		qnums.clear();

		U.assign(bonds);
		U.orthoRand();
		Us.push_back(U);
		H0.assign(bonds);
		Obs.push_back(H0);
		bonds.clear();
	
		bdi1 = bdi;
		bdi1.combine(bdi);
		bdi1.combine(bdi);
		degs = bdi1.degeneracy();
		dim = 0;
		for(map<Qnum,int>::iterator it = degs.begin(); it != degs.end(); ++it)
			dim += it->second;
		if (dim > CHI) {
			qnums.push_back(q30);qnums.push_back(q30);qnums.push_back(q30);qnums.push_back(q30);
			qnums.push_back(q11);qnums.push_back(q11);qnums.push_back(q11);qnums.push_back(q11);
			qnums.push_back(q11);qnums.push_back(q11);qnums.push_back(q11);
			qnums.push_back(q_10);qnums.push_back(q_10);qnums.push_back(q_10);qnums.push_back(q_10);
			qnums.push_back(q_10);qnums.push_back(q_10);qnums.push_back(q_10);
			qnums.push_back(q_31);qnums.push_back(q_31);qnums.push_back(q_31);qnums.push_back(q_31);
			bdi1.assign(BD_IN,qnums);
			bdo1.assign(BD_OUT,qnums);
		}
		else{
			bdo1 = bdo;
			bdo1.combine(bdo);
			bdo1.combine(bdo);
		}
		degs.clear();
		
		bonds.push_back(bdi1);
		bonds.push_back(bdo);
		bonds.push_back(bdo);
		bonds.push_back(bdo);
		W1.assign(bonds);
		W1.orthoRand();
		Ws.push_back(W1);
		bonds.clear();
		
		bonds.push_back(bdi1);
		bonds.push_back(bdi1);
		bonds.push_back(bdo1);
		bonds.push_back(bdo1);
		Rho.assign(bonds);
		Rho.orthoRand();
		Rhos.push_back(Rho);

		bdi = bdi1;
		bdo = bdo1;
	}

}

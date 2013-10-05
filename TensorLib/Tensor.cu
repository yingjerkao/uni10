#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "cublas.h"
#include "cula.h"
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <boost/random.hpp>
using namespace boost;
#ifdef SEED
	mt19937 rng(SEED);
#else
	mt19937 rng(777);
#endif
uniform_01<mt19937> uni01_sampler(rng);
#define THREADMAX 10
#define BLOCKMAX 2
using namespace std;
#include "myDgemm.cu"
#include "Symmetry.cpp"
#define GPUMEM 2000000000
size_t MEM;

//Status of Tensors
#define INIT 1			//initialized
#define HAVELABEL 2		//tensor with label
#define HAVEELEM 4		//tensor with elements
#define DISPOSABLE 	8	//The elements of tensor is disposable through 'clone', 'reshapeClone', 'reshapeElem'. The elements of disposable Tensor is read-only.
#define ELEMFREED 	16	//the memory space of elements is freed
#define HAVESYM		32
#define RBLOCK 		64	//the memory space of elements is freed
#define CBLOCK		128
#define ONBLOCK 	192	//RBLOCK + CBLOCK 
#define ONGPU		(1<<30)

typedef struct{
	int status;		//Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements, 7 initialized completely.
	int bondNum;		//Tensor type
	int64_t elemNum;		//Number of elements in elem[] array
	int *bondDim;		//Dimemsions of the bonds, bonDim[0] is the high order bits in elem array.
	int *label;
	double *elem;		//Array of elements, allocated on Device(GPU)
	Symmetry *sym;
}Tensor;

//Initialization and memory allocatation
//Note that ,in matrix case, tensor's bondDim = {3,2} means it is a 3x2 matrix;
void initTensor(Tensor* T, int BondNumber, int *BondDimension, int type);

//Add new label to a tensor
void addLabel(Tensor* T, int* Label);

#define HtoD 0
#define DtoD 1
//Add elememnts to a tensor.
//Two ways to do it, copy Elements from main memory(type = HtoD) or copy Elements from device memory(type = DtoD)
void addElem(Tensor* T, double* Elements, int type);

//Output a Tensor struct
void outTensorf(Tensor* T, char* Filename);

//Fill a Tensor struct from file, to initialized an uninitilized Tensor.
void inTensorf(Tensor* T, char* Filename);

//Clone the input Tensor.
void clone(Tensor* Tin, Tensor* Tout);

//Reshape the element array of T to the given bond order. order[3] = 1 means the present forth bond should be the second bond.
//This function doesn't change the structure of tensor T; it just outputs the reshaped element array to NewElem.
void reshapeElem(Tensor* T, int* BondOrder, double *NewElem);
void reshapeElemH(Tensor* T, int* order, double *newElem);

//Reshape the Tin tensor to the given bond order. order[3] = 1 means the present forth bond in Tin should be the second bond in Tout.
//This function doesn't change the structure of tensor Tin; it just outputs the reshaped Tensor to Tout.
void reshapeClone(Tensor* Tin, int* order, Tensor* Tout);

//Find the matched labels in Ta and Tb, contracted the corresponding matched bonds. Tc = Ta*Tb
void contraction(Tensor* Ta, Tensor* Tb, Tensor* Tc);

//Reclaim memory
void recycle(Tensor* T);


//Print out a m*n matrix.
void print_matrix(char* desciption, int m, int n, double* matrix);
void print_matrix( char* desc, int m, int n, int* a);

//transpose a matrix represented in C to the one represented in Fortran
void CtoF(int m, int n, double* Mc, double* Mf);
//transpose a matrix represented in Fortran to the one represented in C
void FtoC(int m, int n, double* Mf, double* Mc);

void Matrix_Product(double* A, int At, double* B, int Bt, int m, int k, int n, double* C);

void print_matrix(double* Matrix, int m, int n);

void Substract_Identity(double* Matrix, int size, double factor);
void initTensor(Tensor* T, int bn, int *bd, int type = 0){
	T->bondNum = bn;
	int i;
	T->bondDim = (int*)realloc(T->bondDim, sizeof(int) * bn);
	memcpy((T->bondDim), bd, sizeof(int)*bn);
	T->label = (int*)realloc(T->label, sizeof(int) * bn);
	int oldElemNum = T->elemNum;
	T->elemNum = 1;
	for(i = 0; i < bn; i++)
		T->elemNum *= bd[i];
	if(type==1)
		T->elem = (double*)realloc(T->elem, sizeof(double) * T->elemNum);
	else{
		if(T->status & INIT){
			cudaFree(T->elem);
			MEM -= oldElemNum;
		}
		assert(cudaMalloc((void**)&(T->elem), sizeof(double) * T->elemNum) == cudaSuccess);
		T->status |= ONGPU;
		MEM += T->elemNum;
	}
	T->status |= INIT;
}

void addLabel(Tensor* T, int* lb){
	assert(T->status & INIT);
	memcpy(T->label, lb, sizeof(int)*(T->bondNum));
	T->status |= HAVELABEL;
}

void addElem(Tensor* T, double* elem, int type){
	assert(T->status & INIT);
	assert(T->status & ONGPU);
	assert(!(T->status & ELEMFREED));
	if(type == HtoD)	//elem is a cpu memory pointer
		assert(cudaMemcpy(T->elem, elem, sizeof(double) * (T->elemNum), cudaMemcpyHostToDevice) == cudaSuccess);
	else if(type == DtoD)	//elem is a device memory pointer
		assert(cudaMemcpy(T->elem, elem, sizeof(double) * (T->elemNum), cudaMemcpyDeviceToDevice) == cudaSuccess);
	else{
		fprintf(stderr, "Function addElem() in Tensor.cu:\nIllegal data transfer type\n");
		exit(0);
	}
	T->status |= HAVEELEM;
}

/*-----------------Tensor output format------------------*/
/* The status of Tensor(integer)
 * Number of bonds(integer)
 * Number of elements(integer)
 * Array of bond dimensions(integer array)
 * Array of label, if the tensor has it.(integer array)
 * Array of elements, if the tensor has it.(double array)
 */
void outTensorf(Tensor* T, char* fn){
	assert(T->status & ONGPU);
	FILE *fp = fopen(fn, "w");
	fprintf(fp, "%d\n", T->status);
	fprintf(fp, "%d\n", T->bondNum);
	fprintf(fp, "%d\n", T->elemNum);
	int64_t i;
	for(i = 0; i<T->bondNum; i++)
		fprintf(fp, "%d ",T->bondDim[i]);
	fprintf(fp, "\n");
	if(T->status & HAVELABEL){
		for(i = 0; i<T->bondNum; i++)
			fprintf(fp, "%d ",T->label[i]);
		fprintf(fp, "\n");
	}
	if(T->status & HAVEELEM){
		double *elem = (double*)malloc(T->elemNum * sizeof(double));
		assert(cudaMemcpy(elem, T->elem, T->elemNum * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess);
		for(i = 0; i <T->elemNum; i++)
			fprintf(fp, "%f ",elem[i]);
		fprintf(fp, "\n");
		free(elem);
	}
	fclose(fp);
}

void outTensor(Tensor* T, char* fn){
	assert(T->status & ONGPU);
	FILE *fp = fopen(fn, "w");
	fwrite(&T->status, sizeof(int), 1, fp);
	fwrite(&T->bondNum, sizeof(int), 1, fp);
	fwrite(&T->elemNum, sizeof(int64_t), 1, fp);
	fwrite(T->bondDim, sizeof(int), T->bondNum, fp);
	if(T->status & HAVELABEL)
		fwrite(T->label, sizeof(int), T->bondNum, fp);
	if(T->status & HAVEELEM){
		double *elem = (double*)malloc(T->elemNum * sizeof(double));
		assert(cudaMemcpy(elem, T->elem, T->elemNum * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess);
		fwrite(elem, sizeof(double), T->elemNum, fp);
		free(elem);
	}
	fclose(fp);
}

void inTensorf(Tensor* T, char* fn){
	FILE *fp = fopen(fn, "r");
	assert(fp != NULL);
	if(T->status & INIT)
		cudaFree(T->elem);
	fscanf(fp, "%d", &(T->status));
	T->status |= ONGPU;
	fscanf(fp, "%d", &(T->bondNum));
	fscanf(fp, "%ld", &(T->elemNum));
	T->bondDim = (int*)realloc(T->bondDim, sizeof(int)*(T->bondNum));
	T->label = (int*)realloc(T->label, sizeof(int)*(T->bondNum));
	assert(cudaMalloc((void**)&(T->elem), sizeof(double) * T->elemNum) == cudaSuccess);
	int64_t i;
	for(i = 0; i<T->bondNum; i++)
		fscanf(fp, "%d", &(T->bondDim[i]));
	if(T->status & HAVELABEL)
		for(i = 0; i<T->bondNum; i++)
			fscanf(fp, "%d", &(T->label[i]));
	if(T->status & HAVEELEM){
		double *elem = (double*)malloc(T->elemNum * sizeof(double));
		for(i = 0; i < T->elemNum; i++)
			fscanf(fp, "%lf", &(elem[i]));
		assert(cudaMemcpy(T->elem, elem, T->elemNum * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess);
		free(elem);
	}
	fclose(fp);
}

void inTensor(Tensor* T, char* fn){
	FILE *fp = fopen(fn, "r");
	assert(fp != NULL);
	if(T->status & INIT)
		cudaFree(T->elem);
	fread(&T->status, sizeof(int), 1, fp);
	T->status |= ONGPU;
	fread(&T->bondNum, sizeof(int), 1, fp);
	fread(&T->elemNum, sizeof(int64_t), 1, fp);
	T->bondDim = (int*)realloc(T->bondDim, sizeof(int)*(T->bondNum));
	T->label = (int*)realloc(T->label, sizeof(int)*(T->bondNum));
	assert(cudaMalloc((void**)&(T->elem), sizeof(double) * T->elemNum) == cudaSuccess);
	fread(T->bondDim, sizeof(int), T->bondNum, fp);
	if(T->status & HAVELABEL)
		fread(T->label, sizeof(int), T->bondNum, fp);
	if(T->status & HAVEELEM){
		double *elem = (double*)malloc(T->elemNum * sizeof(double));
		fread(elem, sizeof(double), T->elemNum, fp);
		assert(cudaMemcpy(T->elem, elem, T->elemNum * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess);
		free(elem);
	}
	fclose(fp);
}

void clone(Tensor* Tin, Tensor* Tout){
	assert(Tin->status & ONGPU);
	initTensor(Tout, Tin->bondNum, Tin->bondDim);
	if(Tin->status & HAVELABEL)
		addLabel(Tout, Tin->label);
	if(Tin->status & HAVEELEM){
		addElem(Tout, Tin->elem, DtoD);	
		if(Tin->status & DISPOSABLE){
			cudaFree(Tin->elem);
			Tin->status |= ELEMFREED;
		}
	}
}

void recycle(Tensor* T){
	if(T->status & INIT){
		free(T->bondDim);
		free(T->label);
		if((T->status & ELEMFREED) == 0){
			if(T->status & ONGPU){
				cudaFree(T->elem);
				MEM -= T->elemNum;
			}
			else
				free(T->elem);
		}
	}
	free(T);
}

void reshapeElemH(Tensor* T, int* order, double *newElem){
	assert((T->status & ONGPU) == 0);
	assert(T->status & HAVEELEM);
	int num = T->bondNum;
	int64_t i, j;
	//if order = [0 1 2 3], there's no need to reshape
	int done = 1;
	for(i = 0; i<num; i++)
		if(order[i] != i){
			done = 0;
			break;
		}
	if(done)
		memcpy(newElem, T->elem, sizeof(double)*(T->elemNum));
	else{
		int newDim[num];
		for(i = 0; i<num; i++)
			newDim[order[i]] = T->bondDim[i];
		int64_t newOffset[num];
		int64_t tempOffset[num];
		tempOffset[num-1] = 1;
		for(i = num-2; i >= 0; i--)
			tempOffset[i] = tempOffset[i+1]*newDim[i+1];
		for(i = 0; i<num; i++)
			newOffset[i] = tempOffset[order[i]];
		int64_t offset = 0;	//the offset of newElem.
		int idx[num];	//the indice of oldElem.
		memset(idx, 0, sizeof(int)*num);
		newElem[0] = T->elem[0];
		for(i = 1; i<T->elemNum; i++){
			for(j = num-1; j>=0; j--){
				if((idx[j]+1)/T->bondDim[j]){
					idx[j] = 0;
					offset -= (T->bondDim[j]-1)*newOffset[j];
				}
				else{
					idx[j]++;
					offset += newOffset[j];
					break;
				}
			}
			newElem[offset] = T->elem[i];
		}
	}	
	if(T->status & DISPOSABLE){
		free(T->elem);
		T->status |= ELEMFREED;
	}
}

__global__ void reshape(double* oldElem, int bondNum, int64_t elemNum, int64_t* offset, double* newElem, int64_t threadNum){
	int64_t oldIdx = blockIdx.x * blockDim.x + threadIdx.x;
	int64_t idx;
	int64_t newIdx;
	int cnt = 0;
	while(1){
		idx = oldIdx;
		newIdx = 0;
		if(idx < elemNum){
			for(int i = 0; i < bondNum; i++){
				newIdx += (idx/offset[i]) * offset[bondNum + i];
				idx = idx % offset[i];
			}
			newElem[newIdx] = oldElem[oldIdx];
		}
		cnt++;
		if((elemNum - threadNum * cnt) <= 0)
			break;
		oldIdx += threadNum;
	}

}

void reshapeElem(Tensor* T, int* order, double *newElem){
	assert(T->status & ONGPU);
	assert(T->status & HAVEELEM);
	int num = T->bondNum;
	int i;
	//if order = [0 1 2 3], there's no need to reshape
	int done = 1;
	for(i = 0; i<num; i++)
		if(order[i] != i){
			done = 0;
			break;
		}
	if(done)
		assert(cudaMemcpy(newElem, T->elem, sizeof(double)*(T->elemNum), cudaMemcpyDeviceToDevice) == cudaSuccess);
	else{
		int newDim[num];
		for(i = 0; i < num; i++)
			newDim[order[i]] = T->bondDim[i];
		int64_t OldNewOffset[2 * num];
		int64_t tempOffset[num];
		tempOffset[num-1] = 1;
		for(i = num-2; i >= 0; i--)
			tempOffset[i] = tempOffset[i+1] * newDim[i+1];
		OldNewOffset[num-1] = 1;
		for(i = num-2; i >= 0; i--)
			OldNewOffset[i] = OldNewOffset[i+1] * T->bondDim[i+1];
		for(i = 0; i < num; i++)
			OldNewOffset[num + i] = tempOffset[order[i]];
		int64_t *D_offset;
		assert(cudaMalloc((void**)&D_offset, 2 * sizeof(int64_t) * num) == cudaSuccess);
		assert(cudaMemcpy(D_offset, OldNewOffset, 2 * sizeof(int64_t) * num, cudaMemcpyHostToDevice) == cudaSuccess);
		int blockSize = THREADMAX;
		int gridSize = (T->elemNum + THREADMAX - 1) / THREADMAX;
		if(gridSize > BLOCKMAX){
			gridSize = BLOCKMAX;
		}
		reshape<<<gridSize, blockSize>>>(T->elem, T->bondNum, T->elemNum, D_offset, newElem, gridSize * blockSize);
		cudaFree(D_offset);
	}
	if(T->status & DISPOSABLE){
		cudaFree(T->elem);
		T->status |= ELEMFREED;
		MEM -= T->elemNum;
	}
}

void reshapeClone(Tensor* Tin, int* order, Tensor* Tout){
	assert(Tin->status & ONGPU);
	int i;
	int newDim[Tin->bondNum];
	for(i = 0; i<Tin->bondNum; i++)
		newDim[order[i]] = Tin->bondDim[i];
	initTensor(Tout, Tin->bondNum, newDim);
	if(Tin->status & HAVELABEL){
		for(i = 0; i<Tout->bondNum; i++)
			Tout->label[order[i]] = Tin->label[i];
		Tout->status |= HAVELABEL;
	}
	if(Tin->status & HAVEELEM){
		reshapeElem(Tin, order, Tout->elem);
		Tout->status |= HAVEELEM;
	}
}

void contraction(Tensor* Ta, Tensor* Tb, Tensor* Tc){
	if(Ta->elemNum < Tb->elemNum){
		Tensor* Tmp = Tb;
		Tb = Ta;
		Ta = Tmp;
	}
	assert(Ta->status&HAVEELEM && Tb->status&HAVEELEM);
	//Tensor contraction. Tc = Ta*Tb.
	int M = 1, N = 1, K = 1;	//Ta = MxK, Tb = KxN, Tc = MxN
	int i, j, p;
	int num = 0;	//number of bonds which are to be contracted.
	int aBondNum = Ta->bondNum;
	int bBondNum = Tb->bondNum;
	int aMatch[aBondNum];
	int bMatch[bBondNum];
	//set the reshape order, order[3] = 1 means the present forth bond should be the second bond.
	int aOrder[aBondNum];
	int bOrder[bBondNum];
	memset(aMatch, 0, sizeof(int)*aBondNum);
	memset(bMatch, 0, sizeof(int)*bBondNum);
	p = 1;
	for(i = 0; i<aBondNum; i++)
		for(j = 0; j<bBondNum; j++)
			if(Ta->label[i] == Tb->label[j]){
				//the ith bond of Ta match the jth bond of Tb.
				aMatch[i] = p;	
				bMatch[j] = p;
				num++;
				p++;
				if(Ta->bondDim[i] != Tb->bondDim[j])
						printf("Ta-bondDim[%d] = %d ; Tb-bondDim[%d] = %d\n", i, Ta->bondDim[i], j, Tb->bondDim[j]);
				assert(Ta->bondDim[i] == Tb->bondDim[j]);
			}
	//example for Ta->label=[12 29 7 43]; Tb->label=[43 33 12 19]  =>aMatch = [1 0 0 2], bMatch[2 0 1 0]
	//Generate a new Tensor C = A*B
	int cBondNum = aBondNum+bBondNum-2*num;
	int cBondDim[cBondNum];
	int cLabel[cBondNum];
	j = 0, p = 0;
	for(i = 0; i<aBondNum; i++){
		if(aMatch[i]){
			aOrder[i] = aBondNum-(num-aMatch[i])-1;
			K *= Ta->bondDim[i];
		}
		else{
			aOrder[i] = j;
			cBondDim[p] = Ta->bondDim[i];
			cLabel[p] = Ta->label[i];
			M *= Ta->bondDim[i];
			p++;
			j++;
		}
	}
	j = 0;
	for(i = 0; i<bBondNum; i++){
		if(bMatch[i])
			bOrder[i] = bMatch[i]-1;
		else{
			bOrder[i] = num+j;
			cBondDim[p] = Tb->bondDim[i];
			cLabel[p] = Tb->label[i];
			N *= Tb->bondDim[i];
			p++;
			j++;
		}
	}
	assert(p == cBondNum);
	//printf("MEM = %d\n", MEM);
	//printf("Ta->elemNum = %d\n", Ta->elemNum);
	//printf("Tb->elemNum = %d\n", Tb->elemNum);

	size_t chunk_max = GPUMEM / 5;
	size_t transfer_max = 1.8 * chunk_max;
	int mm_t = 0;
	int mp = 1, mq = 1, mc = 1;
	
	size_t Tc_mem = sizeof(double);
	for(int i = 0; i < cBondNum; i++)
		Tc_mem *= cBondDim[i];
	if(Tc_mem > chunk_max){
		mm_t = 1;	//Tc on CPU
		mc = ceil(float(Tc_mem) / chunk_max);
	}
	else
		mm_t = 0;	//Tc on GPU
	
	initTensor(Tc, cBondNum, cBondDim, mm_t);
	addLabel(Tc, cLabel);
	printf("Tc->elemNum = %d\n", Tc->elemNum);
	
	double *aNewElem, *bNewElem;
	if(Ta->status & ONGPU){
		MEM += Ta->elemNum;
		assert(cudaMalloc((void**)&aNewElem, sizeof(double) * Ta->elemNum) == cudaSuccess);
		reshapeElem(Ta, aOrder, aNewElem);
	}
	else{
		aNewElem = (double*)malloc(sizeof(double) * Ta->elemNum);
		reshapeElemH(Ta, aOrder, aNewElem);
		mp = ceil(float(Ta->elemNum * sizeof(double)) / transfer_max);
		mm_t |= 4;
	}
	if(Tb->status & ONGPU){
		MEM += Tb->elemNum;
		assert(cudaMalloc((void**)&bNewElem, sizeof(double) * Tb->elemNum) == cudaSuccess);
		reshapeElem(Tb, bOrder, bNewElem);
	}
	else{
		bNewElem = (double*)malloc(sizeof(double) * Tb->elemNum);
		reshapeElemH(Tb, bOrder, bNewElem);
		mq = ceil(float(Tb->elemNum * sizeof(double)) / transfer_max);
		mm_t |= 2;
	}
		
	if(mc > mp * mq)
		mp = ceil(float(Tc_mem) / transfer_max / mq);
	//printf("mm_t = %d\n", mm_t);
	//printf("mp = %d\n", mp);
	//printf("mq = %d\n", mq);
	//printf("M = %d, N = %d, K = %d\n", M, N, K);
	//printf("MEM = %d\n", MEM);
	mmtype types[] = {MM_DDD, MM_DDH, MM_DHD, MM_DHH, MM_HDD, MM_HDH, MM_HHD, MM_HHH};	
	myDgemm(mp, mq, M, N, K, aNewElem, bNewElem, Tc->elem, types[mm_t]);
	
	//cublasDgemm('N', 'N', N, M, K, 1, bNewElem, N, aNewElem, K, 0, Tc->elem, N);
	cublasStatus_t status = cublasGetError();
	assert(status == CUBLAS_STATUS_SUCCESS);
	Tc->status |= HAVEELEM;
	if(mm_t & 4)
		free(aNewElem);
	else{
		cudaFree(aNewElem);
		MEM -= Ta->elemNum;
	}
	if(mm_t & 2)
		free(bNewElem);
	else{
		cudaFree(bNewElem);
		MEM -= Tb->elemNum;
	}
	//printf("\n\n");
}

void contraction(Tensor* Ta, Tensor* Tb, Tensor* Tc, int* outLabel){
	contraction(Ta, Tb, Tc);
	int bondNum = Tc->bondNum;
	int order[bondNum];
	int hit = 0;
	for(int i = 0; i < bondNum; i++)
		for(int j = 0; j < bondNum; j++)
			if(Tc->label[i] == outLabel[j]){
				order[i] = j;
				hit ++;
			}
	assert(hit == bondNum);
	double *cNewElem;
	MEM += Tc->elemNum;
	assert(cudaMalloc((void**)&cNewElem, sizeof(double) * Tc->elemNum) == cudaSuccess);
	reshapeElem(Tc, order, cNewElem);
	addElem(Tc, cNewElem, DtoD);
	int oldDim[bondNum];
	for(int i = 0; i<Tc->bondNum; i++)
		oldDim[i] = Tc->bondDim[i];
	for(int i = 0; i<Tc->bondNum; i++)
		Tc->bondDim[order[i]] = oldDim[i];
	addLabel(Tc, outLabel);
	cudaFree(cNewElem);
	MEM -= Tc->elemNum;
}

/*Generate a set of row vectors which form a othonormal basis
  For the incoming Tensor, the number of row <= the number of column
 *Caution!! This function is used only for initialization. 
 *DO NOT used in the main part of your application, since it is not a well optimized function on GPU.
 */
void OrthoRandomize(Tensor *Tin, int row_bond_num){
	assert(Tin->status & ONGPU);
	int M = 1, N = 1;
	for(int i = 0; i < Tin->bondNum; i++){
		if(i < row_bond_num)
			M *= Tin->bondDim[i];
		else
			N *= Tin->bondDim[i];
	}
	assert(M <= N);
	double *random;
	assert(cudaMalloc((void**)&random, Tin->elemNum * sizeof(double)) == cudaSuccess);
	double *cpurandom = (double*)malloc(Tin->elemNum * sizeof(double));
	for (int i=0 ; i < Tin->elemNum ; i++){
		cpurandom[i] = uni01_sampler();
	}
	assert(cudaMemcpy(random, cpurandom, sizeof(double) * (Tin->elemNum), cudaMemcpyHostToDevice) == cudaSuccess);
	Tensor * tmpT = (Tensor*)calloc(1, sizeof(Tensor));
	int tmp_bondDim[Tin->bondNum];
	for(int i = 0; i < Tin->bondNum; i++)	//row_bond_num of Tin is col_bond_num of tmp
		tmp_bondDim[i] = Tin->bondDim[(i + row_bond_num) % Tin->bondNum];
	initTensor(tmpT, Tin->bondNum, tmp_bondDim);
	int min = M; //min = min(M,N)
	int ldA = M, ldu = M, ldvT = min;
	double *S, *u;
	assert(cudaMalloc((void**)&S, min * sizeof(double)) == cudaSuccess);
	assert(cudaMalloc((void**)&u, ldu * M * sizeof(double)) == cudaSuccess);
	//tmpT = u*S*vT
	culaDeviceDgesvd('A', 'S', M, N, random, ldA, S, u, ldu, tmpT->elem, ldvT);
	tmpT->status |= HAVEELEM;
	//reshape from Fortran format to C format
	int rsp_order[Tin->bondNum];
	for(int i = 0; i < Tin->bondNum; i++)
		rsp_order[i] = (i + row_bond_num)%Tin->bondNum;
	reshapeElem(tmpT, rsp_order, Tin->elem);
	Tin->status |= HAVEELEM;
	recycle(tmpT);
	cudaFree(S);
	cudaFree(u);
	cudaFree(random);
	free(cpurandom);
}

//threadNum is the total number of threads, min is the number of elements needed to be set to 1, min >= threadNum
__global__ void _identity(double* elem, int M, int N, int min, int threadNum){
	int64_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	int cnt = 0;
	while(1){
		if(idx < min)
			elem[idx * N + idx] = 1;
		cnt++;
		if((min - threadNum * cnt) <= 0)
			break;
		idx += threadNum;
	}
}
void Identity(Tensor* Tin, int row_bond_num){
	assert(Tin->status & ONGPU);
	int M = 1, N = 1;
	for(int i = 0; i < Tin->bondNum; i++){
		if(i < row_bond_num)
			M *= Tin->bondDim[i];
		else
			N *= Tin->bondDim[i];
	}
	int min;
	if(M < N)	min = M;
	else		min = N;
	
	assert(cudaMemset(Tin->elem, 0, Tin->elemNum * sizeof(double)) == cudaSuccess);
	int blockSize = THREADMAX;
	int gridSize = (min + THREADMAX - 1) / THREADMAX;
	if(gridSize > BLOCKMAX){
		gridSize = BLOCKMAX;
	}
	printf("min = %d\n", min);
	printf("threadNum = %d\n", gridSize * blockSize);
	_identity<<<gridSize, blockSize>>>(Tin->elem, M, N, min, gridSize * blockSize);
	Tin->status |= HAVEELEM;
}

//Tout = Tout + a * Ta
void addTensor(Tensor* T, double a, Tensor* Ta){
	assert(T->status & ONGPU);
	assert(T->status & HAVEELEM);
	assert(Ta->status & HAVEELEM);
	int64_t left = T->elemNum;
	assert(left == Ta->elemNum);
	int64_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		cublasDaxpy(chunk, a, Ta->elem + offset, 1, T->elem + offset, 1);
		offset += chunk;
		left -= INT_MAX;
	}
}

void print_tensor(Tensor* T, int row_bond_num){
	assert(T->status & ONGPU);
	int M = 1, N = 1;
	for(int i = 0; i < T->bondNum; i++){
		if(i < row_bond_num)
			M *= T->bondDim[i];
		else
			N *= T->bondDim[i];
	}
	int i, j;
	printf("\nRank %d Tensor\n", T->bondNum);
	printf("Bonds Dimemsion: ");
	for(i = 0; i < T->bondNum; i++)
		printf("%d ", T->bondDim[i]);
	printf("\n");
	if(T->status & HAVELABEL){
		printf("Label: ");
		for(i = 0; i < T->bondNum; i++)
			printf("%d ", T->label[i]);
		printf("\n");
	}
	double *elem = (double*)malloc(sizeof(double) * T->elemNum);
	cudaMemcpy(elem, T->elem, sizeof(double) * (T->elemNum), cudaMemcpyDeviceToHost);
	printf("\n");
	if(T->status & HAVEELEM)
		for(i = 0; i < M; i++) {
			for(j = 0; j < N; j++)
				printf("%6.3f", elem[i*N+j]);
			printf( "\n" );
			printf( "\n" );
		}
	printf("\n");
	free(elem);
}


/*------For Symmetry------*/
#define TOBLOCK 1
#define OFFBLOCK 0
void addSym(Tensor* T, int row_bond_num, Symmetry* S, int status = 0){
	//status=ONBLOCK means that, the element is block diagonalized now;
	assert(T->status & INIT);
	int M = 1, N = 1;
	for(int i = 0; i < T->bondNum; i++)
		if(i < row_bond_num)
			M *= T->bondDim[i];
		else
			N *= T->bondDim[i];
	assert(S->totRow == M);
	assert(S->totCol == N);
	T->sym = S;
	T->status |= HAVESYM;
	assert(status == 0 || status == RBLOCK || status == CBLOCK || status == ONBLOCK);
	T->status |= status;
}

void addTable(Symmetry* S, int mapNum, int* table, double* factor){
	assert(S->group.size());
	int tableSz = S->totRow > S->totCol ? S->totRow: S->totCol;	//max(S->totRow, S->totCol)
	S->mapNum = mapNum;
	if(S->table != NULL)
		free(S->table);
	S->table = (int*)calloc(mapNum * tableSz, sizeof(int));
	memcpy(S->table, table, sizeof(int)*(mapNum * tableSz));
	if(S->invTable != NULL)
		free(S->invTable);
	S->invTable = (int*)calloc(mapNum * tableSz, sizeof(int));
    for (int i=0; i< mapNum * tableSz; i++)
    	S->invTable[i]  = 0;
	
    #ifdef SANITY_CHECK
	int check[tableSz];
	memset(check, 0, tableSz * sizeof(int));
	for(int i = 0; i < tableSz * mapNum; i++){
		if(mapNum == 1)
			check[table[i] % tableSz] ++;
		else if(factor[i] != 0)
			check[table[i] % tableSz] ++;
	}
	for(int i = 0; i < tableSz; i++){
		if (check[i] < 1)
			printf("check[%d] = %d; too small!!!\n", i, check[i]);
	   	if (check[i] > mapNum)
		    printf("check[%d] = %d; too large!!!\n", i, check[i]);
		assert(check[i] >= 1);
		assert(check[i] <= mapNum);
	}
	//printf("CHECK OK!\n");
	#endif

	if(mapNum == 1)
		for(int i = 0; i < tableSz; i++)
			S->invTable[table[i]] = i;
	else{
		if(S->factor != NULL)
			free(S->factor);
		S->factor = (double*)calloc(mapNum * tableSz, sizeof(double));
		memcpy(S->factor, factor, sizeof(double)*(mapNum * tableSz));
		if(S->invFactor != NULL)
			free(S->invFactor);
		S->invFactor = (double*)calloc(mapNum * tableSz, sizeof(double));
        for (int i=0; i<tableSz*mapNum; i++)
			S->invFactor[i] = 0;

		int count[tableSz];
		int ustime[tableSz];
		memset(count, 0, tableSz * sizeof(int));
		memset(ustime, 0, tableSz * sizeof(int));
		int elemNum = tableSz * mapNum;
		for(int i = 0; i < elemNum; i++)
			if(factor[i] != 0) {
				int index = table[i] % tableSz;
				int i0 = i/mapNum;
				S->invTable[index * mapNum + count[index]] = i0 + ustime[i0] * tableSz;
				ustime[i0]++;
				S->invFactor[index * mapNum + count[index]] = factor[i];
				count[index]++;
			}
	}
}

//Caution!!! This is a unsafe function.
//Make sure you FREE the double pointer you get from it.
double* getBlock(Tensor* Tsrc, int row_bond_num, int groupNo){
	assert(Tsrc->status & HAVEELEM);
	assert(Tsrc->status & HAVESYM);
	assert(Tsrc->status & RBLOCK && Tsrc->status & CBLOCK);
	Symmetry* S = Tsrc->sym;
	assert(groupNo < S->group.size());
	Group g = S->group[groupNo];
	int M = 1, N = 1;
	for(int i = 0; i < Tsrc->bondNum; i++){
		if(i < row_bond_num)
			M *= Tsrc->bondDim[i];
		else
			N *= Tsrc->bondDim[i];
	}
	assert(M == S->totRow);
	assert(N == S->totCol);
	//double *outBlock = (double*)malloc(g.col * g.row * sizeof(double));
	double *outBlock;
	assert(cudaMalloc((void**)&outBlock, sizeof(double) * g.row * g.col) == cudaSuccess);
	for(int i = 0; i < g.row; i++)
		//memcpy(outBlock + i * g.col, Tsrc->elem + (i + g.row_begin) * N + g.col_begin, g.col * sizeof(double));
		cudaMemcpy(outBlock + i * g.col, Tsrc->elem + (i + g.row_begin) * N + g.col_begin, sizeof(double) * g.col, cudaMemcpyDeviceToDevice);

	return outBlock;
}

void putBlock(Tensor* Tout, int row_bond_num, int groupNo, double* src){
	assert(Tout->status & INIT);
	assert(Tout->status & HAVESYM);
	assert(Tout->status & RBLOCK && Tout->status & CBLOCK);
	Symmetry* S = Tout->sym;
	assert(groupNo < S->group.size());
	Group g = S->group[groupNo];
	int M = 1, N = 1;
	for(int i = 0; i < Tout->bondNum; i++){
		if(i < row_bond_num)
			M *= Tout->bondDim[i];
		else
			N *= Tout->bondDim[i];
	}
	assert(M == S->totRow);
	assert(N == S->totCol);
	for(int i = 0; i < g.row; i++)
		assert(cudaMemcpy(Tout->elem + (i + g.row_begin) * N + g.col_begin, src + i * g.col, sizeof(double) * g.col, cudaMemcpyDeviceToDevice) == cudaSuccess);
}

void symIdentity(Tensor* T, int row_bond_num){
	assert(T->status & INIT);
	assert(T->status & HAVESYM);
	assert(cudaMemset(T->elem, 0, T->elemNum * sizeof(double)) == cudaSuccess);
	T->status |= RBLOCK;
	T->status |= CBLOCK;
	Symmetry* S = T->sym;
	int tmpDim[2];
	Tensor* Ttmp = (Tensor*)calloc(1, sizeof(Tensor));
	for(int i = 0; i < S->group.size(); i++){
		Group g = S->group[i];
		tmpDim[0] = g.row;
		tmpDim[1] = g.col;
		initTensor(Ttmp, 2, tmpDim);
		Identity(Ttmp, 1);
		putBlock(T, row_bond_num, i, Ttmp->elem);
	}
	recycle(Ttmp);	
	T->status |= HAVEELEM;
}

void symOrthoRandomize(Tensor* T, int row_bond_num){
	assert(T->status & INIT);
	assert(T->status & HAVESYM);
	assert(cudaMemset(T->elem, 0, T->elemNum * sizeof(double)) == cudaSuccess);
	T->status |= RBLOCK;
	T->status |= CBLOCK;
	Symmetry* S = T->sym;
	int tmpDim[2];
	Tensor* Ttmp = (Tensor*)calloc(1, sizeof(Tensor));
	for(int i = 0; i < S->group.size(); i++){
		Group g = S->group[i];
		tmpDim[0] = g.row;
		tmpDim[1] = g.col;
		initTensor(Ttmp, 2, tmpDim);
		OrthoRandomize(Ttmp, 1);
		putBlock(T, row_bond_num, i, Ttmp->elem);
	}
	recycle(Ttmp);	
	T->status |= HAVEELEM;
}


__global__ void ColMap1(int* table, int M, int N, double* oldElem, double* newElem){
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	if (i<M && j<N)
		newElem[i*N + table[j]] = oldElem[i*N + j];
}

__global__ void ColMap2(int* table, double* factor, int M, int N, int mapNum, double* oldElem, double* newElem){
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	if (i<M && j< N*mapNum && factor[j]!=0){
		int index = table[j] % N;
		int count = table[j] / N;
		int off = i * N + count * N * M;
		newElem[off + index] = factor[j] * oldElem[i*N + j/mapNum];
	}
	__syncthreads();
}

__global__ void ColCollect(int M, int N, int k, double* Elem){
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	if (i<M && j<N){
		int dis = k * M*N;
		int off = i*N;
		Elem[off + j] += Elem[dis + off + j];
	}
	__syncthreads();
}

//how = TOBLOCK: basis transform to block-diagonal basis
//how = OFFBLOCK: basis transform to original basis
void chCol(Tensor* T, int row_bond_num, int how){
	assert((how == OFFBLOCK) || (how == TOBLOCK));
	assert(T->status & HAVEELEM);
	assert((T->sym)->mapNum);
	int M = (T->sym)->totRow;
	int N = (T->sym)->totCol;
	int mapNum = (T->sym)->mapNum;
	double *newElem;
	assert(cudaMalloc((void**)&newElem, sizeof(double) * T->elemNum * mapNum) == cudaSuccess);
	assert(cudaMemset(newElem, 0, sizeof(double) * T->elemNum * mapNum) == cudaSuccess);
	int *table;
	double *factor;
	assert(cudaMalloc((void**)&table, sizeof(int) * N * mapNum) == cudaSuccess);
	assert(cudaMalloc((void**)&factor, sizeof(double) * N * mapNum) == cudaSuccess);
	if(how == OFFBLOCK){
		assert(T->status & CBLOCK);
		assert(cudaMemcpy(table, (T->sym)->invTable, sizeof(int) * N * mapNum, cudaMemcpyHostToDevice) == cudaSuccess);
		if (mapNum > 1)
			assert(cudaMemcpy(factor, (T->sym)->invFactor, sizeof(double) * N * mapNum, cudaMemcpyHostToDevice) == cudaSuccess);
	}
	else{//(how == TOBLOCK)
		assert((T->status & CBLOCK) == 0);
		assert(cudaMemcpy(table, (T->sym)->table, sizeof(int) * N * mapNum, cudaMemcpyHostToDevice) == cudaSuccess);
		if (mapNum > 1)
			assert(cudaMemcpy(factor, (T->sym)->factor, sizeof(double) * N * mapNum, cudaMemcpyHostToDevice) == cudaSuccess);
	}
	dim3 block_size;
	block_size.x = 16;//I don't know what value is better here???????????????????????????????????????????????????????????????
	block_size.y = 16;

	dim3 grid_size;
	grid_size.x = (int)(N*mapNum + block_size.x - 1) / block_size.x;
	grid_size.y = (int)(M + block_size.y - 1) / block_size.y;

	if(mapNum == 1)
		ColMap1<<<grid_size , block_size>>>(table, M, N, T->elem, newElem);
	else{
		ColMap2<<<grid_size , block_size>>>(table, factor, M, N, mapNum, T->elem, newElem);
		grid_size.x = (int)(N + block_size.x -1) / block_size.x;
		for (int k = 1; k<mapNum; k++)
			ColCollect<<<grid_size , block_size>>>(M, N, k, newElem);
	}
	T->status ^= CBLOCK;
	addElem(T, newElem, DtoD);
	cudaFree(newElem);
	cudaFree(table);
	cudaFree(factor);
}

void chRow(Tensor* T, int row_bond_num, int how){
	assert((how == OFFBLOCK) || (how == TOBLOCK));
	assert(T->status & HAVEELEM);
	assert((T->sym)->mapNum);
	double *newElem;
	assert(cudaMalloc((void**)&newElem, sizeof(double) * T->elemNum) == cudaSuccess);
	assert(cudaMemset(newElem, 0, sizeof(double) * T->elemNum ) == cudaSuccess);
	int M = (T->sym)->totRow;
	int N = (T->sym)->totCol;
	int mapNum = (T->sym)->mapNum;
	int *table;
	double* factor;
	if(how == OFFBLOCK){
		assert(T->status & RBLOCK);
		table = (T->sym)->invTable;
		factor = (T->sym)->invFactor;
	}
	else{//(how == TOBLOCK)
		assert((T->status & RBLOCK) == 0);
		table = (T->sym)->table;
		factor = (T->sym)->factor;
	}
	if(mapNum == 1)
		for(int i = 0; i < M; i++)
			assert(cudaMemcpy(newElem + table[i] * N, T->elem + i * N, sizeof(double) * N, cudaMemcpyDeviceToDevice) == cudaSuccess);
	else
		for(int i = 0; i < M; i++){
			int inc = 1;
			for(int k = 0; k < mapNum; k++)
				cublasDaxpy(N, factor[i * mapNum + k], T->elem + i * N, inc, newElem + (table[i * mapNum + k] % M) * N, inc);
		}
	T->status ^= RBLOCK;
	addElem(T, newElem, DtoD);
	cudaFree(newElem);
}

void Block_Check(Tensor* T, int row_bond_number){
    #ifdef BLOCK_CHECK
	assert(T->status & HAVESYM);
	assert(T->status & ONBLOCK);
	double* Ta = (double*)malloc(T->elemNum * sizeof(double));
	assert(cudaMemcpy(Ta, T->elem, sizeof(double) *  T->elemNum, cudaMemcpyDeviceToHost) == cudaSuccess);
		                
	double error = (1.0E-10);
	for(int g = 0; g < (T->sym)->group.size(); g++){
		int row_start = (T->sym)->group[g].row_begin;
		int row_end   = row_start + (T->sym)->group[g].row;
		double temp;  
		for (int i = row_start; i< row_end; i++){
			for (int j = 0; j<(T->sym)->group[g].col_begin; j++){
				temp = Ta[ i * ((T->sym)->totCol) + j ];
				if(temp < 0)
					temp = -temp;
				if (temp > error){
			   		printf("1. Group[%d] : elem[%d][%d] = %e\n", g, i, j, Ta[ i * ((T->sym)->totCol) + j ]);
					print_tensor(T, row_bond_number);
				}
                assert(temp < error);
          	}
			for (int j = (T->sym)->group[g].col_begin + (T->sym)->group[g].col; j<(T->sym)->totCol; j++){
				temp = Ta[ i * ((T->sym)->totCol) + j ];
				if (temp < 0)
					temp = -temp;
				if (temp > error){
				    printf("2. Group[%d] : elem[%d][%d] = %e\n", g, i, j, Ta[ i * ((T->sym)->totCol) + j ]);
					print_tensor(T, row_bond_number);
				}
				assert(temp < error);
			}
		}
	}
	free(Ta);
	#endif
}

void Unitary_Check(Tensor* T, int row_bond_num){
	#ifdef UNITARY_CHECK
	assert(T->status & HAVEELEM);
	double error = 1.0E-10;
	int M = 1;
	int N = 1;
	for(int i = 0; i<row_bond_num; i++)
		M *= T->bondDim[i];
	for(int i = row_bond_num; i<T->bondNum; i++)
	    N *= T->bondDim[i];

	double* Result = (double*)malloc(M*M*sizeof(double));
	double *ResultDev;
	assert(cudaMalloc((void**)&ResultDev, sizeof(double) * M*M) == cudaSuccess);
	Matrix_Product(T->elem, 0, T->elem, 1, M, N, M, ResultDev);
	assert(cudaMemcpy(Result, ResultDev, sizeof(double) * M*M, cudaMemcpyDeviceToHost) == cudaSuccess);
    for (int i=0; i<M; i++)
    	Result[i*M + i] -= 1;
	for (int i=0; i<M*M; i++)
		if (Result[i]<0)
        	Result[i] *= -1;

    for (int i=0; i<M*M; i++){
	    if (Result[i] >= error){
	    	print_tensor(T, row_bond_num);
			printf("T * T dagger = :\n\n");
			print_matrix(Result,M,M);
			printf("elem[%d][%d] = %e\n", i/M, i%M, Result[i]);
		}
		assert (Result[i] < error);
	}
	cudaFree(ResultDev);
	free(Result);
	#endif
}

void Matrix_Product(double* A, int At, double* B, int Bt, int m, int k, int n, double* C){
	int M = n;
	int K = k;
	int N = m;
	double* a = B;
	double* b = A;
	int at = Bt;
	int bt = At;
	char transa;
	char transb;
	double alpha = 1;
	double beta = 0;
	int lda,ldb;
	int ldc = M;
	//at == 0 means no transpose, but we have to transpose it for mkl
	if (at == 0){
		transa = 'N';
		lda = M;
	}
	else if (at == 1){
	    transa = 'T';
		lda = K;
	}
	else
		printf("at must be 0 or 1!!!!!!\n");

	if (bt == 0){
		transb = 'N';
		ldb = K;
	}
	else if (bt == 1){
		transb = 'T';
		ldb = N;
	}
	else
	printf("bt must be 0 or 1!!!!!!\n");
	
	cublasDgemm(transa, transb, M, N, K, alpha, a, lda, b, ldb, beta, C, ldc);
}

void print_matrix(double* Matrix, int m, int n){
	printf("\n\n%dx%d Matrix:\n\n", m, n);
	for (int i=0; i<m; i++){
		for (int j=0; j<n; j++)
			printf("%6.3f", Matrix[i*n + j]);
		printf("\n");
	}
	printf("\n\n");
}

void Substract_Identity(double* Matrix, int size, double factor){
	for (int i=0; i<size; i++)
		Matrix[i*size + i] -= factor;
}

double Trace(double* A, double* B, int N){//This function only take trace of two square matrix
	double* C;
	assert(cudaMalloc((void**)&C, N*N*sizeof(double)) == cudaSuccess);
	Matrix_Product(A, 0, B, 0, N, N, N, C);
	double* CHost = (double*)malloc(N*N*sizeof(double));
	assert(cudaMemcpy(CHost, C, N*N*sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess);
	double TrC = 0;
	for (int i=0; i<N; i++)
		TrC += CHost[(N+1) * i];
	cudaFree(C);
	free(CHost);
	return TrC;
}

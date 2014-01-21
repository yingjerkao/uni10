#include <boost/random.hpp>
#include "myTools.h"
#include <stdio.h>
using namespace boost;
mt19937 lapack_rng(777);
uniform_01<mt19937> lapack_uni01_sampler(lapack_rng);

size_t MEM_USAGE = 0;
size_t GPU_MEM_USAGE = 0;
const size_t GPU_MEM_MAX = 1<<30;
const int THREADMAX = 256;
const int BLOCKMAX = 65535;


void* myMalloc(size_t memsize, int& status){
	void* ptr = NULL;
	if((USEGPU || (status & ONGPU)) && (GPU_MEM_USAGE + memsize <= GPU_MEM_MAX)){
		cudaError_t cuflag = cudaMalloc(&ptr, memsize);
		assert(cuflag == cudaSuccess);
		//printf("GPU: %d used, allocate %d bytes, %x\n", GPU_MEM_USAGE, memsize, ptr);
		GPU_MEM_USAGE += memsize;
		status |= ONGPU;
	}else{
		//printf("CPU: %d used, allocate %d bytes, %x\n", MEM_USAGE, memsize, ptr);
		ptr = malloc(memsize);
		assert(ptr != NULL);
		MEM_USAGE += memsize;
		status &= ~ONGPU;
	}
	return ptr;
}

void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st){
	cudaError_t cuflag;
	if((des_st & ONGPU)){
		if(src_st & ONGPU){
			cuflag = cudaMemcpy(des, src, memsize, cudaMemcpyDeviceToDevice);
			//printf("memcpy %x to %x D2D\n", src, des);
		}
		else{
			cuflag = cudaMemcpy(des, src, memsize, cudaMemcpyHostToDevice);
			//printf("memcpy %x to %x H2D\n", src, des);
		}
		assert(cuflag == cudaSuccess);
	}else{
		if(src_st & ONGPU){
			cuflag = cudaMemcpy(des, src, memsize, cudaMemcpyDeviceToHost);
			//printf("memcpy %x to %x D2H\n", src, des);
			assert(cuflag == cudaSuccess);
		}
		else{
			memcpy(des, src, memsize);
			//printf("memcpy H2H\n");
		}
	}
	return des;
}

void myFree(void* ptr, size_t memsize, int status){
	cudaError_t cuflag;
	assert(ptr != NULL);
	if(status & ONGPU){
		//printf("FREE(%x) %d from GPU, %d used\n", ptr, memsize, GPU_MEM_USAGE);
		cuflag = cudaFree(ptr);
		assert(cuflag == cudaSuccess);
		GPU_MEM_USAGE -= memsize;
	}else{
		//printf("FREE %d from CPU, %d used\n", memsize, MEM_USAGE);
		free(ptr);
		MEM_USAGE -= memsize;
	}
	ptr = NULL;
}
void membzero(void* ptr, size_t memsize, int status){
	if(status & ONGPU)
		cudaMemset(ptr, 0, memsize);
	else
		memset(ptr, 0, memsize);
}

void randomNums(double* elem, int N, int status){
	for(int i = 0; i < N; i++)
		elem[i] = lapack_uni01_sampler();
}

void myTranspose(double* A, int M, int N, double* AT, int status){
	for(int i = 0; i < M; i++)
		for(int j = 0; j < N; j++)
			AT[j * M + i] = A[i * N + j];
}

void myEye(double* elem, int M, int N, int status){
	int min;
	if(M < N)	min = M;
	else		min = N;
	memset(elem, 0, M * N * sizeof(double));
	for(int i = 0; i < min; i++)
		elem[i * N + i] = 1;
}

#define INTSZ 4
__device__ int Bsearch(int* arr, int len, int enc, int tar){
	int low = 0;
	int mid = 0;
	while(low < len){
		mid = low + ((len - low) / 2);
		if(tar < arr[mid * enc])
			len = mid;
		else
			low = mid;
		if(len - low <= 1){
			return low;
		}
	}
	return low;

}
typedef struct{
	unsigned int offset;
	int Cnum;
}_gpuBlock;

typedef struct{
	int bidx;
	int off;
	int dim;
}_gpuRQinfo;

typedef struct{
	int off;
	int dim;
}_gpuCQinfo;

typedef struct{
	int Qdeg;
	int prtF;
}_gpuQinfo;

typedef struct{
	int b1;
	int b2;
}_gpuSwap;

__global__ void reshape(double* des, int des_meta_num, int* des_meta, double* src, int src_meta_num, int* src_meta, int allocThread, int rsp_info_num, int* rsp_info){
	extern __shared__ int meta[];
	int idx_s = threadIdx.x;
	int metaoff = src_meta_num + des_meta_num;
	int metabase = src_meta_num + des_meta_num + rsp_info_num;
	while(idx_s < metabase){
		if(idx_s < src_meta_num)
			meta[idx_s] = src_meta[idx_s];
		else if(idx_s < metaoff)
			meta[idx_s] = des_meta[idx_s - src_meta_num];
		else
			meta[idx_s] = rsp_info[idx_s - metaoff];
		idx_s += blockDim.x;
	}
	__syncthreads();
	idx_s = blockIdx.y * BLOCKMAX * THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	if(idx_s < allocThread){
		int bondNum = meta[0];
		int enc_ot = Bsearch(&meta[8], meta[6], 2, idx_s);
		int Qin_off = meta[7 + enc_ot * 2];
		int idx_s1 = (idx_s - meta[8 + enc_ot * 2]) * ELEM_PER_THREAD;
		int Qin_RQoff = Qin_off / meta[3];	// Qidx / CQdim
		int Qin_CQoff = Qin_off % meta[3];	// Qidx % CQdim
		int Qot_off, Qot_RQoff, Qot_CQoff;
		int sRQdim = meta[2];
		int sCQdim = meta[4];
		metaoff = 7 + 2 * meta[6];
		Qin_RQoff = Bsearch(&meta[metaoff], sRQdim, 1, Qin_RQoff);
		metaoff += sRQdim;
		Qin_CQoff = Bsearch(&meta[metaoff], sCQdim, 1, Qin_CQoff);
		metaoff += sCQdim;
		_gpuBlock* blocks = (_gpuBlock*)&meta[metaoff];
		metaoff += 2 * meta[5];
		_gpuRQinfo* RQinfo = (_gpuRQinfo*)&meta[metaoff];
		metaoff += 3 * sRQdim;
		_gpuCQinfo* CQinfo = (_gpuCQinfo*)&meta[metaoff];
		metaoff += 2 * sCQdim;
		int sBin_cDim = CQinfo[Qin_CQoff].dim;

		int* Qin_Qdims = (int*)&meta[metaoff];	//Qin_Qdims
		metaoff += bondNum;
		int* bondOff = (int*)&meta[metaoff];	//bondOff
		metaoff += bondNum;
		_gpuQinfo* Qnums = (_gpuQinfo*)&meta[metaoff];
		metaoff = src_meta_num + des_meta_num;
		int* Qot_acc = (int*)&meta[metaoff];	//Qot_acc
		metaoff += bondNum;
		int* rsp_outin = (int*)&meta[metaoff];	//rsp_outin
		metaoff += bondNum;
		sRQdim = meta[metaoff];	//number of swaps
				
		int* sBin_sBdims = (int*)&meta[metabase + threadIdx.x * bondNum];
		int* sBot_acc = (int*)&meta[metabase + blockDim.x * bondNum + threadIdx.x * bondNum];
		int* sBin_idxs = (int*)&meta[metabase + 2 * blockDim.x * bondNum + threadIdx.x * bondNum];

		sCQdim = RQinfo[Qin_RQoff].bidx;	//bidx
		int Bin_cDim = blocks[sCQdim].Cnum;	//Bin_cDim
		int Ein_off = blocks[sCQdim].offset + (RQinfo[Qin_RQoff].off * Bin_cDim) + CQinfo[Qin_CQoff].off;
			
		Qot_off = Qin_off;	//Qot_off used as tmp
		for(int b = bondNum - 1; b >= 0; b--){
			sBin_sBdims[b] = Qot_off % Qin_Qdims[b];	//sBin_sBdims = Qin_idxs
			Qot_off /= Qin_Qdims[b];
		}
		Qot_off = 0;
		for(int b = 0; b < bondNum; b++)
			Qot_off += sBin_sBdims[rsp_outin[b]] * Qot_acc[b];

#ifdef FERMIONIC
		_gpuSwap* swaps = (_gpuSwap*)&meta[metaoff + 1];
		int sign = 0;
		for(int i = 0; i < sRQdim; i++)	//sRQdim = number of swaps
			sign ^= (Qnums[bondOff[swaps[i].b1] + sBin_sBdims[swaps[i].b1]].prtF & Qnums[bondOff[swaps[i].b2] + sBin_sBdims[swaps[i].b2]].prtF);
		sign = sign ? -1: 1;
#endif

		Qin_off = sBin_cDim * RQinfo[Qin_RQoff].dim;	//number of element per Qidx set

		sRQdim = meta[src_meta_num + 2];
		sCQdim = meta[src_meta_num + 4];
		Qot_RQoff = Qot_off / meta[src_meta_num + 3];	// Qidx / CQdim
		Qot_CQoff = Qot_off % meta[src_meta_num + 3];	// Qidx % CQdim
		metaoff = src_meta_num + 7 + 2 * meta[src_meta_num + 6];
		Qot_RQoff = Bsearch(&meta[metaoff], sRQdim, 1, Qot_RQoff);
		metaoff += sRQdim;
		Qot_CQoff = Bsearch(&meta[metaoff], sCQdim, 1, Qot_CQoff);
		metaoff += sCQdim;
		blocks = (_gpuBlock*)&meta[metaoff];
		metaoff += 2 * meta[src_meta_num + 5];
		RQinfo = (_gpuRQinfo*)&meta[metaoff];
		metaoff += 3 * sRQdim;
		CQinfo = (_gpuCQinfo*)&meta[metaoff];
		int sBot_cDim = CQinfo[Qot_CQoff].dim;

		sCQdim = RQinfo[Qot_RQoff].bidx;	//bidx
		int Bot_cDim = blocks[sCQdim].Cnum;
		int Eot_off = blocks[sCQdim].offset + (RQinfo[Qot_RQoff].off * Bot_cDim) + CQinfo[Qot_CQoff].off;
		
		sBot_acc[rsp_outin[bondNum - 1]] = 1;
		for(int b = bondNum	- 1; b > 0; b--){
			sBot_acc[rsp_outin[b-1]] = sBot_acc[rsp_outin[b]] * Qnums[bondOff[rsp_outin[b]] + sBin_sBdims[rsp_outin[b]]].Qdeg; 
		}
		for(int b = 0; b < bondNum; b++){
			sBin_sBdims[b] = Qnums[bondOff[b] + sBin_sBdims[b]].Qdeg;
		}
		enc_ot = 0;
		int cnt = idx_s1;
		for(int b = bondNum	- 1; b >= 0; b--){
			sBin_idxs[b] = cnt % sBin_sBdims[b];
			cnt /= sBin_sBdims[b];
			enc_ot += sBin_idxs[b] * sBot_acc[b];
		}
		Qin_RQoff = idx_s1 / sBin_cDim;	//Qin_RQoff = sBin_r
		Qin_CQoff = idx_s1 % sBin_cDim;	//Qin_RQoff = sBin_c
		for(cnt = 0; cnt < ELEM_PER_THREAD; cnt++){
			Qot_RQoff = enc_ot / sBot_cDim;	//Qot_RQoff = sBot_r
			Qot_CQoff = enc_ot % sBot_cDim;	//Qot_CQoff = sBot_c
#ifdef FERMIONIC
			des[Eot_off + (Qot_RQoff * Bot_cDim) + Qot_CQoff] = sign * src[Ein_off + (Qin_RQoff * Bin_cDim) + Qin_CQoff];
#else
			des[Eot_off + (Qot_RQoff * Bot_cDim) + Qot_CQoff] = src[Ein_off + (Qin_RQoff * Bin_cDim) + Qin_CQoff];
#endif
			idx_s1++;
			if(idx_s1 >= Qin_off)
				break;
			Qin_CQoff++;
			if(Qin_CQoff % sBin_cDim == 0){
				Qin_CQoff = 0;
				Qin_RQoff++;
			}
			for(int bend = bondNum - 1; bend >= 0; bend--){
				sBin_idxs[bend]++;
				if(sBin_idxs[bend] < sBin_sBdims[bend]){
					enc_ot += sBot_acc[bend];
					break;
				}
				else{
					enc_ot -= sBot_acc[bend] * (sBin_idxs[bend] - 1);
					sBin_idxs[bend] = 0;
				}
			}
		}
	}

}
void reshapeElem(int bondNum, double* des, int des_metasz, int* des_meta, double* src, int src_metasz, int* src_meta, int allocThread, int rsp_infosz, int* rsp_info){
	int blockSize = THREADMAX;
	int blockNum = (allocThread + THREADMAX - 1) / THREADMAX;
	dim3 gridSize(blockNum % BLOCKMAX, (blockNum + BLOCKMAX - 1) / BLOCKMAX);
	//printf("allocThread: %d\n", allocThread);
	int intsz = sizeof(int);
	reshape<<<gridSize, blockSize, src_metasz + des_metasz + rsp_infosz + 3 * bondNum * THREADMAX * intsz>>>(des, des_metasz / intsz, des_meta, src, src_metasz / intsz, src_meta, allocThread, rsp_infosz / intsz, rsp_info);
}

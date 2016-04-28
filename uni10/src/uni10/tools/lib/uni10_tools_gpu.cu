#include <uni10/tools/uni10_tools.h>
#include <string.h>

namespace uni10{
  
const size_t GPU_MEM_MAX = UNI10_GPU_GLOBAL_MEM * 2 / 3;

void* elemAlloc(size_t memsize, bool& ongpu){
  void* ptr = NULL;
  if(GPU_MEM_USAGE + memsize <= GPU_MEM_MAX){
    cudaError_t cuflag = cudaMalloc(&ptr, memsize);
    assert(cuflag == cudaSuccess);
    GPU_MEM_USAGE += memsize;
    ongpu = true;
  }else{
    ptr = malloc(memsize);
    assert(ptr != NULL);
    MEM_USAGE += memsize;
    ongpu = false;
  }
  //printf("ongpu = %d, GPU_MEM_USAGE = %u, allocate %u\n", ongpu, GPU_MEM_USAGE, memsize);
  return ptr;
}

void* elemAllocForce(size_t memsize, bool ongpu){
  void* ptr = NULL;
  if(ongpu){
    cudaError_t cuflag = cudaMalloc(&ptr, memsize);
    assert(cuflag == cudaSuccess);
    GPU_MEM_USAGE += memsize;
  }
  else{
    ptr = malloc(memsize);
    assert(ptr != NULL);
    MEM_USAGE += memsize;
  }
  return ptr;
}

void* elemCopy(void* des, const void* src, size_t memsize, bool des_ongpu, bool src_ongpu){
	cudaError_t cuflag;
	if((des_ongpu)){
		if(src_ongpu){
			cuflag = cudaMemcpy(des, src, memsize, cudaMemcpyDeviceToDevice);
		}
		else{
			//printf("mvGPU\n");
			cuflag = cudaMemcpy(des, src, memsize, cudaMemcpyHostToDevice);
		}
		assert(cuflag == cudaSuccess);
	}else{
		if(src_ongpu){
			//printf("mvCPU\n");
			cuflag = cudaMemcpy(des, src, memsize, cudaMemcpyDeviceToHost);
			assert(cuflag == cudaSuccess);
		}
		else{
			memcpy(des, src, memsize);
		}
	}
	return des;
}

void elemFree(void* ptr, size_t memsize, bool ongpu){
	cudaError_t cuflag;
	assert(ptr != NULL);
	if(ongpu){
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
void elemBzero(void* ptr, size_t memsize, bool ongpu){
	if(ongpu)
		cudaMemset(ptr, 0, memsize);
	else
		memset(ptr, 0, memsize);
}

__global__ void gpu_rand(double* elem, size_t N){
	size_t idx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int r = (1664525 * ((1664525 * idx + 1013904223) % UINT_MAX) + 1013904223) % UINT_MAX;
	if(idx < N)
		elem[idx] = double(r) / UINT_MAX;
}

void elemRand(double* elem, size_t N, bool ongpu){
	if(ongpu){
		size_t blockNum = (N + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
		dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
		gpu_rand<<<gridSize, UNI10_THREADMAX>>>(elem, N);
	}
	else{
		for(size_t i = 0; i < N; i++)
			elem[i] = ((double)rand()) / RAND_MAX; //lapack_uni01_sampler();
	}
}

__global__ void _setDiag(double* elem, double* diag_elem, size_t M, size_t N, size_t diag_N){
	size_t idx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	if(idx < diag_N && idx < M && idx < N)
		elem[idx * N + idx] = diag_elem[idx];
}

__global__ void _getDiag(double* elem, double* diag_elem, size_t M, size_t N, size_t diag_N){
	size_t idx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	if(idx < diag_N && idx < M && idx < N)
		diag_elem[idx] = elem[idx * N + idx];
}

void setDiag(double* elem, double* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu){
	if((ongpu)){
		size_t blockNum = (N + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
		dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
		if(diag_ongpu){
			_setDiag<<<gridSize, UNI10_THREADMAX>>>(elem, diag_elem, M, N, diag_N);
		}
		else{
			size_t memsize = diag_N * sizeof(double);
			double* src_elem;
			cudaError_t cuflag = cudaMalloc(&src_elem, memsize);
			assert(cuflag == cudaSuccess);
			cuflag = cudaMemcpy(src_elem, diag_elem, memsize, cudaMemcpyHostToDevice);
			//printf("mvGPU");
			assert(cuflag == cudaSuccess);
			_setDiag<<<gridSize, UNI10_THREADMAX>>>(elem, src_elem, M, N, diag_N);
			cudaFree(src_elem);
		}
	}else{
		double* src_elem;
		if(diag_ongpu){
			size_t memsize = diag_N * sizeof(double);
			src_elem = (double*) malloc(memsize);
			cudaError_t cuflag = cudaMemcpy(src_elem, diag_elem, memsize, cudaMemcpyDeviceToHost);
			//printf("mvCPU");
			assert(cuflag == cudaSuccess);
		}
		else
			src_elem = diag_elem;
		int min = M < N ? M : N;
		min = min < diag_N ? min : diag_N;
		for(size_t i = 0; i < min; i++)
			elem[i * N + i] = src_elem[i];
	}
}

void getDiag(double* elem, double* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu){
	if((ongpu)){
		size_t blockNum = (N + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
		dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
		if(diag_ongpu){
			_getDiag<<<gridSize, UNI10_THREADMAX>>>(elem, diag_elem, M, N, diag_N);
		}
		else{
			size_t memsize = diag_N * sizeof(double);
			double* tmp_elem;
			cudaError_t cuflag = cudaMalloc(&tmp_elem, memsize);
			assert(cuflag == cudaSuccess);
			_getDiag<<<gridSize, UNI10_THREADMAX>>>(elem, tmp_elem, M, N, diag_N);
			cuflag = cudaMemcpy(diag_elem, tmp_elem, memsize, cudaMemcpyHostToDevice);
			//printf("mvGPU");
			assert(cuflag == cudaSuccess);
			cudaFree(tmp_elem);
		}
	}else{
		double* tmp_elem;
		size_t memsize = diag_N * sizeof(double);
		if(diag_ongpu)
			tmp_elem = (double*)malloc(memsize);
		else
			tmp_elem = diag_elem;
		int min = M < N ? M : N;
		min = min < diag_N ? min : diag_N;
		for(size_t i = 0; i < min; i++)
			tmp_elem[i] = elem[i * N + i];
		if(diag_ongpu){
			cudaError_t cuflag = cudaMemcpy(diag_elem, tmp_elem, memsize, cudaMemcpyDeviceToHost);
			//printf("mvCPU");
			assert(cuflag == cudaSuccess);
		}
	}
}

void* mvGPU(void* elem, size_t memsize, bool& ongpu){
	if(!ongpu)
		if(GPU_MEM_USAGE + memsize <= GPU_MEM_MAX){
			void* newElem = elemAlloc(memsize, ongpu);
			elemCopy(newElem, elem, memsize, ongpu, false);
			elemFree(elem, memsize, false);
			elem = newElem;
		}
	return elem;
}

void* mvCPU(void* elem, size_t memsize, bool& ongpu){
	if(ongpu){
		double *newElem = (double*)malloc(memsize);
		elemCopy(newElem, elem, memsize, false, true);
		elemFree(elem, memsize, true);
		MEM_USAGE += memsize;
		ongpu = false;
		elem = newElem;
	}
	return elem;
}

void syncMem(void** elemA, void** elemB, size_t memsizeA, size_t memsizeB, bool& ongpuA, bool& ongpuB){
	if((!ongpuA) || (!ongpuB)){
		size_t memsize = 0;
		if(!ongpuA)
			memsize += memsizeA;
		if(!ongpuB)
			memsize += memsizeB;
		if(GPU_MEM_USAGE + memsize <= GPU_MEM_MAX){
			if(!ongpuA)
				*elemA = mvGPU(*elemA, memsizeA, ongpuA);
			if(!ongpuB)
				*elemB = mvGPU(*elemB, memsizeB, ongpuB);
		}
		else{
			if(ongpuA)
				*elemA = mvCPU(*elemA, memsizeA, ongpuA);
			if(ongpuB)
				*elemB = mvCPU(*elemB, memsizeB, ongpuB);
		}
	}
}
void shrinkWithoutFree(size_t memsize, bool ongpu){
	if(ongpu)
		GPU_MEM_USAGE -= memsize;
	else
		MEM_USAGE -= memsize;
}

__global__ void _reshapeElem(double* oldElem, int bondNum, size_t elemNum, size_t* offset, double* newElem){
	size_t oldIdx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	size_t idx = oldIdx;
	size_t newIdx = 0;
	if(idx < elemNum){
		for(int i = 0; i < bondNum; i++){
			newIdx += (idx/offset[i]) * offset[bondNum + i];
			idx = idx % offset[i];
		}
		newElem[newIdx] = oldElem[oldIdx];
	}
}

void reshapeElem(double* oldElem, int bondNum, size_t elemNum, size_t* offset, double* newElem){
	size_t* D_offset;
	assert(cudaMalloc((void**)&D_offset, 2 * sizeof(size_t) * bondNum) == cudaSuccess);
	assert(cudaMemcpy(D_offset, offset, 2 * sizeof(size_t) * bondNum, cudaMemcpyHostToDevice) == cudaSuccess);
	size_t blockNum = (elemNum + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
	dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
	_reshapeElem<<<gridSize, UNI10_THREADMAX>>>(oldElem, bondNum, elemNum, D_offset, newElem);
}

double getElemAt(size_t idx, double* elem, bool ongpu){
	if(ongpu){
		double val;
		assert(cudaMemcpy(&val, &(elem[idx]), sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess);
		return val;
	}
	else
		return elem[idx];
}

void setElemAt(size_t idx, double val, double* elem, bool ongpu){
	if(ongpu){
		assert(cudaMemcpy(&(elem[idx]), &val, sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess);
	}
	else
		elem[idx] = val;
}

double  elemMax(double* elem, size_t elemNum, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

double  elemAbsMax(double* elem, size_t elemNum, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

std::complex<double> getElemAt(size_t idx, std::complex<double>* elem, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void setElemAt(size_t idx, std::complex<double> val, std::complex<double>* elem, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void elemRand(std::complex<double>* elem, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void elemCast(std::complex<double>* des, double* src, size_t N, bool des_ongpu, bool src_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void elemCast(double *des, std::complex<double> *src, size_t N, bool des_ongpu, bool src_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void setDiag(std::complex<double>* elem, std::complex<double>* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu){
	
  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void getDiag(std::complex<double>* elem, std::complex<double>* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void reshapeElem(std::complex<double>* oldElem, int bondNum, size_t elemNum, size_t* offset, std::complex<double>* newElem){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

};	/* namespace uni10 */

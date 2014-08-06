#include <stdio.h>
void printDevProp(FILE *fp)
{
	int dev;
	cudaGetDevice(&dev);
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, dev);
    fprintf(fp, "\n\n***Device Property***\n");
    fprintf(fp, "Device Number:                 %d\n",  dev);
    fprintf(fp, "Major revision number:         %d\n",  devProp.major);
    fprintf(fp, "Minor revision number:         %d\n",  devProp.minor);
    fprintf(fp, "Name:                          %s\n",  devProp.name);
    fprintf(fp, "Total global memory:           %u\n",  devProp.totalGlobalMem);
    fprintf(fp, "Total shared memory per block: %u\n",  devProp.sharedMemPerBlock);
    fprintf(fp, "Total registers per block:     %d\n",  devProp.regsPerBlock);
    fprintf(fp, "Warp size:                     %d\n",  devProp.warpSize);
    fprintf(fp, "Maximum memory pitch:          %u\n",  devProp.memPitch);
    fprintf(fp, "Maximum threads per block:     %d\n",  devProp.maxThreadsPerBlock);
    for (int i = 0; i < 3; ++i)
    fprintf(fp, "Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
    for (int i = 0; i < 3; ++i)
    fprintf(fp, "Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
    fprintf(fp, "Clock rate:                    %d\n",  devProp.clockRate);
    fprintf(fp, "Total constant memory:         %u\n",  devProp.totalConstMem);
    fprintf(fp, "Texture alignment:             %u\n",  devProp.textureAlignment);
    fprintf(fp, "Concurrent copy and execution: %s\n",  (devProp.deviceOverlap ? "Yes" : "No"));
    fprintf(fp, "Number of multiprocessors:     %d\n",  devProp.multiProcessorCount);
    fprintf(fp, "Kernel execution timeout:      %s\n",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
    return;
}
int main(){
	FILE* fp = fopen("device.txt", "w");
	printDevProp(fp);
	return 0;
}

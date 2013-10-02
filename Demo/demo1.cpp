#include "network.h"

int main(int argc, char* argv[]){
	Tensor* T = (Tensor*)calloc(1, sizeof(Tensor));
	int bondNum = 3;
	int bondDim[] = {4, 2, 3};
	int label[] = {-1, -2, -3};
	initTensor(T, bondNum, bondDim);
	addLabel(T, label);
	print_tensor(T, 1);
}

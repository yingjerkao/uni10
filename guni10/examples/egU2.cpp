#include <iostream>
#include <uni10.hpp>

int main(){
	// Construct spin 1 Heisenberg model by reading in the tensor which is written out in example egU1
	uni10::UniTensor H_U1("egU1_H_U1");

	// Get the block of quantum number q0 as a matrix "block0"
	uni10::Qnum q0(0);
	uni10::Matrix block0 = H_U1.getBlock(q0);
	std::cout<<block0;
	// Randomly assign "block0" and put it back to H_U1
	block0.randomize();
	H_U1.putBlock(q0, block0);
	//std::cout<<H_U1;

	// Permute bonds by its labels, the default labels are [0 1 2 3]
	int permuted_label[] = {1, 2, 3, 0};
	// Permute bonds to which with label [1, 2, 3, 0] and leaving 1 bond as in-coming bonds.
	H_U1.permute(permuted_label, 1);
	//std::cout<<H_U1;

	// combine the two bonds with label 2 and 3
	std::vector<int> combined_label;
	combined_label.push_back(2);
	combined_label.push_back(3);	// combined_label = [2, 3]
	H_U1.combineBond(combined_label);
	std::cout<< H_U1;

	return 0;
}



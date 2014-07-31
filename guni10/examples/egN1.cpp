#include <iostream>
#include <uni10.hpp>

int main(){
	// Read in the tensor H_U1 which is written out in example egU1 and W, WT in example egU3
	uni10::UniTensor H_U1("egU1_H_U1");
	uni10::UniTensor W("egU3_W");
	uni10::UniTensor WT("egU3_WT");

	// Create network by reading in network file "egN1_network"
	uni10::Network net("egN1_network");
	// Put tensors to the Network net
	net.putTensor("H", &H_U1);
	net.putTensor("W", &W);
	net.putTensor("WT", &WT);

	// Perform contractions inside the tensor network
	std::cout<<net.launch();
	// Print out the network
	std::cout<<net;

	return 0;
}



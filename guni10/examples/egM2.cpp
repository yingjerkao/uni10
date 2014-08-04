#include <iostream>
#include <uni10.hpp>

int main(){
	uni10::Matrix M(4, 5);
	M.randomize();
	std::cout<<M;
	// carry out SVD
	std::vector<uni10::Matrix> rets = M.svd();
	
	// write matrice out to file
	rets[0].save("mat_U");
	rets[1].save("mat_Sigma");
	rets[2].save("mat_VT");

	uni10::Matrix U(rets[0].row(), rets[0].col(), rets[0].isDiag());
	uni10::Matrix S(rets[1].row(), rets[1].col(), rets[1].isDiag());
	uni10::Matrix VT(rets[2].row(), rets[2].col(), rets[2].isDiag());

	// read in the matrice we just write out
	U.load("mat_U");
	S.load("mat_Sigma");
	VT.load("mat_VT");
	std::cout<< S;
	std::cout<< U * S * VT;
}

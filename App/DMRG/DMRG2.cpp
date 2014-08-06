#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
const int M = 20;

UniTensor mirror(const UniTensor& T);

int main(){
	/*** Initialization ***/
	Qnum q0(0);
	Bond bdi(BD_IN, 2);
	Bond bdo(BD_OUT, 2);
	vector<Bond> bond4;
	bond4.push_back(bdi);
	bond4.push_back(bdi);
	bond4.push_back(bdo);
	bond4.push_back(bdo);

	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

	UniTensor H0(bond4);
	H0.addRawElem(H_elem);

	vector<Bond> bond2;
	bond2.push_back(bdi);
	bond2.push_back(bdo);

	double sig_z_elem[] = {0.5,    0,
					         0, -0.5};
	double sig_p_elem[] = {  0, 1.0,
		                     0,   0};
	double sig_m_elem[] = {  0, 0,
             		       1.0, 0};

	UniTensor sig_z(bond2);
	UniTensor sig_p(bond2);
	UniTensor sig_m(bond2);
    sig_z.addRawElem(sig_z_elem);
	sig_p.addRawElem(sig_p_elem);
	sig_m.addRawElem(sig_m_elem);
	UniTensor I2(bond2);
	I2.eye();
	UniTensor Sz = sig_z;
	UniTensor Sp = sig_p;
	UniTensor Sm = sig_m;
	UniTensor H(bond2);
	H.set_zero();
	int m = 2;
	/*** END initilization ***/

	for(int it = 0; it < 30; it++){
		/*** super block construction ***/
		UniTensor Hs = otimes(H, I2);	//H system
		Hs += otimes(Sz, sig_z) + 0.5 * otimes(Sp, sig_m) + 0.5 * otimes(Sm, sig_p);
		UniTensor Im2 = Hs;
		Im2.eye();
		UniTensor I2m = mirror(Im2);
		UniTensor He = mirror(Hs);	//H environment
		UniTensor Im = H;
		Im.eye();
		UniTensor SB = otimes(Hs, I2m);
		SB += otimes(Im2, He);
		SB += otimes(otimes(Im, H0), Im);
		/*** END construction ***/
		Matrix mat_SB = SB.getBlock(q0);
		vector<Matrix> rets = mat_SB.diagonalize();
		cout<<"N = "<< (it + 2) * 2 <<", m = " << m << setprecision(10) << ", E = " << rets[0][0] << ", e = " << rets[0][0] / ((it + 2) * 2) <<endl;
		Matrix mat_GS(1, mat_SB.row(), rets[1].elem());
		Matrix mat_GST = mat_GS;
		mat_GST.transpose();
		UniTensor Rho = SB;
		Rho.putBlock(q0, mat_GST * mat_GS);
		Rho.partialTrace(3, 7);
		Rho.partialTrace(2, 6);
		Matrix mat_Rho = Rho.getBlock(q0);
		rets = mat_Rho.diagonalize();
		m = rets[0].col();
		/*
		for(m = 0; m < rets[0].col(); m++)
			if(rets[0][rets[0].col() - m - 1] < 1E-7)
				break;
		*/
		m = m < M ? m : M;
		Matrix mat_TL(m, rets[1].col());	//renormalization matrix
		for(int r = 0; r < m; r++)
			memcpy(mat_TL.elem() + (r * rets[1].col()), rets[1].elem() + (rets[1].row() - r - 1) * rets[1].col(), rets[1].col() * sizeof(double));
		Matrix mat_TR = mat_TL;	
		mat_TR.transpose();
		bdi.assign(BD_IN, m);
		bdo.assign(BD_OUT, m);
		vector<Bond> bonds;
		bonds.push_back(bdi);
		bonds.push_back(bdo);
		H.assign(bonds);
		H.putBlock(q0, mat_TL * Hs.getBlock(q0) * mat_TR);
		Sz.assign(bonds);
		Sp.assign(bonds);
		Sm.assign(bonds);
		Sz.putBlock(q0, mat_TL * otimes(Im, sig_z).getBlock(q0) * mat_TR);
		Sp.putBlock(q0, mat_TL * otimes(Im, sig_p).getBlock(q0) * mat_TR);
		Sm.putBlock(q0, mat_TL * otimes(Im, sig_m).getBlock(q0) * mat_TR);
		//rets = H.getBlock(q0).diagonalize();
		//SP += tensorProduct(Sz, Sz) + tensorProduct(Sp, Sp) + tensorProduct(Sm, Sm);
	}
}

UniTensor mirror(const UniTensor& T){
	vector<int> labels = T.label();
	vector<int> mlabels(labels.size());
	for(int l = 0; l < labels.size(); l++){
		if(l < T.inBondNum())
			mlabels[l] = labels[T.inBondNum() - l - 1];
		else
			mlabels[l] = labels[T.inBondNum() + T.bondNum() - l - 1];
	}
	UniTensor mT = T;
	mT.permute(mlabels, T.inBondNum());
	return mT;
}

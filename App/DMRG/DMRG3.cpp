#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
const int M = 40;

UniTensor mirror(const UniTensor& T);

int main(){
	/*** Initialization ***/
	Qnum q0(0);
	Qnum q1(1);
	Qnum qm1(-1);
	Qnum q2(2);
	Qnum qm2(-2);
	Qnum q3(3);
	Qnum qm3(-3);
	vector<Qnum> qnums;
	qnums.push_back(q1); qnums.push_back(qm1);
	Bond bdi_sy(BD_IN, qnums);
	Bond bdo_sy(BD_OUT, qnums);
	//Bond bdi(BD_IN, 2);
	//Bond bdo(BD_OUT, 2);
	vector<Bond> bonds;
	bonds.push_back(bdi_sy);
	bonds.push_back(bdi_sy);
	bonds.push_back(bdo_sy);
	bonds.push_back(bdo_sy);

	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};
	double sig_z_elem[] = { 0.5,    0,
					          0, -0.5};
	double sig_p_elem[] = { 0, 1.0,
		                    0,   0};
	double sig_m_elem[] = {  0, 0,
             		       1.0, 0};

	UniTensor H0(bonds);
	H0.addRawElem(H_elem);

	Bond bdi(BD_IN, 2);
	Bond bdo(BD_OUT, 2);
	bonds.clear();
	bonds.push_back(bdi);
	bonds.push_back(bdo);
	UniTensor sig_z(bonds);
	UniTensor sig_p(bonds);
	UniTensor sig_m(bonds);
    sig_z.addRawElem(sig_z_elem);
	sig_p.addRawElem(sig_p_elem);
	sig_m.addRawElem(sig_m_elem);

	bonds.clear();
	bonds.push_back(bdi_sy);
	bonds.push_back(bdo_sy);
	UniTensor H(bonds);
	H.set_zero();
	UniTensor I2(bonds);
	I2.eye();
	UniTensor Sz = sig_z;
	UniTensor Sp = sig_p;
	UniTensor Sm = sig_m;
	/*** END initilization ***/

	for(int it = 0; it < 40; it++){
		UniTensor Im = H;
		Im.eye();
		UniTensor Hs_tmp = outer(H, I2);	//H system
		cout<<Hs_tmp;
		exit(0);
		Hs_tmp += outer(Sz, sig_z) + 0.5 * outer(Sp, sig_m) + 0.5 * outer(Sm, sig_p);
		UniTensor Hs = Hs_tmp;	
		UniTensor Im2 = Hs;
		Im2.eye();
		UniTensor I2m = mirror(Im2);
		UniTensor He = mirror(Hs);	//H environment
		UniTensor SB = outer(Hs, I2m);
		SB += outer(Im2, He);
		SB += outer(outer(Im, H0), Im);
		SB.printRawElem();
		Matrix mat_SB = SB.getBlock(q0);
		vector<Matrix> rets = mat_SB.diagonalize();
		cout<<"N = "<< (it + 2) * 2 << ", E = " << rets[0][0] / ((it + 2) * 2) <<endl;
		Matrix mat_GS(1, mat_SB.row(), rets[1].elem());
		Matrix mat_GST = mat_GS;
		mat_GST.transpose();
		UniTensor Rho = SB;
		Rho.putBlock(q0, mat_GST * mat_GS);
		Rho.partialTrace(3, 7);
		Rho.partialTrace(2, 6);
		Matrix mat_Rho = Rho.getBlock(q0);
		rets = mat_Rho.diagonalize();
		int m = rets[0].col();
		m = m < M ? m : M;
		Matrix mat_TL(m, rets[1].col());	//renormalization matrix
		for(int r = 0; r < m; r++)
			memcpy(mat_TL.elem() + (r * rets[1].col()), rets[1].elem() + (rets[1].row() - r - 1) * rets[1].col(), rets[1].col() * sizeof(double));
		Matrix mat_TR = mat_TL;	
		mat_TR.transpose();
		bdi.assign(BD_IN, m);
		bdo.assign(BD_OUT, m);
		bonds.clear();
		bonds.push_back(bdi);
		bonds.push_back(bdo);
		H.assign(bonds);
		H.putBlock(q0, mat_TL * Hs.getBlock(q0) * mat_TR);
		Sz.assign(bonds);
		Sp.assign(bonds);
		Sm.assign(bonds);
		Sz.putBlock(q0, mat_TL * outer(Im, sig_z).getBlock(q0) * mat_TR);
		Sp.putBlock(q0, mat_TL * outer(Im, sig_p).getBlock(q0) * mat_TR);
		Sm.putBlock(q0, mat_TL * outer(Im, sig_m).getBlock(q0) * mat_TR);
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

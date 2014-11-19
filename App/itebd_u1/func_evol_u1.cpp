#include <algorithm>
#include <mkl.h>
#include <uni10.hpp>

#include "func_net_u1.h"

//============================

/* return the exp(beta * OP) in UniTensor form */
uni10::UniTensor expBO(double beta, uni10::UniTensor op){

	std::vector<uni10::Bond> bd = op.bond();
	uni10::UniTensor eBO(bd);

	std::vector<uni10::Qnum> qnums = op.blockQnum();

	for (std::vector<int>::size_type i = 0; i != qnums.size(); i++)
		eBO.putBlock( qnums[i], uni10::takeExp(beta, op.getBlock(qnums[i])) );
	
	return eBO;
}

//============================

bool compare_function(double a, double b) {

	return (a > b);
}

//===========================================

void truncUpdate(uni10::UniTensor& wf, uni10::UniTensor& ga, uni10::UniTensor& gb, uni10::UniTensor& la, uni10::UniTensor& lb){

	/* retrieve bond dimensions from ga */
        std::vector<uni10::Bond> bd = ga.bond();
	int D = bd[0].dim();
        int d = bd[1].dim();
	std::vector<uni10::Bond> bd2 = gb.bond();

	if (wf.bondNum() != wf.inBondNum()*2)	wf.permute(wf.label(), wf.bondNum()/2);	// permute wf to square matrix


	/* retrieve qnums from wf */
	int qn = wf.blockNum();
	std::vector<uni10::Qnum> qnums = wf.blockQnum();


	/* do SVD on all blocks, save the results in vector */
	std::vector< std::vector<uni10::Matrix> > wf_svd;

	for (int i = 0; i < qn; ++i)
        	wf_svd.push_back( wf.getBlock(qnums[i]).svd() );


	/* sort all lambda elements */
	std::vector<double> lam_all;
	for (int i = 0; i < qn; ++i) {
		int block_dim = wf_svd[i][1].col();
        	for (int j = 0; j < block_dim; ++j) {
			lam_all.push_back( wf_svd[i][1].at(j, j) );
		}
	}
	std::sort(lam_all.begin(), lam_all.end(), compare_function);


	/* normalization factor */
	double scale;
	for (int i = 0; i < D; ++i) {
		scale += (lam_all[i] * lam_all[i]);
	}
	scale = std::sqrt(scale);


	/* assign new qnums for la and update its elements */	
	std::vector<uni10::Bond> bd_la;
	std::vector<uni10::Qnum> qnums_la;
	
	for (int i = 0; i < qn; ++i) {
		int block_dim = wf_svd[i][1].col();
        	for (int j = 0; j < block_dim; ++j) {
			if ( wf_svd[i][1].at(j, j) >= lam_all[D-1] ) {
				qnums_la.push_back( qnums[i] );
			}
			else {
				wf_svd[i][1].resize(j, j);
				j = block_dim;
			}
		}
		wf_svd[i][1] *= (1.0 / scale);
	}

	uni10::Bond la_in(uni10::BD_IN, qnums_la);
	uni10::Bond la_out(uni10::BD_OUT, qnums_la);
	
	bd_la.push_back( la_in );
	bd_la.push_back( la_out );

	uni10::UniTensor lam(bd_la);

	for (int i = 0; i < qn; ++i) {
		if ( wf_svd[i][1].col() > 0 ) {
			lam.putBlock( qnums[i], wf_svd[i][1] );
		}
	}

	/* assign new qnums for ga gb and update their elements */
	std::vector<uni10::Bond> bd_ga;
	bd_ga.push_back( bd[0] );
	bd_ga.push_back( bd[1] );
	bd_ga.push_back( bd_la[1] );
	uni10::UniTensor gamA(bd_ga);

	std::vector<uni10::Bond> bd_gb;
	bd_gb.push_back( bd_la[0] );
	bd_gb.push_back( bd2[1] );
	bd_gb.push_back( bd2[2] );
	uni10::UniTensor gamB(bd_gb);

        gamA.permute(gamA.label(), 2);
	gamB.permute(gamB.label(), 1);
	for (int i = 0; i < qn; ++i) {
		if ( wf_svd[i][1].col() > 0 ) {
			int dim = wf_svd[i][0].col();
			int block_dim = wf_svd[i][1].col();
			gamA.putBlock( qnums[i], wf_svd[i][0].resize( dim, block_dim ) );
			gamB.putBlock( qnums[i], wf_svd[i][2].resize( block_dim, dim ) );
		}
	}
        //gamA.permute(gamA.label(), 1);
        //gamB.permute(gamB.label(), 2);

	/* derive lb^-1 */
	uni10::UniTensor lb_inv( lb.bond() );
	for ( int i = 0; i < lb.blockNum(); ++i ) {
		uni10::Matrix tmp = lb.getBlock( lb.blockQnum()[i] );
		int dim = tmp.col();
		uni10::Matrix inv_lambda(dim, dim, true);

		for (int j = 0; j < dim; ++j)
			inv_lambda.at(j, j) = 1.0 / tmp.at(j, j);

		lb_inv.putBlock( lb.blockQnum()[i], inv_lambda );
	}
	
	gamA = netLG( gamA, lb_inv );
	gamB = netGL( gamB, lb_inv );

	ga = gamA;
	gb = gamB;
	la = lam;	

	lam_all.clear();
	wf_svd.clear();
}


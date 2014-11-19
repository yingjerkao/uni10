#include <uni10.hpp>
#include <string>

//======================

uni10::UniTensor initGamma(std::vector<uni10::Qnum> qnums_phys, int dim, bool AB){

	int physdim = (int) qnums_phys.size();	// for e.g., qnums_phys = (1, -1)

	std::vector<uni10::Qnum> qnumsA;
	std::vector<uni10::Qnum> qnumsB;

	for (int i = 0; i < physdim; ++i)
		for (int j = 0; j < physdim; ++j)
			qnumsB.push_back( qnums_phys[i] * qnums_phys[j] );	// (2, 0, 0, -2)

	for (int i = 0; i < physdim; ++i)
		for (int j = 0; j < physdim * physdim; ++j)
			qnumsA.push_back( qnums_phys[i] * qnumsB[j] );	// (3, 1, 1, 1, -1, -1, -1, -3)
	
	for (int i = 0; i < physdim * physdim; ++i)
			qnumsB.push_back( qnumsB[i] );	// (2, 0, 0, -2, 2, 0, 0, -2)

	int q_len = qnumsA.size();

	for (int i = 0; i < (dim - q_len); ++i)
		qnumsA.push_back( qnumsA[ (i % q_len) ] );

	for (int i = 0; i < (dim - q_len); ++i)
		qnumsB.push_back( qnumsB[ (i % q_len) ] );

	std::vector<uni10::Bond> bonds_gamma;

	if (AB) {
		uni10::Bond bd_in(uni10::BD_IN, qnumsB);
        	uni10::Bond bd_phys(uni10::BD_OUT, qnums_phys);
        	uni10::Bond bd_out(uni10::BD_OUT, qnumsA);
        	bonds_gamma.push_back(bd_in);
        	bonds_gamma.push_back(bd_phys);
        	bonds_gamma.push_back(bd_out);
	}
	else {
		uni10::Bond bd_in(uni10::BD_IN, qnumsA);
        	uni10::Bond bd_phys(uni10::BD_OUT, qnums_phys);
        	uni10::Bond bd_out(uni10::BD_OUT, qnumsB);
        	bonds_gamma.push_back(bd_in);
        	bonds_gamma.push_back(bd_phys);
        	bonds_gamma.push_back(bd_out);
	}

	uni10::UniTensor gamma(bonds_gamma);

	return gamma;
}

//=====================

uni10::UniTensor initLambda(std::vector<uni10::Qnum> qnums_phys, int dim, bool AB){
	
	int physdim = (int) qnums_phys.size();

	std::vector<uni10::Qnum> qnumsA;
	std::vector<uni10::Qnum> qnumsB;

	for (int i = 0; i < physdim; ++i)
		for (int j = 0; j < physdim; ++j)
			qnumsB.push_back( qnums_phys[i] * qnums_phys[j] );

	for (int i = 0; i < physdim; ++i)
		for (int j = 0; j < physdim * physdim; ++j)
			qnumsA.push_back( qnums_phys[i] * qnumsB[j] );
	
	for (int i = 0; i < physdim * physdim; ++i)
			qnumsB.push_back( qnumsB[i] );

	int q_len = qnumsA.size();

	for (int i = 0; i < (dim - q_len); ++i)
		qnumsA.push_back( qnumsA[ (i % q_len) ] );

	for (int i = 0; i < (dim - q_len); ++i)
		qnumsB.push_back( qnumsB[ (i % q_len) ] );

	std::vector<uni10::Bond> bonds_lambda;

	if (AB) {
		uni10::Bond bd_in(uni10::BD_IN, qnumsA);
        	uni10::Bond bd_out(uni10::BD_OUT, qnumsA);
       		bonds_lambda.push_back(bd_in);
        	bonds_lambda.push_back(bd_out);
	}
	else {
		uni10::Bond bd_in(uni10::BD_IN, qnumsB);
        	uni10::Bond bd_out(uni10::BD_OUT, qnumsB);
        	bonds_lambda.push_back(bd_in);
        	bonds_lambda.push_back(bd_out);
	}

        uni10::UniTensor lambda(bonds_lambda);

	return lambda;
}
/*
//=====================

uni10::UniTensor loadGamma(const std::string& fname, int dim){

	uni10::UniTensor ga(fname);

	std::vector<uni10::Bond> bd = ga.bond();
        int D = bd[0].dim();
	int d = bd[1].dim();

	if (dim == D){
		return ga;
	}
	else{
		uni10::Qnum q(0);
		uni10::UniTensor gamma = initGamma(0, dim, d);

		int lab_3b[3];
		lab_3b[0] = 0;
		lab_3b[1] = 100;
		lab_3b[2] = 1;
		ga.addLabel(lab_3b);
		gamma.addLabel(lab_3b);

		lab_3b[0] = 0;
                lab_3b[1] = 1;
                lab_3b[2] = 100;
		ga.permute(lab_3b, 1);
		gamma.permute(lab_3b, 1);

		uni10::Matrix tmp = gamma.getBlock(q);
        	//for (int i=0; i < ((dim<D)? dim : D); i++){
                	//for (int j=0; j < ((dim<D)? d*dim : d*D); j++){
		for (int i=0; i < dim; i++){
                        for (int j=0; j < d*dim ; j++){
				if (i<D && j<d*D)
                        		tmp.at(i, j) = ga.getBlock(q).at(i, j);
				else
					tmp.at(i, j) = 1e-64;
                	}
        	}
		gamma.putBlock(q, tmp);

		lab_3b[0] = 0;
                lab_3b[1] = 100;
                lab_3b[2] = 1;
                gamma.permute(lab_3b, 1);

		return gamma;
	}
}


//=====================

uni10::UniTensor loadLambda(const std::string& fname, int dim){

	uni10::UniTensor la(fname);

        std::vector<uni10::Bond> bd = la.bond();
        int D = bd[0].dim();

        if (dim == D){
                return la;
        }
	else{
                uni10::Qnum q(0);
                uni10::UniTensor lambda = initLamda(0, dim);

                uni10::Matrix tmp(dim, dim, true); 
                //for (int i=0; i < ((dim<D)? dim : D); i++)
		for (int i=0; i<dim; i++){
			if (i<D)
				tmp.at(i, i) = la.getBlock(q).at(i, i);
			else
				tmp.at(i, i) = 1e-64;
		}
                lambda.putBlock(q, tmp);

                return lambda;
        }
}
*/

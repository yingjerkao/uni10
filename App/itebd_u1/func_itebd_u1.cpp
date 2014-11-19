#include <iostream>
#include <cstdlib>
#include <mkl.h>
#include <uni10.hpp>

#include "func_init_u1.h"
#include "func_net_u1.h"
#include "func_evol_u1.h"

//============================

void itebd(uni10::UniTensor hamiltonian, int D, double dt_init, int ite_max, bool load_file){

	std::vector<uni10::Bond> bd = hamiltonian.bond();
	int d = bd[0].dim();	// physical dimension

	// deploy gamma, lambda tensors
	std::vector<uni10::UniTensor> gamma;
	std::vector<uni10::UniTensor> lambda;
	std::vector<uni10::UniTensor> expV;

/*
	if (load_file) {	// load gamma lambda from file
		gamma.push_back(loadGamma("gammaA", D));
                gamma.push_back(loadGamma("gammaB", D));
                lambda.push_back(loadLambda("lambdaA", D));
                lambda.push_back(loadLambda("lambdaB", D));
	}
*/	//else {	// randomize gamma lambda
		std::srand(5566);       // random seed

		for (int i=0; i<2; ++i){

			gamma.push_back( initGamma( bd[0].Qlist(), D, i==0 ) );
			gamma[i].randomize();

			lambda.push_back( initLambda( bd[0].Qlist(), D, i==0 ) );
			lambda[i].randomize();
        	}
	//}

	uni10::UniTensor eBH = expBO( -1.0 * dt_init, hamiltonian );       // generate exp(-Hdt)

	int Nite;
	for (Nite = 0; Nite < ite_max; ++Nite) {

		// ====== evolve and update =======
		for (int n=0; n<=1; n++){

			int p = (n+1)%2;
			uni10::UniTensor ket = netLGLGL(gamma[n], gamma[p], lambda[n], lambda[p]);
			ket = tensorOp(ket, eBH);
			truncUpdate(ket, gamma[n], gamma[p], lambda[n], lambda[p]);
		}
	}

	//====== calculate energy ======
	for (int n=0; n<=1; n++){

		int p = (n+1)%2;
		uni10::UniTensor ket = netLGLGL(gamma[n], gamma[p], lambda[n], lambda[p]);
		expV.push_back( expect(ket, hamiltonian) );
	}

	double avgV = 0.5 * (expV[0][0] + expV[1][0]);
	std::cout << "\ndt = " << dt_init << std::endl << "# of iteration = " << Nite << "\n";
	std::cout << "Energy = " << std::setprecision(10) << avgV << "\n\n";

	gamma[0].save("gammaA");
	gamma[1].save("gammaB");
	lambda[0].save("lambdaA");
	lambda[1].save("lambdaB");

	gamma.clear();
	lambda.clear();
	expV.clear();
}


//============================

/* iTEBD func ( hamiltonian, bond dimension, initial dt, lower bound for dt, max # of Trotter steps, load gamma/lambda from file ) */
void itebdGS(uni10::UniTensor hamiltonian, int D, double dt_init, double dt_min, int ite_max, bool load_file){

	double _dt = -1.0 * dt_init;      // initial -dt
        double totT = 0;        // total time evolved
        double last_expV = 0;   // expV of previous round

        int ckpt = 0;   // loop condition check factor
        int Nite = 0;   // # of iteration

	std::vector<uni10::Bond> bd = hamiltonian.bond();
	int d = bd[0].dim();	// physical dimension

	// deploy gamma, lambda tensors
	std::vector<uni10::UniTensor> gamma;
	std::vector<uni10::UniTensor> lambda;
	std::vector<uni10::UniTensor> expV;

/*
	if (load_file) {	// load gamma lambda from file
		gamma.push_back(loadGamma("gammaA", D));
                gamma.push_back(loadGamma("gammaB", D));
                lambda.push_back(loadLambda("lambdaA", D));
                lambda.push_back(loadLambda("lambdaB", D));
	}
*/	//else {	// randomize gamma lambda
		std::srand(5566);       // random seed

		for (int i=0; i<2; ++i){

			gamma.push_back( initGamma( bd[0].Qlist(), D, i==0 ) );
			gamma[i].randomize();

			lambda.push_back( initLambda( bd[0].Qlist(), D, i==0 ) );
			lambda[i].randomize();
        	}
	//}


	while ( fabs(_dt) >= dt_min ) {	// loop until dt < dt_min

		ckpt = 0;
		Nite = 0;

		uni10::UniTensor eBH = expBO( _dt, hamiltonian );	// generate exp(-Hdt)

		while (ckpt<1){	// loop until check point condition (ckpt=1) satisfied

			// ====== evolve and update =======

			for (int n=0; n<=1; n++){

				int p = (n+1)%2;
				uni10::UniTensor ket = netLGLGL(gamma[n], gamma[p], lambda[n], lambda[p]);
				ket = tensorOp(ket, eBH);
				truncUpdate(ket, gamma[n], gamma[p], lambda[n], lambda[p]);
			}


if ( ( Nite > 0 && Nite%10 == 0 ) || Nite == ite_max ){
			//====== calculate energy ======
			for (int n=0; n<=1; n++){

				int p = (n+1)%2;
				uni10::UniTensor ket = netLGLGL(gamma[n], gamma[p], lambda[n], lambda[p]);
				expV.push_back( expect(ket, hamiltonian) );
			}

		/*	if ( Nite%50==0 )
				std::cout << fabs(_dt) << "\t" << std::setprecision(10) << expV[0][0] << "\t" << expV[1][0] << "\n"; 
		*/
			// ====== output ======

			totT += fabs(_dt);
			double avgV = 0.5 * (expV[0][0] + expV[1][0]);

			if ( Nite == ite_max || 
				( Nite>0 && fabs(last_expV - avgV)<( 1e-6 * fabs(_dt) ) ) ) {
        	                std::cout << "\ndt = " << -1 * _dt << std::endl << "# of iteration = " << Nite << std::endl;
                	        std::cout << "Energy = " << std::setprecision(10) << avgV << "\n\n";
                        	ckpt = 1;
	                }

			last_expV = avgV;
			expV.clear();
}

			Nite += 1;
		}

		_dt = _dt / 2.0;
	}

	//gamma[0].permute(gamma[0].label(), 2);
	//gamma[1].permute(gamma[1].label(), 1);

	gamma[0].save("gammaA");
	gamma[1].save("gammaB");
	lambda[0].save("lambdaA");
	lambda[1].save("lambdaB");

	gamma.clear();
	lambda.clear();
}

//======================================================


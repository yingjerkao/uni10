#include <iostream>
#include <sstream>
#include <uni10.hpp>

using namespace std;
using namespace uni10;


int main(int argc, char* argv[])
{

	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
		std::cerr << "Allowed options:" << std::endl;
		std::cerr << "-H itf arg  generate Hamiltonian of transverse Ising model with transverse field = arg" << std::endl;
		std::cerr << "-H xx  generate Hamiltonian of xx model" << std::endl;
		std::cerr << "-sx  generate Sx (sigma x) operator" << std::endl;
		std::cerr << "-sz  generate Sz (sigma z) operator" << std::endl;
		std::cerr << "-id  generate identity operator" << std::endl;
		return 1;
	}

	std::vector<Qnum> qnums;
	std::vector<Qnum> qnum2;
	std::vector<Qnum> qnum_2;
        std::vector<Bond> bonds;
	std::vector<Bond> bondsP;
	std::vector<Bond> bondsM;

	// define quantum numbers
	uni10::Qnum q1(1);	
	uni10::Qnum q_1(-1);
	uni10::Qnum q2(2);
	uni10::Qnum q_2(-2);
	
	qnums.push_back( q1 );
	qnums.push_back( q_1 );

	qnum2.push_back( q2 );
	qnum_2.push_back( q_2 );
	
	// define bond parameters
	uni10::Bond bd_in( BD_IN, qnums );
	uni10::Bond bd_out( BD_OUT, qnums );

	uni10::Bond bd_sup2( BD_OUT, qnum2 );
	uni10::Bond bd_sup_2( BD_OUT, qnum_2 );	

	bonds.push_back( bd_in );
	bonds.push_back( bd_out );

	bondsP.push_back( bd_in );
	bondsP.push_back( bd_sup2 );
	bondsP.push_back( bd_out );

	bondsM.push_back( bd_in );
	bondsM.push_back( bd_sup_2 );
	bondsM.push_back( bd_out );

	// define elemental matrices
	double Sp[] = {0, 1,\
                       0, 0};

        double Sm[] = {0, 0,\
                       1, 0};

	double Sx[] = {0, 1,\
                       1, 0};

	double Sz[] = {1, 0,\
                       0,-1};

        double Id[] = {1, 0,\
                       0, 1};


	if (std::string(argv[1]) == "-H") {

		if (std::string(argv[2]) == "itf") { // Make sure we aren't at the end of argv!			

			if (argc < 4) {
				std::cerr << "Usage: " << argv[0] << " -H itf [strength of transverse field]" << std::endl;
				return 1;
			}

			double h;	// transverse field
			std::stringstream(argv[3]) >> h;

		        uni10::UniTensor tsz1( bonds );
		        uni10::UniTensor tsz2( bonds );

		        int label1[2] = {0, 1};
		        int label2[2] = {2, 3};

		        tsz1.setLabel( label1 );
		        tsz2.setLabel( label2 );

		        tsz1.setRawElem(Sz);
		        tsz2.setRawElem(Sz);

		        uni10::UniTensor hamiltonian = -1 * tsz1 * tsz2;

			uni10::UniTensor tsx( bonds );
			uni10::UniTensor tid( bonds );
		        tsx.setLabel( label1 );
		        tid.setLabel( label2 );
			tsx.setRawElem(Sx);
		        tid.setRawElem(Id);

		        hamiltonian += 0.5 * h * tsx * tid;

			tid.assign( bonds );
			tsx.assign( bonds );
			tid.setLabel( label1 );
		        tsx.setLabel( label2 );
			tsx.setRawElem(Sx);
		        tid.setRawElem(Id);

		        hamiltonian += 0.5 * h * tid * tsx;

		        uni10::UniTensor op = hamiltonian;

			op.save("hamiltonian");
		        //op.save("ham-itf");

			std::cout << op;
		}

		else if (std::string(argv[2]) == "xxz") {

			uni10::UniTensor tsp( bondsP );
			uni10::UniTensor tsm( bondsM );

			int label1[3] = {0, 100, 1};
			int label2[3] = {2, 101, 3};

			std::vector<int> combine;
			combine.push_back(100);
			combine.push_back(1);
			combine.push_back(101);

			tsp.setLabel( label1 );
			tsm.setLabel( label2 );

			tsp.setRawElem(Sp);
			tsm.setRawElem(Sm);

			uni10::UniTensor spsm = tsp * tsm;
			spsm.combineBond( combine );
			//std::cout << spsm << "ckpt 1\n";

			tsm.assign( bondsM );
			tsp.assign( bondsP );

			tsm.setLabel( label1 );
		        tsp.setLabel( label2 );
			tsm.setRawElem(Sm);
		        tsp.setRawElem(Sp);

			uni10::UniTensor smsp = tsm * tsp;
			smsp.combineBond( combine );
			//std::cout << smsp << "ckpt 2\n";
			
			uni10::UniTensor hamiltonian = spsm + smsp;
			hamiltonian *= 0.5;

			int lab_4b[4] = {0, 2, 1, 3};
			hamiltonian.setLabel(lab_4b);
			//std::cout << "ckpt 3\n";

			uni10::UniTensor tsz1( bonds );
			uni10::UniTensor tsz2( bonds );

			int lab1[2] = {0, 1}; 
			int lab2[2] = {2, 3}; 
			tsz1.setLabel( lab1 );
			tsz2.setLabel( lab2 );

			tsz1.setRawElem(Sz);
			tsz2.setRawElem(Sz);

			hamiltonian += tsz1 * tsz2 * 0.25;

			uni10::UniTensor op = hamiltonian;

			op.save("hamiltonian");
			//op.save("ham-xx");

			std::cout << op;
		}

		else { // Uh-oh, there was no argument to the -H option.
			std::cerr << "-H model Hamiltonian. Valid choices: itf, xx" << std::endl;
			return 1;
		}
	}

	else if (std::string(argv[1]) == "-id") {

        	uni10::UniTensor id1( bonds );

        	int label1[2] = {0, 1};
        	id1.setLabel( label1 );

        	id1.setRawElem(Id);
		uni10::UniTensor op = id1;
		op.save("identity");

		std::cout << op;	
	}

	else if (std::string(argv[1]) == "-sz") {

		uni10::UniTensor op( bonds, "Sz");
		op.setRawElem(Sz);
		op.save("sz");

		std::cout << op;
	}

	else if (std::string(argv[1]) == "-sx") {

        	uni10::UniTensor op( bonds, "Sx");
	        op.setRawElem(Sx);

        	op.save("sx");

		std::cout << op;
	}

	else {

		std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
		std::cerr << "Allowed options:" << std::endl;
		std::cerr << "-H itf arg  generate Hamiltonian of transverse Ising model with transverse field = arg" << std::endl;
		std::cerr << "-H xx  generate Hamiltonian of xx model" << std::endl;
		std::cerr << "-sx  generate Sx (sigma x) operator" << std::endl;
		std::cerr << "-sz  generate Sz (sigma z) operator" << std::endl;
		std::cerr << "-id  generate identity operator" << std::endl;
		return 1;
	}


	return 0;
}


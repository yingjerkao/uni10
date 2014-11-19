#include <iostream>
#include <sstream>

#include "itebd_u1.h"

//============================

int main(int argc, char* argv[]){

	if (argc < 2) {
                std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
		std::cerr << "Allowed options:" << std::endl;
		std::cerr << "-H arg  model Hamiltonian" << std::endl;
		std::cerr << "-w  load Wavefunction (gamma's and lambda's) from current folder" << std::endl;
		std::cerr << "-m arg  Max bond dimension" << std::endl;
		std::cerr << "-t arg  initial evolution Timestep" << std::endl;
		std::cerr << "-s arg  number of evolution Steps" << std::endl;
		std::cerr << "-gs  evolve to Ground State using default method" << std::endl;
                return 1;
        }

	int scenario = 0;

	int bd_dim = 8;		// bond dimension
	std::string ham;	// hamiltonian
	double dt_init = 0.1;	// initial dt
	int ite_max = -1;	// number of evolution
	bool load_file = false;	// load gamma lambda from file?

	for (int i = 1; i < argc; ++i) {

		if (std::string(argv[i]) == "-H") {
                        if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                                ham = std::string(argv[i+1]);
                        }
                        else { // Uh-oh, there was no argument to the -H option.
                                std::cerr << "-H option requires the filename of the Hamiltonian." << std::endl;
                                return 1;
                        }
                }
		else if (std::string(argv[i]) == "-w") {
                        load_file = true;
                }
		else if (std::string(argv[i]) == "-m") {
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				std::istringstream(argv[i+1]) >> bd_dim;
			}
			else { // Uh-oh, there was no argument to the -m option.
				std::cerr << "-m option requires a positive integer bond dimension." << std::endl;
				return 1;
			}  
		}
		else if (std::string(argv[i]) == "-t") {
                        if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                                std::stringstream(argv[i+1]) >> dt_init;
                        }
                        else { // Uh-oh, there was no argument to the -t option.
                                std::cerr << "-t option requires a positive initial dt." << std::endl;
                                return 1;
                        }
                } 
		else if (std::string(argv[i]) == "-s") {
                        if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                                std::istringstream(argv[i+1]) >> ite_max;
                        }
                        else { // Uh-oh, there was no argument to the -s option.
                                std::cerr << "-s option requires a positive integer number of evolution." << std::endl;
                                return 1;
                        }
                }
		else if (std::string(argv[i]) == "-gs") {
                        scenario = 1;
                }
	}


	uni10::UniTensor hamiltonian(ham);

	std::cout << "\nStarting iTEBD with: bond dimension = " << bd_dim << ", initial dt = " << dt_init << "\n";

	if (scenario)
		itebdGS(hamiltonian, bd_dim, dt_init, 1e-8, ite_max, load_file);
	else
		itebd(hamiltonian, bd_dim, dt_init, ite_max, load_file);


	return 0;
}


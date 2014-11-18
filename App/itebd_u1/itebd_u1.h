#ifndef ITEBD_H
#define ITEBD_H

#include <mkl.h>
#include <uni10.hpp>

#include "func_init_u1.h"
#include "func_net_u1.h"
#include "func_evol_u1.h"

void itebd(uni10::UniTensor hamiltonian, int D, double dt_init, int ite_max, bool load_file);
void itebdGS(uni10::UniTensor hamiltonian, int D, double dt_init, double dt_min, int ite_max, bool load_file);

#endif

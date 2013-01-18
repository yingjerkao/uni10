#define DIAG_NUM 2
/*Did definition*/
#define DIAG_UT 0
#define DIAG_AC 1
/*End Did definition*/
#ifdef GPU
# include "Tensor.cu"
#else
# include "Tensor.cpp"
#endif
#include "netBasic.cpp"
#include "myDias.c"

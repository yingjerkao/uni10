#define DIAG_NUM 0
/*Did definition*/
/*End Did definition*/
#ifdef GPU
# include "Tensor.cu"
#else
# include "Tensor.cpp"
#endif
#include "netBasic.cpp"

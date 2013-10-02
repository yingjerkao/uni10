LINK := -lm -lmkl -lmkl_lapack -lguide -lpthread
MKLLIB=/opt/intel/mkl/10.1.1.019/lib/em64t
INC=/opt/intel/mkl/10.1.1.019/include

CUBLASLIB=/usr/local/cuda/lib64
CUBLASINC=/usr/local/cuda/include
CULALIB=/usr/local/cula/lib64
CULAINC=/usr/local/cula/include
G_LINK := -lcula_core -lcula_lapack -lcublas -lcudart -arch=sm_35 -lpthread

MKLFILE = $(MKLROOT)/lib/intel64
TENLIB := ../TensorLib



rsp: reshape.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)
demo1: demo1.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)
gfus: fusion.cu
	nvcc $< -o $@ -g $(G_LINK) -L$(CULALIB) -L$(CUBLASLIB) -I$(CULAINC) -I$(CUBLASINC) -I$(TENLIB)
fus: fusion.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)
	./$@


gmm: dgemm.cu
	nvcc $< -o $@ -g $(G_LINK) -L$(CULALIB) -L$(CUBLASLIB) -I$(CULAINC) -I$(CUBLASINC) -I$(TENLIB)
gasd: Ascend.cu ../TensorLib/Tensor.cu
	nvcc $< -o $@ -g $(G_LINK) -L$(CULALIB) -L$(CUBLASLIB) -I$(CULAINC) -I$(CUBLASINC) -I$(TENLIB)
uni: Unitary.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)
unidm: Uni_diagram.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)

cmm: dgemm.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)
asd: Ascend.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)


direct: Direct.cpp
	g++ $< -o $@ $(LINK) -L$(MKLLIB) -I$(INC) -I$(TENLIB)
rand: random.cpp
	g++ $< -o $@ $(LINK) -L$(MKLLIB) -I$(INC) -I$(TENLIB)
jj: J1J2.cpp
	g++ $< -o $@ $(LINK) -L$(MKLLIB) -I$(INC) -I$(TENLIB)
sym: symmetry.cpp
	icpc $< -o $@ -mkl=parallel -I$(TENLIB)


gpu_unidm: Uni_diagram.cu
	nvcc $< -o $@ -g $(G_LINK) -L$(CULALIB) -L$(CUBLASLIB) -I$(CULAINC) -I$(CUBLASINC) -I$(TENLIB)
gpu_direct: Direct.cu
	nvcc $< -o $@ -g $(G_LINK) -L$(CULALIB) -L$(CUBLASLIB) -I$(CULAINC) -I$(CUBLASINC) -I$(TENLIB)

gpu_rand: random.cu
	nvcc $< -o $@ -g $(G_LINK) -L$(CULALIB) -L$(CUBLASLIB) -I$(CULAINC) -I$(CUBLASINC) -I$(TENLIB)
gpu_jj: J1J2.cu
	nvcc $< -o $@ -g $(G_LINK) -L$(CULALIB) -L$(CUBLASLIB) -I$(CULAINC) -I$(CUBLASINC) -I$(TENLIB)
cln:
	rm -f fus gfus asd gasd uni unidm direct gpu_direct jj gpu_jj rand gpu_rand gpu_unidm sym

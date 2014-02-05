#include "../../numeric/myLapack.h"
#include "../Matrix.h"

std::ostream& operator<< (std::ostream& os, const Matrix_t& m){
	os << std::endl << m.Rnum << " x " << m.Cnum << " = " << m.elemNum;
	if(m.diag)
		os << ", Diagonal";
	os <<std::endl << std::endl;
	for(int i = 0; i < m.Rnum; i++){
		for(int j = 0; j < m.Cnum; j++)
			if(m.diag){
				if(i == j)
					os << std::setw(7) << std::fixed << std::setprecision(3) << m.elem[i];
				else
					os << std::setw(7) << std::fixed << std::setprecision(3) << 0.0;
			}
			else
				os << std::setw(7) << std::fixed << std::setprecision(3) << m.elem[i * m.Cnum + j];
		os << std::endl << std::endl;
	}
	return os;
}

Matrix_t::Matrix_t(const Matrix_t& _m): Rnum(_m.Rnum), Cnum(_m.Cnum), elemNum(_m.elemNum), diag(_m.diag), elem(NULL){
	if(elemNum){
		elem = (double*)malloc(elemNum * sizeof(double));
		memcpy(elem, _m.elem, elemNum * sizeof(double));
	}
}

Matrix_t::Matrix_t(int _Rnum, int _Cnum, double* _elem, bool _diag): Rnum(_Rnum), Cnum(_Cnum), elemNum(_Rnum * _Cnum), diag(_diag), elem(NULL){
	if(_diag)
		elemNum = _Rnum < _Cnum ? _Rnum : _Cnum;
	if(elemNum){
		elem = (double*)malloc(elemNum * sizeof(double));
		memcpy(elem, _elem, elemNum * sizeof(double));
	}
}

Matrix_t::Matrix_t(int _Rnum, int _Cnum, bool _diag): Rnum(_Rnum), Cnum(_Cnum), elemNum(_Rnum * _Cnum), diag(_diag), elem(NULL){
	if(_diag)
		elemNum = _Rnum < _Cnum ? _Rnum : _Cnum;
	if(elemNum){
		elem = (double*)malloc(elemNum * sizeof(double));
		memset(elem, 0, elemNum * sizeof(double));
	}
}

Matrix_t& Matrix_t::operator=(const Matrix_t& _m){
	Rnum = _m.Rnum;
	Cnum = _m.Cnum;
	elemNum = _m.elemNum;
	diag = _m.diag;
	elem = (double*)realloc(elem, Rnum * Cnum * sizeof(double));
	memcpy(elem, _m.elem, Rnum * Cnum * sizeof(double));
	return *this;
}

Matrix_t::~Matrix_t(){
	if(elem != NULL)
		free(elem);
}

int Matrix_t::row()const{
	return Rnum;
}

int Matrix_t::col()const{
	return Cnum;
}
size_t Matrix_t::getElemNum()const{
	return elemNum;
}

Matrix_t operator* (const Matrix_t& Ma, const Matrix_t& Mb){
	assert(Ma.Cnum == Mb.Rnum);
	if((!Ma.diag) && (!Mb.diag)){
		Matrix_t Mc(Ma.Rnum, Mb.Cnum);	
		myDgemm(Ma.elem, Mb.elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.elem);
		return Mc;
	}
	else if(Ma.diag && (!Mb.diag)){
		Matrix_t Mc(Ma.Rnum, Mb.Cnum);	
		for(int i = 0; i < Ma.elemNum; i++)
			for(int j = 0; j < Mb.Cnum; j++)
				Mc.elem[i * Mb.Cnum + j] = Ma.elem[i] * Mb.elem[i * Mb.Cnum + j];
		return Mc;
	}
	else if((!Ma.diag) && Mb.diag){
		Matrix_t Mc(Ma.Rnum, Mb.Cnum);	
		for(int i = 0; i < Ma.Rnum; i++)
			for(int j = 0; j < Mb.elemNum; j++)
				Mc.elem[i * Mb.Cnum + j] = Ma.elem[i * Ma.Cnum + j] * Mb.elem[j];
		return Mc;
	}
	else{
		Matrix_t Mc(Ma.Rnum, Mb.Cnum, true);	
		for(int i = 0; i < Ma.Rnum; i++)
			Mc.elem[i] = Ma.elem[i] * Mb.elem[i];
		return Mc;
	}
}
bool operator== (const Matrix_t& m1, const Matrix_t& m2){
	double diff;
	if(m1.elemNum == m2.elemNum){	
		for(int i = 0; i < m1.elemNum; i++){
			diff = fabs(m1.elem[i] - m2.elem[i]);
			if(diff > 1E-6)
				return false;
		}
	}
	else
		return false;
	return true;
}


Matrix_t& Matrix_t::operator*= (const Matrix_t& Mb){
	return *this = *this * Mb;
}

std::vector<Matrix_t> Matrix_t::diagonalize(){
	assert(Rnum == Cnum);
	assert(!diag);
	std::vector<Matrix_t> outs;
	Matrix_t Eig(Rnum, Cnum, true);
	Matrix_t EigV(Rnum, Cnum);
	syDiag(elem, Rnum, Eig.elem, EigV.elem);
	outs.push_back(Eig);
	outs.push_back(EigV);	
	return outs;
}

std::vector<Matrix_t> Matrix_t::svd(){
	assert(!diag);
	std::vector<Matrix_t> outs;
	int min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
	Matrix_t U(Rnum, min);
	Matrix_t S(min, min, true);
	Matrix_t VT(min, Cnum);
	myDgesvd(elem, Rnum, Cnum, U.elem, S.elem, VT.elem);
	outs.push_back(U);
	outs.push_back(S);
	outs.push_back(VT);
	return outs;
}

void Matrix_t::orthoRand(){
	assert(!diag);
	orthoRandomize(elem, Rnum, Cnum);
}

void Matrix_t::bzero(){
	if(elemNum)
		memset(elem, 0, elemNum * sizeof(double));
}

Matrix_t operator*(const Matrix_t& Ma, double a){
	Matrix_t Mb(Ma);
	vecScal(a, Mb.elem, Mb.elemNum);
	return Mb;
}

void Matrix_t::operator*= (double a){
	vecScal(a, elem, elemNum);
}

Matrix_t operator+(const Matrix_t& Ma, const Matrix_t& Mb){
	Matrix_t Mc(Ma);
	vecAdd(Mb.elem, Mc.elem, Ma.elemNum);
	return Mc;
}

void Matrix_t::operator+= (const Matrix_t& Mb){
	vecAdd(Mb.elem, elem, elemNum);
}

void Matrix_t::transpose(){
	if(!diag){
		double* oldElem = (double*)malloc(elemNum * sizeof(double));
		memcpy(oldElem, elem, elemNum * sizeof(double));
		myTranspose(oldElem, Rnum, Cnum, elem, 0);
		free(oldElem);
	}
	int tmp = Rnum;
	Rnum = Cnum;
	Cnum = tmp;
}
double Matrix_t::trace(){
	assert(Rnum == Cnum);
	double sum = 0;
	if(diag)
		for(int i = 0; i < elemNum; i++)
			sum += elem[i];
	else
		for(int i = 0; i < Rnum; i++)
			sum += elem[i * Cnum + i];
	return sum;
}
void Matrix_t::save(const std::string fname){
	FILE *fp = fopen(fname.c_str(), "w");
	assert(fp != NULL);
	fwrite(elem, sizeof(double), elemNum, fp);
	fclose(fp);
}
void Matrix_t::load(const std::string fname){
	FILE *fp = fopen(fname.c_str(), "r");
	assert(fp != NULL);
	fread(elem, sizeof(double), elemNum, fp);
	fclose(fp);
}
double& Matrix_t::operator[](size_t idx){
	assert(idx < elemNum);
	return elem[idx];
}

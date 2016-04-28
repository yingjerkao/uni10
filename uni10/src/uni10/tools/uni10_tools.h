/****************************************************************************
#define UNI10_DTYPE_STRING "double"
*  @file uni10_tools.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Header file for useful string functions and wrapper functions
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef UNI10_TOOLS_H
#define UNI10_TOOLS_H
#include <cstdint>
#include <string>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <sstream>
#include <complex>
#include <uni10/data-structure/uni10_struct.h>
namespace uni10{

extern size_t MEM_USAGE;
extern size_t GPU_MEM_USAGE;

const size_t UNI10_GPU_GLOBAL_MEM = ((size_t)5) * 1<<30;
const int UNI10_THREADMAX = 1024;
const int UNI10_BLOCKMAX = 65535;
void* elemAlloc(size_t memsize, bool& ongpu);
void* elemAllocForce(size_t memsize, bool ongpu);
void* elemCopy(void* des, const void* src, size_t memsize, bool des_ongpu, bool src_ongpu);
void elemFree(void* ptr, size_t memsize, bool ongpu);
void elemBzero(void* ptr, size_t memsize, bool ongpu);
void elemRand(double* elem, size_t N, bool ongpu);
std::vector<_Swap> recSwap(std::vector<int>& ord, std::vector<int>& ordF);
std::vector<_Swap> recSwap(std::vector<int>& ord);	//Given the reshape order out to in.
void setDiag(double* elem, double* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu);
void getDiag(double* elem, double* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu);
void* mvGPU(void* elem, size_t memsize, bool& ongpu);
void* mvCPU(void* elem, size_t memsize, bool& ongpu);
void syncMem(void** elemA, void** elemB, size_t memsizeA, size_t memsizeB, bool& ongpuA, bool& ongpuB);
void shrinkWithoutFree(size_t memsize, bool ongpu);
void reshapeElem(double* oldElem, int bondNum, size_t elemNum, size_t* offset, double* newElem);
double getElemAt(size_t idx, double* elem, bool ongpu);
void setElemAt(size_t idx, double val, double* elem, bool ongpu);
void propogate_exception(const std::exception& e, const std::string& func_msg);
std::string exception_msg(const std::string& msg);
double elemMax(double *elem, size_t ElemNum, bool ongpu);
double elemAbsMax(double *elem, size_t ElemNum, bool ongpu);
/***** Complex version *****/
std::complex<double> getElemAt(size_t idx, std::complex<double>* elem, bool ongpu);
void setElemAt(size_t idx, std::complex<double> val, std::complex<double>* elem, bool ongpu);
void elemRand(std::complex<double>* elem, size_t N, bool ongpu);
void elemCast(std::complex<double>* des, double* src, size_t N, bool des_ongpu, bool src_ongpu);
void elemCast(double *des, std::complex<double> *src, size_t N, bool des_ongpu, bool src_ongpu);
void setDiag(std::complex<double>* elem, std::complex<double>* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu);
void getDiag(std::complex<double>* elem, std::complex<double>* diag_elem, size_t M, size_t N, size_t diag_N, bool ongpu, bool diag_ongpu);
void reshapeElem(std::complex<double>* oldElem, int bondNum, size_t elemNum, size_t* offset, std::complex<double>* newElem);

// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}
// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}
// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

static inline void throwTypeError(rflag tp){
	if(tp == RNULL){
		std::ostringstream err;
		err<<"Set Matrix type with wrong flag. Please use RTYPE instead of RNULL." << std::endl 
			<< "In the file Block.cpp, line(" << __LINE__ << ")";
		throw std::runtime_error(exception_msg(err.str()));
	}
}

static inline void throwTypeError(cflag tp){
	if(tp == CNULL){
		std::ostringstream err;
		err<<"Set Matrix type with wrong flag. Please use CTYPE instead of CNULL." << std::endl 
			<< "In the file Block.cpp, line(" << __LINE__ << ")";
		throw std::runtime_error(exception_msg(err.str()));
	}
}

};	/* namespace uni10 */
#endif /* UNI10_TOOLS_H */

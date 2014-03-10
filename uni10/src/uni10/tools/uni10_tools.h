#ifndef UNI10_TOOLS_H
#define UNI10_TOOLS_H
#include <stdint.h>
#include <limits.h>
#include <string>
#include <assert.h>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <uni10/data-structure/uni10_struct.h>
namespace uni10{
void* myMalloc(void* ptr, size_t memsize, int& status);
void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st);
void myFree(void* ptr, size_t memsize, int status);
void membzero(void* ptr, size_t memsize, int status);
void randomNums(double* elem, int N, int status);
std::vector<_Swap> recSwap(int* ord, int n, int* ordF);
std::vector<_Swap> recSwap(int* _ord, int n);	//Given the reshape order out to in. 

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
};	/* namespace uni10 */	

#endif /* UNI10_TOOLS_H */

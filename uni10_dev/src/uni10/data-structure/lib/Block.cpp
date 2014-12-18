/****************************************************************************
*  @file Block.cpp
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao
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
*  @brief Implementation file of Block class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/data-structure/Block.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#include <uni10/tensor-network/Matrix.h>
//using namespace uni10::datatype;
namespace uni10{
Block::Block(): Rnum(0), Cnum(0), m_elemNum(0), diag(false), ongpu(false), m_elem(NULL){}
Block::Block(size_t _Rnum, size_t _Cnum, bool _diag): Rnum(_Rnum), Cnum(_Cnum), m_elemNum(_Rnum * _Cnum), diag(_diag), ongpu(false), m_elem(NULL){}
Block::Block(const Block& _b): Rnum(_b.Rnum), Cnum(_b.Cnum), m_elemNum(_b.m_elemNum), diag(_b.diag), ongpu(_b.ongpu), m_elem(_b.m_elem){}

size_t Block::row()const{
	return Rnum;
}

size_t Block::col()const{
	return Cnum;
}
size_t Block::elemNum()const{
	return m_elemNum;
}
double* Block::getElem()const{
	return m_elem;
}
void Block::save(const std::string& fname)const{
  try{
    FILE *fp = fopen(fname.c_str(), "w");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in writing to file '"<<fname<<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    double* elem = m_elem;
    if(ongpu){
      elem = (double*)malloc(m_elemNum * sizeof(double));
      elemCopy(elem, m_elem, m_elemNum * sizeof(double), false, ongpu);
    }
    fwrite(elem, sizeof(double), m_elemNum, fp);
    fclose(fp);
    if(ongpu)
      free(elem);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Block::save(std::string&):");
  }
}

std::vector<Matrix> Block::eigh()const{
  std::vector<Matrix> outs;
  try{
    if(!(Rnum == Cnum)){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a non-square matrix.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(diag){
      std::ostringstream err;
      err<<"Cannot perform eigenvalue decomposition on a diagonal matrix. Need not to do so.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    Matrix Eig(Rnum, Cnum, true, ongpu);
    Matrix EigV(Rnum, Cnum, false, ongpu);
    syDiag(m_elem, Rnum, Eig.m_elem, EigV.m_elem, ongpu);
    outs.push_back(Eig);
    outs.push_back(EigV);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Matrix::eigh():");
  }
	return outs;
}


Block::~Block(){}
std::ostream& operator<< (std::ostream& os, const Block& b){
  try{
    os << b.Rnum << " x " << b.Cnum << " = " << b.m_elemNum;
    if(b.diag)
      os << ", Diagonal";
    if(b.ongpu)
      os<< ", onGPU";
    os <<std::endl << std::endl;
    double* elem;
    if(b.ongpu){
      elem = (double*)malloc(b.m_elemNum * sizeof(double));
      elemCopy(elem, b.m_elem, b.m_elemNum * sizeof(double), false, b.ongpu);
    }
    else
      elem = b.m_elem;
    for(size_t i = 0; i < b.Rnum; i++){
      for(size_t j = 0; j < b.Cnum; j++)
        if(b.diag){
          if(i == j)
            os << std::setw(7) << std::fixed << std::setprecision(3) << elem[i];
          else
            os << std::setw(7) << std::fixed << std::setprecision(3) << 0.0;
        }
        else
          os << std::setw(7) << std::fixed << std::setprecision(3) << elem[i * b.Cnum + j];
      os << std::endl << std::endl;
    }
    if(b.ongpu)
      free(elem);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator<<(std::ostream&, uni10::Matrix&):");
  }
  return os;
}
};	/* namespace uni10 */

/****************************************************************************
*  @file Block.cpp
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
*  @brief Implementation file of Block class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/tensor-network/CMatrix.h>

typedef double Real;
typedef std::complex<double> Complex;

namespace uni10{
  
  /********************* going move **************************/	    
 
  Block::Block(muType _tp, size_t _Rnum, size_t _Cnum, bool _diag): m_type(_tp), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL), cm_elem(NULL){}
  Real* Block::getRealElem()const{return m_elem;}
  Complex* Block::getComplexElem()const{return cm_elem;}

  /*********************  OPERATOR **************************/	    
  std::ostream& operator<< (std::ostream& os, const Block& b){
    try{
      os << b.Rnum << " x " << b.Cnum << " = " << b.elemNum();
      if(b.typeID() == 0){
        os <<std::endl << std::endl;
        os << "This matrix is EMPTY" <<std::endl << std::endl;
      } 
      if(b.typeID() == 1){
        os << ", REAL";
        if(b.diag)
          os << ", Diagonal";
        if(b.ongpu)
          os<< ", onGPU";
        os <<std::endl << std::endl;
        Real* elem;
        if(b.ongpu){
          elem = (Real*)malloc(b.elemNum() * sizeof(Real));
          elemCopy(elem, b.m_elem, b.elemNum() * sizeof(Real), false, b.ongpu);
        }
        else
          elem = b.m_elem;
        for(size_t i = 0; i < b.Rnum; i++){
          for(size_t j = 0; j < b.Cnum; j++)
            if(b.diag){
              //#ifdef UNI10_COMPLEX
              if(i == j)
                os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i];
              else
                os << std::setw(17) << std::fixed << std::setprecision(3) << 0.0;
            }
            else
              os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i * b.Cnum + j];
          os << std::endl << std::endl;
        }
        if(b.ongpu)
          free(elem);
      }
      if(b.typeID() == 2){
        os << ", COMPLEX";
        if(b.diag)
          os << ", Diagonal";
        if(b.ongpu)
          os<< ", onGPU";
        os <<std::endl << std::endl;
        Complex* elem;
        if(b.ongpu){
          elem = (Complex*)malloc(b.elemNum() * sizeof(Complex));
          elemCopy(elem, b.cm_elem, b.elemNum() * sizeof(Complex), false, b.ongpu);
        }
        else
          elem = b.cm_elem;
        for(size_t i = 0; i < b.Rnum; i++){
          for(size_t j = 0; j < b.Cnum; j++)
            if(b.diag){
              if(i == j)
                os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i];
              else
                os << std::setw(17) << std::fixed << std::setprecision(3) << 0.0;
            }
            else
              os << std::setw(17) << std::fixed << std::setprecision(3) << elem[i * b.Cnum + j];
          os << std::endl << std::endl;
        }
        if(b.ongpu)
          free(elem);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator<<(std::ostream&, uni10::Matrix&):");
    }
    return os;
  }

  bool operator== (const Block& m1, const Block& m2){
    try{
      if(m1.m_type != m2.m_type)
        return false;
      else{
        double diff;
        if(m1.elemNum() == m2.elemNum()){
          for(size_t i = 0; i < m1.elemNum(); i++){
            diff = (m1.m_type == RL) ? std::abs(m1.m_elem[i] - m2.m_elem[i]) : std::abs(m1.cm_elem[i] -m2.cm_elem[i]);
            if(diff > 1E-12)
              return false;
          }
        }
        else
          return false;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator==(uni10::Matrix&, uni10::Matrix&):");
    }
    return true;
  }

  Matrix operator* (const Block& Ma, const Block& Mb){
    try{
      if(!(Ma.Cnum == Mb.Rnum)){
        std::ostringstream err;
        err<<"The dimensions of the two matrices do not match for matrix multiplication.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      // R*R
      if(Ma.typeID() == 1 && Mb.typeID() == 1){
        if((!Ma.diag) && (!Mb.diag)){
          Matrix Mc(RTYPE, Ma.Rnum, Mb.Cnum);
          matrixMul(Ma.m_elem, Mb.m_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.m_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
          return Mc;
        }
        else if(Ma.diag && (!Mb.diag)){
          Matrix Mc(Mb);
          Mc.resize(Ma.Rnum, Mb.Cnum);
          diagRowMul(Mc.m_elem, Ma.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else if((!Ma.diag) && Mb.diag){
          Matrix Mc(Ma);
          Mc.resize(Ma.Rnum, Mb.Cnum);
          diagColMul(Mc.m_elem, Mb.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else{
          Matrix Mc(RTYPE, Ma.Rnum, Mb.Cnum, true);
          Mc.set_zero(RTYPE);
          size_t min = std::min(Ma.elemNum(), Mb.elemNum());
          elemCopy(Mc.m_elem, Ma.m_elem, min * sizeof(Real), Mc.ongpu, Ma.ongpu);
          vectorMul(Mc.m_elem, Mb.m_elem, min, Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }else if(Ma.typeID() == 2 && Mb.typeID() == 2){
        if((!Ma.diag) && (!Mb.diag)){
          Matrix Mc(CTYPE, Ma.Rnum, Mb.Cnum);
          matrixMul(Ma.cm_elem, Mb.cm_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.cm_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
          return Mc;
        }
        else if(Ma.diag && (!Mb.diag)){
          Matrix Mc(Mb);
          Mc.resize(Ma.Rnum, Mb.Cnum);
          diagRowMul(Mc.cm_elem, Ma.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else if((!Ma.diag) && Mb.diag){
          Matrix Mc(Ma);
          Mc.resize(Ma.Rnum, Mb.Cnum);
          diagColMul(Mc.cm_elem, Mb.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else{
          Matrix Mc(CTYPE, Ma.Rnum, Mb.Cnum, true);
          Mc.set_zero(CTYPE);
          size_t min = std::min(Ma.elemNum(), Mb.elemNum());
          elemCopy(Mc.cm_elem, Ma.cm_elem, min * sizeof(Complex), Mc.ongpu, Ma.ongpu);
          vectorMul(Mc.cm_elem, Mb.cm_elem, min, Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }else if(Ma.typeID() == 1 && Mb.m_type == 2){
        Matrix _Ma(Ma);
        RtoC(_Ma);
        if((!_Ma.diag) && (!Mb.diag)){
          Matrix Mc(CTYPE, _Ma.Rnum, Mb.Cnum);
          matrixMul(_Ma.cm_elem, Mb.cm_elem, _Ma.Rnum, Mb.Cnum, _Ma.Cnum, Mc.cm_elem, _Ma.ongpu, Mb.ongpu, Mc.ongpu);
          return Mc;
        }
        else if(_Ma.diag && (!Mb.diag)){
          Matrix Mc(Mb);
          Mc.resize(_Ma.Rnum, Mb.Cnum);
          diagRowMul(Mc.cm_elem, _Ma.cm_elem, Mc.Rnum, Mc.Cnum, _Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else if((!_Ma.diag) && Mb.diag){
          Matrix Mc(_Ma);
          Mc.resize(_Ma.Rnum, Mb.Cnum);
          diagColMul(Mc.cm_elem, Mb.cm_elem, Mc.Rnum, Mc.Cnum, _Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else{
          Matrix Mc(CTYPE, _Ma.Rnum, Mb.Cnum, true);
          Mc.set_zero(CTYPE);
          size_t min = std::min(_Ma.elemNum(), Mb.elemNum());
          elemCopy(Mc.cm_elem, _Ma.cm_elem, min * sizeof(Complex), Mc.ongpu, _Ma.ongpu);
          vectorMul(Mc.cm_elem, Mb.cm_elem, min, Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      if(Ma.typeID() == 2 && Mb.typeID() == 1){
        Matrix _Mb(Mb);
        RtoC(_Mb);
        if((!Ma.diag) && (!_Mb.diag)){
          Matrix Mc(CTYPE, Ma.Rnum, _Mb.Cnum);
          matrixMul(Ma.cm_elem, _Mb.cm_elem, Ma.Rnum, _Mb.Cnum, Ma.Cnum, Mc.cm_elem, Ma.ongpu, _Mb.ongpu, Mc.ongpu);
          return Mc;
        }
        else if(Ma.diag && (!_Mb.diag)){
          Matrix Mc(_Mb);
          Mc.resize(Ma.Rnum, _Mb.Cnum);
          diagRowMul(Mc.cm_elem, Ma.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else if((!Ma.diag) && _Mb.diag){
          Matrix Mc(Ma);
          Mc.resize(Ma.Rnum, _Mb.Cnum);
          diagColMul(Mc.cm_elem, _Mb.cm_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
          return Mc;
        }
        else{
          Matrix Mc(CTYPE, Ma.Rnum, _Mb.Cnum, true);
          Mc.set_zero(CTYPE);
          size_t min = std::min(Ma.elemNum(), _Mb.elemNum());
          elemCopy(Mc.cm_elem, Ma.cm_elem, min * sizeof(Complex), Mc.ongpu, Ma.ongpu);
          vectorMul(Mc.cm_elem, _Mb.cm_elem, min, Mc.ongpu, _Mb.ongpu);
          return Mc;
        }
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::Matrix&, uni10::Matrix&):");
      return Block();
    }
  }
  
  Matrix operator*(const Block& Ma, double a){
    try{
      Matrix Mb(Ma);
      if(Ma.typeID() != 0)
        (Ma.typeID() == 1) ? vectorScal(a, Mb.m_elem, Mb.elemNum(), Mb.isOngpu()): vectorScal(a, Mb.cm_elem, Mb.elemNum(), Mb.isOngpu());
      return Mb;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::Matrix&, double):");
      return Block();
    }
  }
  Matrix operator*(double a, const Block& Ma){return Ma * a;}

  Matrix operator*(const Block& Ma, const std::complex<double>& a){
    try{
      Matrix Mb(Ma);
      if(a.imag() == 0){
        double _a = a.real();
        return _a * Mb; 
      }else{
        RtoC(Mb);
        vectorScal(a, Mb.cm_elem, Mb.elemNum(), Mb.isOngpu());
      }
      return Mb;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator*(uni10::Matrix&, Complex):");
      return Block();
    }
  }
  Matrix operator*(const std::complex<double>& a, const Block& Ma){return Ma * a;}
  
  Matrix operator+(const Block& Ma, const Block& Mb){
    /*if ((Ma.diag && !Mb.diag) || (!Ma.diag && Mb.diag) ){
      std::ostringstream err;
      err<<"Fatal error(code = T1). Addition of a diagonal and a regular matrices is not implemented.";
      throw std::runtime_error(exception_msg(err.str())); ;
      };*/

    try{
      if(Ma.typeID() == 1 && Mb.typeID() == 1){
        if (Ma.diag && !Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(RTYPE, Ma.Rnum,Ma.Cnum);
          setDiag(Mc.m_elem,Ma.m_elem,Mc.Rnum,Mc.Cnum,Ma.Rnum,Mc.ongpu,Ma.ongpu);
          vectorAdd(Mc.m_elem, Mb.m_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        } else if (Mb.diag && !Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;
          Matrix Mc(RTYPE, Mb.Rnum, Mb.Cnum);
          setDiag(Mc.m_elem,Mb.m_elem,Mc.Rnum,Mc.Cnum,Mb.Cnum,Mc.ongpu,Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.m_elem, Ma.m_elem, Mc.elemNum(), Mc.ongpu, Ma.ongpu);
          return Mc;
        } else {
          Matrix Mc(Ma);
          vectorAdd(Mc.m_elem, Mb.m_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }else if(Ma.typeID() == 2 && Mb.typeID() == 2){
        if (Ma.diag && !Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(CTYPE, Ma.Rnum,Ma.Cnum);
          setDiag(Mc.cm_elem,Ma.cm_elem,Mc.Rnum,Mc.Cnum,Ma.Rnum,Mc.ongpu,Ma.ongpu);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        } else if (Mb.diag && !Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;
          Matrix Mc(CTYPE, Mb.Rnum,Mb.Cnum);
          setDiag(Mc.cm_elem,Mb.cm_elem,Mc.Rnum,Mc.Cnum,Mb.Cnum,Mc.ongpu,Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.cm_elem, Ma.cm_elem, Mc.elemNum(), Mc.ongpu, Ma.ongpu);
          return Mc;
        } else {
          Matrix Mc(Ma);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }else if(Ma.typeID() == 1 && Mb.typeID() == 22){
        Matrix _Ma(Ma);
        RtoC(_Ma);
        if (_Ma.diag && !Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(CTYPE, _Ma.Rnum,_Ma.Cnum);
          setDiag(Mc.cm_elem,_Ma.cm_elem,Mc.Rnum,Mc.Cnum,_Ma.Rnum,Mc.ongpu,_Ma.ongpu);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        } else if (Mb.diag && !_Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;
          Matrix Mc(CTYPE, Mb.Rnum,Mb.Cnum);
          setDiag(Mc.cm_elem,Mb.cm_elem,Mc.Rnum,Mc.Cnum,Mb.Cnum,Mc.ongpu,Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.cm_elem, _Ma.cm_elem, Mc.elemNum(), Mc.ongpu, _Ma.ongpu);
          return Mc;
        } else {
          Matrix Mc(_Ma);
          vectorAdd(Mc.cm_elem, Mb.cm_elem, Mc.elemNum(), Mc.ongpu, Mb.ongpu);
          return Mc;
        }
      }
      if(Ma.typeID() == 2 && Mb.typeID() == 1){
        Matrix _Mb(Mb);
        RtoC(_Mb);
        if (Ma.diag && !_Mb.diag) {
          //std::cout << "Ma is diagonal."<<std::endl;
          Matrix Mc(CTYPE, Ma.Rnum,Ma.Cnum);
          setDiag(Mc.cm_elem,Ma.cm_elem,Mc.Rnum,Mc.Cnum,Ma.Rnum,Mc.ongpu,Ma.ongpu);
          vectorAdd(Mc.cm_elem, _Mb.cm_elem, Mc.elemNum(), Mc.ongpu, _Mb.ongpu);
          return Mc;
        } else if (_Mb.diag && !Ma.diag) {
          //std::cout << "Mb is diagonal."<<std::endl;
          Matrix Mc(CTYPE, _Mb.Rnum,_Mb.Cnum);
          setDiag(Mc.cm_elem,_Mb.cm_elem,Mc.Rnum,Mc.Cnum,_Mb.Cnum,Mc.ongpu,_Mb.ongpu);
          //std::cout << Mc ;
          vectorAdd(Mc.cm_elem, Ma.cm_elem, Mc.elemNum(), Mc.ongpu, Ma.ongpu);
          return Mc;
        } else {
          Matrix Mc(Ma);
          vectorAdd(Mc.cm_elem, _Mb.cm_elem, Mc.elemNum(), Mc.ongpu, _Mb.ongpu);
          return Mc;
        }
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator+(uni10::Matrix&, uni10::Matrix&):");
      return Block();
    }
  }
  
  Complex Block::operator[](size_t idx)const{
    try{
      if(!(idx < elemNum())){
        std::ostringstream err;
        err<<"Index exceeds the number of elements("<<elemNum()<<").";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 0){
        std::ostringstream err;
        err<<"This matrix is EMPTY";
        throw std::runtime_error(exception_msg(err.str()));
      }
      else if(typeID() == 1)
        return Complex(getElemAt(idx, m_elem, ongpu), 0);
      else if(typeID() == 2) 
        return getElemAt(idx, cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::operator[](size_t):");
      return 0;
    }
  }
 
  /*********************  NO TYPE **************************/	    

  Block::Block():r_flag(RNULL), c_flag(CNULL), Rnum(0), Cnum(0), diag(false), ongpu(false), m_elem(NULL), cm_elem(NULL){}

  Block::Block(size_t _Rnum, size_t _Cnum, bool _diag): r_flag(RNULL), c_flag(CNULL), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL), cm_elem(NULL){}
  
  Block::Block(int _typeID, size_t _Rnum, size_t _Cnum, bool _diag): r_flag(RNULL), c_flag(CNULL), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL), cm_elem(NULL){
    if(_typeID == 1)
      r_flag = RTYPE; 
    if(_typeID == 2)
      c_flag = CTYPE; 
  }
  
  Block::Block(const Block& _b): r_flag(_b.r_flag), c_flag(_b.c_flag), Rnum(_b.Rnum), Cnum(_b.Cnum), diag(_b.diag), ongpu(_b.ongpu), m_elem(_b.m_elem), cm_elem(_b.cm_elem){}

  Block::~Block(){}

  size_t Block::row()const{return Rnum;}

  size_t Block::col()const{return Cnum;}

  bool Block::isDiag()const{return diag;}

  bool Block::isOngpu()const{return ongpu;}

  size_t Block::elemNum()const{
    if(diag)
      return (Rnum < Cnum ? Rnum : Cnum);
    else
      return Rnum * Cnum;
  }

  int Block::typeID()const{
    return r_flag + c_flag;
  }
  
  void Block::savePrototype(const std::string& fname)const{
    FILE *fp = fopen(fname.c_str(), "w");
    if(!(fp != NULL)){
      std::ostringstream err;
      err<<"Error in writing to file '"<<fname<<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(typeID() == 0){ 
      fwrite(&r_flag, sizeof(r_flag), 1, fp);
      fwrite(&c_flag, sizeof(c_flag), 1, fp);
      fwrite(&Rnum, sizeof(Rnum), 1, fp);
      fwrite(&Cnum, sizeof(Cnum), 1, fp);
      fwrite(&diag, sizeof(diag), 1, fp);
      fwrite(&ongpu, sizeof(ongpu), 1, fp);
    }
    fclose(fp);
  } 
  
  void Block::save(const std::string& fname)const{
    try{
      if(typeID() == 0) 
        savePrototype(fname);
      if(typeID() == 1) 
        save(RTYPE, fname);
      if(typeID() == 2) 
        save(CTYPE, fname);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::save(std::string&):");
    }
  }

  std::vector<Matrix> Block::qr()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QR decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform QR decomposition on EMPTY matrix . Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return qr(RTYPE);
      if(typeID() == 2)
        return qr(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::qr():");
    }
    return outs;
  }
  

  std::vector<Matrix> Block::rq()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform RQ decomposition on EMPTY matrix . Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return rq(RTYPE);
      if(typeID() == 2)
        return rq(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::rq():");
    }
    return outs;
  }

  std::vector<Matrix> Block::ql()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QL decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform QL decomposition on EMPTY matrix . Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return ql(RTYPE);
      if(typeID() == 2)
        return ql(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::ql():");
    }
    return outs;
  }
  
  

  std::vector<Matrix> Block::lq()const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      } 
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform LQ decomposition on EMPTY matrix . Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return lq(RTYPE);
      if(typeID() == 2)
        return lq(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::lq():");
    }
    return outs;
  }

  std::vector<Matrix> Block::svd()const{
    std::vector<Matrix> outs;
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform SVD decomposition on EMPTY matrix . Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return svd(RTYPE);
      if(typeID() == 2)
        return svd(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::svd():");
    }
    return outs;
  }
  
  
  Real Block::norm()const{
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"This matirx is empty" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 1)
        return vectorNorm(m_elem, elemNum(), 1, ongpu);
      if(typeID() == 2)
        return vectorNorm(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::norm():");
      return 0;
    }
  }
  
  Matrix Block::inverse()const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform inversion on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 0){
        std::ostringstream err;
        err<<"This matirx is empty" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
        throw std::runtime_error(exception_msg(err.str()));
      }
      Matrix invM(*this);
      assert(ongpu == invM.isOngpu());
      if(typeID() == 1)
        matrixInv(invM.m_elem, Rnum, invM.diag, invM.ongpu);
      if(typeID() == 2)
        matrixInv(invM.cm_elem, Rnum, invM.diag, invM.ongpu);
      return invM;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::inverse():");
      return Matrix();
    }
  }
  
  Matrix Block::getDiag()const{
    try{
      if(diag)
        return *this;
      else{
        int _typeID = typeID();
        if(_typeID == 0){
          std::ostringstream err;
          err<<"This matirx is empty" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
          throw std::runtime_error(exception_msg(err.str()));
        }
        if(_typeID == 1)
          getDiag(RTYPE);
        if(_typeID == 2)
          getDiag(CTYPE);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::getDiag():");
      return Matrix();
    }
  }
  
  Complex Block::trace()const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform trace on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform trace on an EMPTY matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      else if(typeID() == 1)
        return Complex(trace(RTYPE), 0);
      else if(typeID() == 2)
        return trace(CTYPE);
    }catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::trace():");
      return 0;
    }
  }

  Complex Block::sum()const{
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform trace on an EMPTY matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      else if(typeID() == 1)
        return Complex(vectorSum(m_elem, elemNum(), 1, ongpu), 0);
      else if(typeID() == 2)
        return vectorSum(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::sum():");
      return 0;
    }
  }
  
  std::vector<Matrix> Block::eig()const{
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform eigenvalue decomposition on an EMPTY matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      else if(typeID() == 1)
        return eig(RTYPE);
      else if(typeID() == 2)
        return eig(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eig():");
    }
  }
  
  std::vector<Matrix> Block::eigh()const{
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform eigenvalue decomposition on an EMPTY matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      else if(typeID() == 1)
        return eigh(RTYPE);
      else if(typeID() == 2)
        return eigh(CTYPE);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eigh():");
    }
  }
 

  Complex Block::at(size_t r, size_t c)const{
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"This matrix is EMPTY.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      else if(typeID() == 1)
        return Complex(at(RTYPE, r, c), 0);
      else if(typeID() == 2)
        return at(CTYPE, r, c);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::at(size_t, size_t):");
      return 0;
    }
  }
  /*********************NEE*****************************/
  /*********************  REAL **********************/

  Block::Block(rflag _tp, size_t _Rnum, size_t _Cnum, bool _diag): r_flag(_tp), c_flag(CNULL), m_type(EMPTY), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL), cm_elem(NULL){}
	    
  Real* Block::getElem(rflag _tp)const{return m_elem;}
  
  void Block::save(rflag _tp, const std::string& fname)const{
    try{
      FILE *fp = fopen(fname.c_str(), "w");
      if(!(fp != NULL)){
        std::ostringstream err;
        err<<"Error in writing to file '"<<fname<<"'.";
        throw std::runtime_error(exception_msg(err.str()));
      }
  //    fwrite(&m_type, sizeof(m_type), 1, fp);
      fwrite(&r_flag, sizeof(r_flag), 1, fp);
      fwrite(&c_flag, sizeof(c_flag), 1, fp);
      fwrite(&Rnum, sizeof(Rnum), 1, fp);
      fwrite(&Cnum, sizeof(Cnum), 1, fp);
      fwrite(&diag, sizeof(diag), 1, fp);
      fwrite(&ongpu, sizeof(ongpu), 1, fp);
      Real* elem = m_elem;
      if(ongpu){
        elem = (Real*)malloc(elemNum() * sizeof(Real));
        elemCopy(elem, m_elem, elemNum() * sizeof(Real), false, ongpu);
      }
      fwrite(elem, sizeof(Real), elemNum(), fp);
      if(ongpu)
        free(elem);
      fclose(fp);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::save(std::string&):");
    }
  }
  
  std::vector<Matrix> Block::qr(rflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QR decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(RTYPE, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(RTYPE, Cnum, Cnum, false, ongpu));
      if(!diag)
          matrixQR(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem);
      else{
        size_t min = std::min(Rnum, Cnum);
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        for(int i = 0; i < min; i++)
          tmpR[i*min+i] = m_elem[i];
        matrixQR(tmpR, min, min, outs[0].m_elem, outs[1].m_elem);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::qr():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::rq(rflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(RTYPE, Rnum, Rnum, false, ongpu)); //r
      outs.push_back(Matrix(RTYPE, Rnum, Cnum, false, ongpu)); //q
      if(!diag){
          matrixRQ(m_elem, Rnum, Cnum, outs[1].m_elem, outs[0].m_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        for(int i = 0; i < min; i++)
          tmpR[i*min+i] = m_elem[i];
        matrixRQ(tmpR, min, min, outs[1].m_elem, outs[0].m_elem);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::rq():");
    }
    return outs;
  }

  std::vector<Matrix> Block::lq(rflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      } 
      outs.push_back(Matrix(RTYPE, Rnum, Rnum, false, ongpu));
      outs.push_back(Matrix(RTYPE, Rnum, Cnum, false, ongpu));
      if(!diag){
        matrixLQ(m_elem, Rnum, Cnum, outs[1].m_elem, outs[0].m_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        for(int i = 0; i < min; i++)
          tmpR[i*min+i] = m_elem[i];
        matrixLQ(tmpR, min, min, outs[1].m_elem, outs[0].m_elem);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::lq():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::ql(rflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QL decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(RTYPE, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(RTYPE, Cnum, Cnum, false, ongpu));
      if(!diag){
        matrixQL(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        for(int i = 0; i < min; i++)
          tmpR[i*min+i] = m_elem[i];
        matrixQL(tmpR, min, min, outs[0].m_elem, outs[1].m_elem);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::ql():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::svd(rflag _tp)const{
    std::vector<Matrix> outs;
    try{
    /*  
      if(diag){
        std::ostringstream err;
        err<<"Cannot perform singular value decomposition on a diagonal matrix. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    */
      size_t min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
      //GPU_NOT_READY
      outs.push_back(Matrix(RTYPE, Rnum, min, false, ongpu));
      outs.push_back(Matrix(RTYPE, min, min, true, ongpu));
      outs.push_back(Matrix(RTYPE, min, Cnum, false, ongpu));
      if(!diag){
          matrixSVD(m_elem, Rnum, Cnum, outs[0].m_elem, outs[1].m_elem, outs[2].m_elem, ongpu);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Real* tmpR = (Real*)calloc(min*min, sizeof(Real));
        for(int i = 0; i < min; i++)
          tmpR[i*min+i] = m_elem[i];
        matrixSVD(tmpR, min, min, outs[0].m_elem, outs[1].m_elem, outs[2].m_elem, ongpu);
        free(tmpR);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::svd():");
    }
    return outs;
  }
  
  Real Block::norm(rflag _tp)const{
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"This matirx is empty" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
        throw std::runtime_error(exception_msg(err.str()));
      }
      return vectorNorm(m_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::norm():");
      return 0;
    }
  }
  
  Matrix Block::inverse(rflag _tp)const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform inversion on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      Matrix invM(*this);
      assert(ongpu == invM.isOngpu());
      matrixInv(invM.m_elem, Rnum, invM.diag, invM.ongpu);
      return invM;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::inverse():");
      return Matrix();
    }
  }
  
  Matrix Block::getDiag(rflag _tp)const{
    try{
      if(diag)
        return *this;
      else{
        Matrix D(RTYPE, Rnum, Cnum, true, ongpu);
        ::uni10::getDiag(m_elem, D.getElem(RTYPE), Rnum, Cnum, D.elemNum(), ongpu, D.isOngpu());
        return D;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::getDiag():");
      return Matrix();
    }
  }
  
  Real Block::trace(rflag _tp)const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform trace on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(diag)
        return vectorSum(m_elem, elemNum(), 1, ongpu);
      else
        return vectorSum(m_elem, Cnum, Cnum + 1, ongpu);
    }catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::trace():");
      return 0;
    }
  }
  
  Real Block::sum(rflag _tp)const{
    try{
      return vectorSum(m_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::sum():");
      return 0;
    }
  }
  
  std::vector<Matrix> Block::eig(rflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform eigenvalue decomposition on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(diag){
        std::ostringstream err;
        err<<"Cannot perform eigenvalue decomposition on a diagonal matrix. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      //GPU_NOT_READY
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      eigDecompose(m_elem, Rnum, outs[0].cm_elem, outs[1].cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eig():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::eigh(rflag _tp)const{
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
      //GPU_NOT_READY
      outs.push_back(Matrix(RTYPE, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(RTYPE, Rnum, Cnum, false, ongpu));
      Matrix Eig(RTYPE, Rnum, Cnum, true, ongpu);
      eigSyDecompose(m_elem, Rnum, Eig.m_elem, outs[1].m_elem, ongpu);
      outs[0] = Eig;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eigh():");
    }
    return outs;
  }

  Real Block::at(rflag _tp, size_t r, size_t c)const{
    try{
      if(!((r < Rnum) && (c < Cnum))){
        std::ostringstream err;
        err<<"The input indices are out of range.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(diag){
        if(!(r == c && r < elemNum())){
          std::ostringstream err;
          err<<"The matrix is diagonal, there is no off-diagonal element.";
          throw std::runtime_error(exception_msg(err.str()));
        }
        return getElemAt(r, m_elem, ongpu);
      }
      else
        return getElemAt(r * Cnum + c, m_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::at(size_t, size_t):");
      return 0;
    }
  }
  /*********************REE*****************************/
  /*********************  COMPLEX **********************/

  Block::Block(cflag _tp, size_t _Rnum, size_t _Cnum, bool _diag): r_flag(RNULL), c_flag(_tp), m_type(EMPTY), Rnum(_Rnum), Cnum(_Cnum), diag(_diag), ongpu(false), m_elem(NULL), cm_elem(NULL){}
  
  Complex* Block::getElem(cflag _tp)const{return cm_elem;}
  
  void Block::save(cflag _tp, const std::string& fname)const{
    try{
      FILE *fp = fopen(fname.c_str(), "w");
      if(!(fp != NULL)){
        std::ostringstream err;
        err<<"Error in writing to file '"<<fname<<"'.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      fwrite(&r_flag, sizeof(r_flag), 1, fp);
      fwrite(&c_flag, sizeof(c_flag), 1, fp);
      fwrite(&Rnum, sizeof(Rnum), 1, fp);
      fwrite(&Cnum, sizeof(Cnum), 1, fp);
      fwrite(&diag, sizeof(diag), 1, fp);
      fwrite(&ongpu, sizeof(ongpu), 1, fp);
      Complex* elem = cm_elem;
      if(ongpu){
        elem = (Complex*)malloc(elemNum() * sizeof(Complex));
        elemCopy(elem, cm_elem, elemNum() * sizeof(Complex), false, ongpu);
      }
      fwrite(elem, sizeof(Complex), elemNum(), fp);
      if(ongpu)
        free(elem);
      fclose(fp);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::save(std::string&):");
    }
  }
  
  std::vector<Matrix> Block::qr(cflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QR decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(CTYPE, Cnum, Cnum, false, ongpu));
      if(!diag)
          matrixQR(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem);
      else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(int i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixQR(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::qr():");
    }
    return outs;
  }
 
  std::vector<Matrix> Block::rq(cflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(CTYPE, Rnum, Rnum, false, ongpu)); //r
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu)); //q
      if(!diag){
        matrixRQ(cm_elem, Rnum, Cnum, outs[1].cm_elem, outs[0].cm_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(int i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixRQ(tmpC, min, min, outs[1].cm_elem, outs[0].cm_elem);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::rq():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::lq(cflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum > Cnum){
        std::ostringstream err;
        err<<"Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      } 
      outs.push_back(Matrix(CTYPE, Rnum, Rnum, false, ongpu));
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      if(!diag){
        matrixLQ(cm_elem, Rnum, Cnum, outs[1].cm_elem, outs[0].cm_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(int i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixLQ(tmpC, min, min, outs[1].cm_elem, outs[0].cm_elem);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::lq():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::ql(cflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(Rnum < Cnum){
        std::ostringstream err;
        err<<"Cannot perform QL decomposition when Rnum < Cnum. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      outs.push_back(Matrix(CTYPE, Cnum, Cnum, false, ongpu));
      if(!diag){
        matrixQL(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(int i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixQL(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::ql():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::svd(cflag _tp)const{
    std::vector<Matrix> outs;
    try{
    /*  
      if(diag){
        std::ostringstream err;
        err<<"Cannot perform singular value decomposition on a diagonal matrix. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
    */
      size_t min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
      //GPU_NOT_READY
      outs.push_back(Matrix(CTYPE, Rnum, min, false, ongpu));
      outs.push_back(Matrix(CTYPE, min, min, true, ongpu));
      outs.push_back(Matrix(CTYPE, min, Cnum, false, ongpu));
      if(!diag){
        matrixSVD(cm_elem, Rnum, Cnum, outs[0].cm_elem, outs[1].cm_elem, outs[2].cm_elem, ongpu);
      }else{
        size_t min = std::min(Rnum, Cnum);
        Complex* tmpC = (Complex*)calloc(min*min , sizeof(Complex));
        for(int i = 0; i < min; i++)
          tmpC[i*min+i] = cm_elem[i];
        matrixSVD(tmpC, min, min, outs[0].cm_elem, outs[1].cm_elem, outs[2].cm_elem, ongpu);
        free(tmpC);
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::svd():");
    }
    return outs;
  }
  
  Real Block::norm(cflag _tp)const{
    try{
      if(typeID() == 0){
        std::ostringstream err;
        err<<"This matirx is empty" << std::endl << "In the file Block.cpp, line(" << __LINE__ << ")";
        throw std::runtime_error(exception_msg(err.str()));
      }
      return vectorNorm(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::norm():");
      return 0;
    }
  }
  
  Matrix Block::inverse(cflag _tp)const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform inversion on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      Matrix invM(*this);
      assert(ongpu == invM.isOngpu());
      matrixInv(invM.cm_elem, Rnum, invM.diag, invM.ongpu);
      return invM;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::inverse():");
      return Matrix();
    }
  }
  
  Matrix Block::getDiag(cflag _tp)const{
    try{
      if(diag)
        return *this;
      else{
        Matrix D(CTYPE, Rnum, Cnum, true, ongpu);
        ::uni10::getDiag(cm_elem, D.getElem(CTYPE), Rnum, Cnum, D.elemNum(), ongpu, D.isOngpu());
        return D;
      }
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::getDiag():");
      return Matrix();
    }
  }
  
  Complex Block::trace(cflag _tp)const{
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform trace on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(diag)
        return vectorSum(cm_elem, elemNum(), 1, ongpu);
      else
        return vectorSum(cm_elem, Cnum, Cnum + 1, ongpu);
    }catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::trace():");
      return 0;
    }
  }
  
  Complex Block::sum(cflag _tp)const{
    try{
      return vectorSum(cm_elem, elemNum(), 1, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::sum():");
      return 0;
    }
  }
  
  std::vector<Matrix> Block::eig(cflag _tp)const{
    std::vector<Matrix> outs;
    try{
      if(!(Rnum == Cnum)){
        std::ostringstream err;
        err<<"Cannot perform eigenvalue decomposition on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(diag){
        std::ostringstream err;
        err<<"Cannot perform eigenvalue decomposition on a diagonal matrix. Nothing to do.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      //GPU_NOT_READY
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      eigDecompose(cm_elem, Rnum, outs[0].cm_elem, outs[1].cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eig():");
    }
    return outs;
  }
  
  std::vector<Matrix> Block::eigh(cflag _tp)const{
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
      //GPU_NOT_READY
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, true, ongpu));
      outs.push_back(Matrix(CTYPE, Rnum, Cnum, false, ongpu));
      Matrix Eig(RTYPE, Rnum, Cnum, true, ongpu);
      eigSyDecompose(cm_elem, Rnum, Eig.m_elem, outs[1].cm_elem, ongpu);
      outs[0] = Eig;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::eigh():");
    }
    return outs;
  }
  
  Complex Block::at(cflag _tp, size_t r, size_t c)const{
    try{
      if(!((r < Rnum) && (c < Cnum))){
        std::ostringstream err;
        err<<"The input indices are out of range.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(diag){
        if(!(r == c && r < elemNum())){
          std::ostringstream err;
          err<<"The matrix is diagonal, there is no off-diagonal element.";
          throw std::runtime_error(exception_msg(err.str()));
        }
        return getElemAt(r, cm_elem, ongpu);
      }
      else
        return getElemAt(r * Cnum + c, cm_elem, ongpu);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Block::at(size_t, size_t):");
      return 0;
    }
  }
  /*********************CEE*****************************/
  /*****************************************************/
  
  std::ostream& operator<< (std::ostream& os, const muType& tp){
    try{
      if(tp == RL)
        os << "The matrix type is REAL.";
      if(tp == CX)
        os << "The matrix type is COMPLEX.";
      if(tp == EMPTY)
        os << "This matrix is EMPTY.";
      os << std::endl << std::endl;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function operator<<(std::ostream&, uni10::muType&):");
    }
    return os;
  }


  Real* Block::getElem()const{return m_elem;}
  


  muType Block::getType()const{
    return m_type;
  }

  
  


  
  
  
  
  
  
  

  
  size_t lanczosEigh(rflag _tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol){
    try{
      if(!(ori_mat.Rnum == ori_mat.Cnum)){
        std::ostringstream err;
        err<<"Cannot perform Lanczos algorithm to find the lowest eigen value and eigen vector on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(!(ori_mat.Rnum == psi.elemNum())){
        std::ostringstream err;
        err<<"Error in Lanczos initial vector psi. The vector dimension does not match with the number of the columns.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(ori_mat.ongpu && !psi.ongpu){
        if(!psi.toGPU()){
          std::ostringstream err;
          err<<"Error when allocating GPU global memory.";
          throw std::runtime_error(exception_msg(err.str()));
        }
      }
      size_t iter = max_iter;
      if(!lanczosEV(ori_mat.m_elem, psi.m_elem, ori_mat.Rnum, iter, err_tol, E0, psi.m_elem, ori_mat.ongpu)){
        std::ostringstream err;
        err<<"Lanczos algorithm fails in converging.";;
        throw std::runtime_error(exception_msg(err.str()));
      }
      return iter;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::lanczosEigh(double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
      return 0;
    }
  }
  size_t lanczosEigh(cflag _tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol){
    try{
      if(!(ori_mat.Rnum == ori_mat.Cnum)){
        std::ostringstream err;
        err<<"Cannot perform Lanczos algorithm to find the lowest eigen value and eigen vector on a non-square matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(!(ori_mat.Rnum == psi.elemNum())){
        std::ostringstream err;
        err<<"Error in Lanczos initial vector psi. The vector dimension does not match with the number of the columns.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      if(ori_mat.ongpu && !psi.ongpu){
        if(!psi.toGPU()){
          std::ostringstream err;
          err<<"Error when allocating GPU global memory.";
          throw std::runtime_error(exception_msg(err.str()));
        }
      }
      size_t iter = max_iter;
      if(!lanczosEV(ori_mat.cm_elem, psi.cm_elem, ori_mat.Rnum, iter, err_tol, E0, psi.cm_elem, ori_mat.ongpu)){
        std::ostringstream err;
        err<<"Lanczos algorithm fails in converging.";;
        throw std::runtime_error(exception_msg(err.str()));
      }
      return iter;
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::lanczosEigh(double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
      return 0;
    }
  }

  size_t lanczosEigh(Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol){
    try{
      if(ori_mat.typeID() == 0){
        std::ostringstream err;
        err<<"Cannot perform lanczos decomposition on an EMPTY matrix.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      else if(ori_mat.typeID() == 1)
        return lanczosEigh(RTYPE, ori_mat, E0, psi, max_iter, err_tol);
      else if(ori_mat.typeID() == 2)
        return lanczosEigh(CTYPE, ori_mat, E0, psi, max_iter, err_tol);
    }
    catch(const std::exception& e){
      propogate_exception(e, "In function Matrix::lanczosEigh(double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
      return 0;
    }
  }

};	/* namespace uni10 */
#ifdef Block
#undef Block
#endif
#ifdef Matrix
#undef Matrix
#endif
#ifdef Real
#undef Real
#endif

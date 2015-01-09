/****************************************************************************
*  @file CMakeLists.txt
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
*  @brief Main specification file for CMake
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
namespace uni10 {
std::ostream& operator<< (std::ostream& os, const Matrix& m) {
    try {
        os << m.Rnum << " x " << m.Cnum << " = " << m.m_elemNum;
        if(m.diag)
            os << ", Diagonal";
        if(m.ongpu)
            os<< ", onGPU";
        os <<std::endl << std::endl;
        double* elem;
        if(m.ongpu) {
            elem = (double*)malloc(m.m_elemNum * sizeof(double));
            elemCopy(elem, m.m_elem, m.m_elemNum * sizeof(double), false, m.ongpu);
        }
        else
            elem = m.m_elem;
        for(size_t i = 0; i < m.Rnum; i++) {
            for(size_t j = 0; j < m.Cnum; j++)
                if(m.diag) {
                    if(i == j)
                        os << std::setw(7) << std::fixed << std::setprecision(3) << elem[i];
                    else
                        os << std::setw(7) << std::fixed << std::setprecision(3) << 0.0;
                }
                else
                    os << std::setw(7) << std::fixed << std::setprecision(3) << elem[i * m.Cnum + j];
            os << std::endl << std::endl;
        }
        if(m.ongpu)
            free(elem);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator<<(std::ostream&, uni10::Matrix&):");
    }
    return os;
}

Matrix::Matrix(): Rnum(0), Cnum(0), m_elemNum(0), diag(false), m_elem(NULL), ongpu(false) {
}
Matrix::Matrix(const Matrix& _m): Rnum(_m.Rnum), Cnum(_m.Cnum), m_elemNum(_m.m_elemNum), diag(_m.diag), m_elem(NULL), ongpu(false) {
    try {
        if(m_elemNum) {
            m_elem = (double*)elemAlloc(m_elemNum * sizeof(double), ongpu);
            elemCopy(m_elem, _m.m_elem, m_elemNum * sizeof(double), ongpu, _m.ongpu);
        }
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In copy constructor Matrix::Matrix(uni10::Matrix&):");
    }
}

void Matrix::init(bool _ongpu) {
    if(diag)
        m_elemNum = Rnum < Cnum ? Rnum : Cnum;
    if(m_elemNum) {
        if(_ongpu)  // Try to allocate GPU memory
            m_elem = (double*)elemAlloc(m_elemNum * sizeof(double), ongpu);
        else {
            m_elem = (double*)elemAllocForce(m_elemNum * sizeof(double), false);
            ongpu = false;
        }
    }
}

void Matrix::init(const double* _elem, bool src_ongpu) {
    init(true);
    elemCopy(m_elem, _elem, m_elemNum * sizeof(double), ongpu, src_ongpu);
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const double* _elem, bool _diag, bool src_ongpu): Rnum(_Rnum), Cnum(_Cnum), m_elemNum(_Rnum * _Cnum), diag(_diag), m_elem(NULL), ongpu(false) {
    try {
        init(_elem, src_ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, double*, bool=false):");
    }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, const std::vector<double>& _elem, bool _diag, bool src_ongpu): Rnum(_Rnum), Cnum(_Cnum), m_elemNum(_Rnum * _Cnum), diag(_diag), m_elem(NULL), ongpu(false) {
    try {
        init(&_elem[0], src_ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, std::vector<double>&, bool=false):");
    }
}

Matrix::Matrix(size_t _Rnum, size_t _Cnum, bool _diag, bool _ongpu): Rnum(_Rnum), Cnum(_Cnum), m_elemNum(_Rnum * _Cnum), diag(_diag), m_elem(NULL), ongpu(false) {
    try {
        init(_ongpu);
        if(m_elemNum)
            elemBzero(m_elem, m_elemNum * sizeof(double), ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In constructor Matrix::Matrix(size_t, size_t, bool=false):");
    }
}

Matrix& Matrix::operator=(const Matrix& _m) {
    try {
        Rnum = _m.Rnum;
        Cnum = _m.Cnum;
        m_elemNum = _m.m_elemNum;
        diag = _m.diag;
        if(m_elem != NULL)
            elemFree(m_elem, m_elemNum * sizeof(double), ongpu);
        m_elem = (double*)elemAlloc(m_elemNum * sizeof(double), ongpu);
        elemCopy(m_elem, _m.m_elem, m_elemNum * sizeof(double), ongpu, _m.ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::operator=(uni10::Matrix&):");
    }
    return *this;
}

Matrix::~Matrix() {
    try {
        if(m_elem != NULL)
            elemFree(m_elem, m_elemNum * sizeof(double), ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In destructor Matrix::~Matrix():");
    }
}

size_t Matrix::row()const {
    return Rnum;
}

size_t Matrix::col()const {
    return Cnum;
}
size_t Matrix::elemNum()const {
    return m_elemNum;
}

Matrix operator* (const Matrix& Ma, const Matrix& Mb) {
    try {
        if(!(Ma.Cnum == Mb.Rnum)) {
            std::ostringstream err;
            err<<"The dimensions of the two matrices do not match for matrix multiplication.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if((!Ma.diag) && (!Mb.diag)) {
            Matrix Mc(Ma.Rnum, Mb.Cnum);
            matrixMul(Ma.m_elem, Mb.m_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.m_elem, Ma.ongpu, Mb.ongpu, Mc.ongpu);
            return Mc;
        }
        else if(Ma.diag && (!Mb.diag)) {
            Matrix Mc(Mb);
            diagMM(Ma.m_elem, Mc.m_elem, Mc.Rnum, Mc.Cnum, Ma.ongpu, Mc.ongpu);
            return Mc;
        }
        else if((!Ma.diag) && Mb.diag) {
            Matrix Mc(Ma.Rnum, Mb.Cnum);
            for(size_t i = 0; i < Ma.Rnum; i++)
                for(size_t j = 0; j < Mb.m_elemNum; j++)
                    Mc.m_elem[i * Mb.Cnum + j] = Ma.m_elem[i * Ma.Cnum + j] * Mb.m_elem[j];
            return Mc;
        }
        else {
            Matrix Mc(Ma.Rnum, Mb.Cnum, true);
            for(size_t i = 0; i < Ma.Rnum; i++)
                Mc.m_elem[i] = Ma.m_elem[i] * Mb.m_elem[i];
            return Mc;
        }
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator*(uni10::Matrix&, uni10::Matrix&):");
        return Matrix();
    }
}

bool operator== (const Matrix& m1, const Matrix& m2) {
    try {
        double diff;
        if(m1.m_elemNum == m2.m_elemNum) {
            for(size_t i = 0; i < m1.m_elemNum; i++) {
                diff = fabs(m1.m_elem[i] - m2.m_elem[i]);
                if(diff > 1E-6)
                    return false;
            }
        }
        else
            return false;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator==(uni10::Matrix&, uni10::Matrix&):");
    }
    return true;
}


Matrix& Matrix::operator*= (const Matrix& Mb) {
    try {
        if(!ongpu)
            m_elem = (double*)mvGPU(m_elem, m_elemNum * sizeof(double), ongpu);
        *this = *this * Mb;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::operator*=(uni10::Matrix&):");
    }
    return *this;
}

void Matrix::setElem(const std::vector<double>& elem, bool _ongpu) {
    try {
        setElem(&elem[0], _ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::setElem(std::vector<double>&, bool=false):");
    }
}
void Matrix::setElem(const double* elem, bool _ongpu) {
    try {
        elemCopy(m_elem, elem, m_elemNum * sizeof(double), ongpu, _ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::setElem(double*, bool=false):");
    }
}

std::vector<Matrix> Matrix::eigh()const {
    std::vector<Matrix> outs;
    try {
        if(!(Rnum == Cnum)) {
            std::ostringstream err;
            err<<"Cannot perform eigenvalue decomposition on a non-square matrix.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if(diag) {
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
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::eigh():");
    }
    return outs;
}

std::vector<Matrix> Matrix::svd()const {
    std::vector<Matrix> outs;
    try {
        if(diag) {
            std::ostringstream err;
            err<<"Cannot perform singular value decomposition on a diagonal matrix. Need not to do so.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        size_t min = Rnum < Cnum ? Rnum : Cnum; //min = min(Rnum,Cnum)
        Matrix U(Rnum, min, false, ongpu);
        Matrix S(min, min, true, ongpu);
        Matrix VT(min, Cnum, false, ongpu);
        assert(U.isOngpu() == ongpu && VT.isOngpu() == ongpu);
        matrixSVD(m_elem, Rnum, Cnum, U.m_elem, S.m_elem, VT.m_elem, ongpu);
        outs.push_back(U);
        outs.push_back(S);
        outs.push_back(VT);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::svd():");
    }
    return outs;
}

size_t Matrix::lanczosEigh(double& E0, Matrix& psi, size_t max_iter, double err_tol) {
    try {
        if(!(Rnum == Cnum)) {
            std::ostringstream err;
            err<<"Cannot perform Lanczos algorithm to find the lowest eigen value and eigen vector on a non-square matrix.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if(!(Rnum == psi.elemNum())) {
            std::ostringstream err;
            err<<"Error in Lanczos initial vector psi. The vector dimension does not match with the number of the columns.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if(ongpu && !psi.ongpu)
            psi.m_elem = (double*)mvGPU(psi.m_elem, psi.m_elemNum * sizeof(double), psi.ongpu);
        size_t iter = max_iter;
        if(!lanczosEV(m_elem, psi.m_elem, Rnum, iter, err_tol, E0, psi.m_elem, ongpu)) {
            std::ostringstream err;
            err<<"Lanczos algorithm fails in converging.";;
            throw std::runtime_error(exception_msg(err.str()));
        }
        return iter;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::lanczosEigh(double& E0, uni10::Matrix&, size_t=200, double=5E-15):");
        return 0;
    }
}

void Matrix::randomize() {
    try {
        if(!ongpu)
            m_elem = (double*)mvGPU(m_elem, m_elemNum * sizeof(double), ongpu);
        elemRand(m_elem, m_elemNum, ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::randomize():");
    }
}


void Matrix::orthoRand() {
    try {
        if(!ongpu)
            m_elem = (double*)mvGPU(m_elem, m_elemNum * sizeof(double), ongpu);
        if(!diag) {
            orthoRandomize(m_elem, Rnum, Cnum, ongpu);
        }
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::orthoRand():");
    }
}

void Matrix::identity() {
    try {
        diag = true;
        m_elemNum = Rnum < Cnum ? Rnum : Cnum;
        if(m_elem != NULL)
            elemFree(m_elem, m_elemNum * sizeof(double), ongpu);
        m_elem = (double*)elemAlloc(m_elemNum * sizeof(double), ongpu);
        double* elemI = (double*)malloc(m_elemNum * sizeof(double));
        for(int i = 0; i < m_elemNum; i++)
            elemI[i] = 1;
        this->setElem(elemI, false);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::identity():");
    }
}

void Matrix::set_zero() {
    try {
        if(m_elemNum)
            elemBzero(m_elem, m_elemNum * sizeof(double), ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::set_zero():");
    }
}

Matrix operator*(const Matrix& Ma, double a) {
    try {
        Matrix Mb(Ma);
        vectorScal(a, Mb.m_elem, Mb.m_elemNum, Mb.ongpu);
        return Mb;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator*(uni10::Matrix&, double):");
        return Matrix();
    }
}

Matrix& Matrix::operator*= (double a) {
    try {
        if(!ongpu)
            m_elem = (double*)mvGPU(m_elem, m_elemNum * sizeof(double), ongpu);
        vectorScal(a, m_elem, m_elemNum, ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::operator*=(double):");
    }
    return *this;
}

Matrix operator+(const Matrix& Ma, const Matrix& Mb) {
    try {
        Matrix Mc(Ma);
        vectorAdd(Mc.m_elem, Mb.m_elem, Mc.m_elemNum, Mc.ongpu, Mb.ongpu);
        return Mc;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function operator+(uni10::Matrix&, uni10::Matrix&):");
        return Matrix();
    }
}

Matrix& Matrix::operator+= (const Matrix& Mb) {
    try {
        if(!ongpu)
            m_elem = (double*)mvGPU(m_elem, m_elemNum * sizeof(double), ongpu);
        vectorAdd(m_elem, Mb.m_elem, m_elemNum, ongpu, Mb.ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::operator+=(uni10::Matrix&):");
    }
    return *this;
}

Matrix& Matrix::transpose() {
    try {
        if(!ongpu)
            m_elem = (double*)mvGPU(m_elem, m_elemNum * sizeof(double), ongpu);
        if(!diag) {
            double* transElem;
            size_t memsize = m_elemNum * sizeof(double);
            transElem = (double*)elemAllocForce(memsize, ongpu);
            setTranspose(m_elem, Rnum, Cnum, transElem, ongpu);
            if(m_elem != NULL)
                elemFree(m_elem, memsize, ongpu);
            m_elem = transElem;
        }
        size_t tmp = Rnum;
        Rnum = Cnum;
        Cnum = tmp;
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::transpose():");
    }
    return *this;
}

Matrix& Matrix::resize(size_t row, size_t col) {
    try {
        if(diag) {
            size_t elemNum = row < col ? row : col;
            if(elemNum > m_elemNum) {
                bool des_ongpu;
                double* elem = (double*)elemAlloc(elemNum * sizeof(double), des_ongpu);
                elemBzero(elem, elemNum * sizeof(double), des_ongpu);
                elemCopy(elem, m_elem, m_elemNum * sizeof(double), des_ongpu, ongpu);
                if(m_elem != NULL)
                    elemFree(m_elem, m_elemNum * sizeof(double), ongpu);
                m_elem = elem;
                ongpu = des_ongpu;
            }
            else
                shrinkWithoutFree((m_elemNum - elemNum) * sizeof(double), ongpu);
            Rnum = row;
            Cnum = col;
            m_elemNum = elemNum;
        }
        else {
            if(col == Cnum) {
                size_t elemNum = row * col;
                if(row > Rnum) {
                    bool des_ongpu;
                    double* elem = (double*)elemAlloc(elemNum * sizeof(double), des_ongpu);
                    elemBzero(elem, elemNum * sizeof(double), des_ongpu);
                    elemCopy(elem, m_elem, m_elemNum * sizeof(double), des_ongpu, ongpu);
                    if(m_elem != NULL)
                        elemFree(m_elem, m_elemNum * sizeof(double), ongpu);
                    m_elem = elem;
                    ongpu = des_ongpu;
                }
                else
                    shrinkWithoutFree((m_elemNum - elemNum) * sizeof(double), ongpu);
                Rnum = row;
                m_elemNum = elemNum;
            }
            else {
                size_t data_row = row < Rnum ? row : Rnum;
                size_t data_col = col < Cnum ? col : Cnum;
                bool des_ongpu;
                double* elem = (double*)elemAlloc(row * col * sizeof(double), des_ongpu);
                elemBzero(elem, row * col * sizeof(double), des_ongpu);
                for(size_t r = 0; r < data_row; r++)
                    elemCopy(&(elem[r * col]), &(m_elem[r * Cnum]), data_col * sizeof(double), des_ongpu, ongpu);
                if(m_elem != NULL)
                    elemFree(m_elem, m_elemNum * sizeof(double), ongpu);
                m_elem = elem;
                ongpu = des_ongpu;
                Rnum = row;
                Cnum = col;
                m_elemNum = row * col;
            }
        }
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::resize(size_t, size_t):");
    }
    return *this;
}

double Matrix::norm() {
    try {
        return vectorNorm(m_elem, m_elemNum, 1, ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::norm():");
        return 0;
    }
}
double Matrix::sum() {
    try {
        return vectorSum(m_elem, m_elemNum, 1, ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::sum():");
        return 0;
    }
}
double Matrix::trace() {
    try {
        if(!(Rnum == Cnum)) {
            std::ostringstream err;
            err<<"Cannot perform trace on a non-square matrix.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        if(diag)
            return vectorSum(m_elem, m_elemNum, 1, ongpu);
        else
            return vectorSum(m_elem, Rnum, Cnum + 1, ongpu);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::trace():");
        return 0;
    }
}
void Matrix::save(const std::string& fname) {
    try {
        FILE *fp = fopen(fname.c_str(), "w");
        if(!(fp != NULL)) {
            std::ostringstream err;
            err<<"Error in writing to file '"<<fname<<"'.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        double* elem = m_elem;
        if(ongpu) {
            elem = (double*)malloc(m_elemNum * sizeof(double));
            elemCopy(elem, m_elem, m_elemNum * sizeof(double), false, ongpu);
        }
        fwrite(elem, sizeof(double), m_elemNum, fp);
        fclose(fp);
        if(ongpu)
            free(elem);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::save(std::string&):");
    }
}

void Matrix::load(const std::string& fname) {
    try {
        FILE *fp = fopen(fname.c_str(), "r");
        if(!(fp != NULL)) {
            std::ostringstream err;
            err<<"Error in opening file '" << fname <<"'.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        double* elem = m_elem;
        if(ongpu)
            elem = (double*)malloc(m_elemNum * sizeof(double));
        fread(elem, sizeof(double), m_elemNum, fp);
        fclose(fp);
        if(ongpu) {
            elemCopy(m_elem, elem, m_elemNum * sizeof(double), ongpu, false);
            free(elem);
        }
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::load(std::string&):");
    }
}

double& Matrix::operator[](size_t idx) {
    try {
        if(!(idx < m_elemNum)) {
            std::ostringstream err;
            err<<"Index exceeds the number of the matrix elements("<<m_elemNum<<").";
            throw std::runtime_error(exception_msg(err.str()));
        }
        m_elem = (double*)mvCPU(m_elem, m_elemNum * sizeof(double), ongpu);
        return m_elem[idx];
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::opeartor[](size_t):");
        return m_elem[0];
    }
}

double* Matrix::getElem()const {
    return m_elem;
}

double* Matrix::getHostElem() {
    try {
        if(ongpu) {
            m_elem = (double*)mvCPU(m_elem, m_elemNum * sizeof(double), ongpu);
        }
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::getHostElem():");
    }
    return m_elem;
}

double& Matrix::at(size_t r, size_t c) {
    try {
        if(!((r < Rnum) && (c < Cnum))) {
            std::ostringstream err;
            err<<"The input indices are out of range.";
            throw std::runtime_error(exception_msg(err.str()));
        }
        m_elem = (double*)mvCPU(m_elem, m_elemNum * sizeof(double), ongpu);
        if(diag) {
            if(!(r == c && r < m_elemNum)) {
                std::ostringstream err;
                err<<"The matrix is diagonal, there is no off-diagonal element.";
                throw std::runtime_error(exception_msg(err.str()));
            }
            return m_elem[r];
        }
        else
            return m_elem[r * Cnum + c];
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function Matrix::at(size_t, size_t):");
        return m_elem[0];
    }
}

bool Matrix::toGPU() {
    if(!ongpu)
        m_elem = (double*)mvGPU(m_elem, m_elemNum * sizeof(double), ongpu);
    return ongpu;
}

Matrix takeExp(double a, const Matrix& mat) {
    try {
        std::vector<Matrix> rets = mat.eigh();
        Matrix UT(rets[1]);
        UT.transpose();
        vectorExp(a, rets[0].getElem(), rets[0].row(), rets[0].isOngpu());
        return UT * (rets[0] * rets[1]);
    }
    catch(const std::exception& e) {
        propogate_exception(e, "In function takeExp(double, uni10::Matrix&):");
        return Matrix();
    }
}

};  /* namespace uni10 */

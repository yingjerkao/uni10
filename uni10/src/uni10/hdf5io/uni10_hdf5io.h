/****************************************************************************
*  @file uni10_lapack.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2016
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
*  @brief Header file for hdf5io class
*  @author Chen-Yen Lai
*  @date 2016-01-23
*  @since 0.9.2
*
*****************************************************************************/
#ifndef UNI10_HDF5IO_H
#define UNI10_HDF5IO_H
#include <string>
#include <complex>
#include <vector>
#include <H5Cpp.h>
#include <uni10/datatype/Qnum.h>
#include <uni10/data-structure/Bond.h>
#include <uni10/data-structure/Block.h>

namespace uni10{

typedef struct uni10_qnum_hdf5 {
    int m_U1;
    parityType m_prt;
    parityFType m_prtF;
} uni10_qnum_hdf5;

typedef struct uni10_rcflag_hdf5 {
    rflag r_flag;
    cflag c_flag;
} uni10_rcflag_hdf5;

class HDF5IO: public H5::H5File {
private:
    std::string FileName;
    H5::Group Uni10TypeGroup;
    H5::CompType ComplexDataType;
    const H5::CompType initCompexDataType(void);
    const H5::EnumType initParityEnumType(void);
    const H5::EnumType initParityFEnumType(void);
    H5::CompType QnumCompType;
    const H5::CompType initQnumCompType(void);
    H5::EnumType BondEnumType;
    const H5::EnumType initBondEnumType(void);
    const H5::EnumType initBlockRflagEnumType(void);
    const H5::EnumType initBlockCflagEnumType(void);
    H5::CompType RCflagType;
    const H5::CompType initRCflagType(void);
    bool fileExists(const std::string &FileName);

public:
    HDF5IO (const std::string &FileName);
    HDF5IO (const std::string &FileName, const bool force);
    virtual ~HDF5IO (void);
    inline H5::H5File getFile(){return *this;};
    H5::Group getGroup(const std::string &GroupName);

    void saveNumber(const std::string& GroupName, const std::string& Name, int x);
    void saveNumber(const std::string& GroupName, const std::string& Name, unsigned long x);
    void saveNumber(const std::string& GroupName, const std::string& Name, double x);
    void saveNumber(const std::string& GroupName, const std::string& Name, std::complex<double> C);

    void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<int>& V);
    void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<size_t>& V);
    void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<double>& V);
    void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<std::complex<double> >& V);

    void saveRawBuffer(const std::string& GroupName, const std::string& Name, const size_t dim, const int* x);
    void saveRawBuffer(const std::string& GroupName, const std::string& Name, const size_t dim, const double* x);
    void saveRawBuffer(const std::string& GroupName, const std::string& Name, const size_t dim, const std::complex<double>* x);

    int loadInt(const std::string& GroupName, const std::string& Name);
    size_t loadUlong(const std::string& GroupName, const std::string& Name);
    double loadReal(const std::string& GroupName, const std::string& Name);
    std::complex<double> loadComplex(const std::string& GroupName, const std::string& Name);
    void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<int>& V);
    void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<size_t>& V);
    void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<double>& V);
    void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<std::complex<double> >& V);

    void loadRawBuffer(const std::string& GroupName, const std::string& Name, size_t& dim, int*& x);
    void loadRawBuffer(const std::string& GroupName, const std::string& Name, size_t& dim, double*& x);
    void loadRawBuffer(const std::string& GroupName, const std::string& Name, size_t& dim, std::complex<double>*& x);

    // uni10 properties
    void saveBond(const std::string& GroupName, const std::string& Name, const bondType& _bt);
    void saveFlag(const std::string& GroupName, const std::string& Name, const rflag& _rf, const cflag& _cf);
    void saveQnum(const std::string& GroupName, const std::string& Name, const int &_U1, const parityType &_pt, const parityFType &_ptf);
    void loadBond(const std::string& GroupName, const std::string& Name, bondType& _bt);
    void loadFlag(const std::string& GroupName, const std::string& Name, rflag& _rf, cflag& _cf);
    void loadQnum(const std::string& GroupName, const std::string& Name, int &_U1, parityType &_pt, parityFType &_ptf);
};

}; /* end of namespace uni10 */
#endif  /* end of include guard: UNI10_HDF5IO_H */

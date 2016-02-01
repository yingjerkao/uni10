/****************************************************************************
*  @file CMakeLists.txt
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
*  @brief Implementation file of hdf5io class
*  @author Chen-Yen Lai
*  @date 2016-01-23
*  @since 0.1.0
*
*****************************************************************************/
#include <iostream>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <cstdlib>
#include <uni10/hdf5io/uni10_hdf5io.h>

namespace uni10{

/* For the HDF5IO class */
const H5::CompType HDF5IO::initCompexDataType(){
    H5::CompType Type(sizeof(std::complex<double>));
    // NOTE: New HDF5 with native complex datatypes ??
    Type.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
    Type.insertMember("imag", sizeof(double), H5::PredType::NATIVE_DOUBLE);
    try {
        Type.commit(*this, "uni10Types/complex");
    } catch(H5::DataTypeIException){};
    return Type;
}

const H5::EnumType HDF5IO::initParityEnumType(){
    H5::EnumType Type( sizeof(int) );
    parityType pt = PRT_EVEN;
    Type.insert("PRT_EVEN", &pt);
    pt = PRT_ODD;
    Type.insert("PRT_ODD", &pt);
    try {
        Type.commit(*this, "uni10Types/parityType");
    } catch(H5::DataTypeIException){};
    return Type;
}

const H5::EnumType HDF5IO::initParityFEnumType(){
    H5::EnumType Type( sizeof(int) );
    parityFType pt = PRTF_EVEN;
    Type.insert("PRTF_EVEN", &pt);
    pt = PRTF_ODD;
    Type.insert("PRTF_ODD", &pt);
    try {
        Type.commit(*this, "uni10Types/parityFType");
    } catch(H5::DataTypeIException){};
    return Type;
}

const H5::CompType HDF5IO::initQnumCompType(){
    H5::CompType Type( sizeof(uni10_qnum_hdf5) );
    H5::EnumType ParityEnumType = initParityEnumType();
    H5::EnumType ParityFEnumType = initParityFEnumType();
    Type.insertMember("U1", HOFFSET(uni10_qnum_hdf5, m_U1), H5::PredType::NATIVE_INT);
    Type.insertMember("parity", HOFFSET(uni10_qnum_hdf5, m_prt), ParityEnumType);
    Type.insertMember("parityF", HOFFSET(uni10_qnum_hdf5, m_prtF), ParityFEnumType);
    try {
        Type.commit(*this, "uni10Types/qnum");
    } catch(H5::DataTypeIException){};
    return Type;
}

/*NOTE: Can I combine all these enum types initilization into single function?
        Can I build a compond type for bond and qnum?
*/
const H5::EnumType HDF5IO::initBondEnumType(){
    H5::EnumType Type( sizeof(int) );
    bondType bond = BD_IN;
    Type.insert("In", &bond);
    bond = BD_OUT;
    Type.insert("Out", &bond);
    try {
        Type.commit(*this, "uni10Types/bondType");
    } catch(H5::DataTypeIException){};
    return Type;
}

const H5::EnumType HDF5IO::initBlockRflagEnumType(){
    H5::EnumType Type( sizeof(int) );
    rflag rf = RNULL;
    Type.insert("RNULL", &rf);
    rf = RTYPE;
    Type.insert("RTYPE", &rf);
    try {
        Type.commit(*this, "uni10Types/rflag");
    } catch(H5::DataTypeIException){};
    return Type;
}

const H5::EnumType HDF5IO::initBlockCflagEnumType(){
    H5::EnumType Type( sizeof(int) );
    cflag cf = CNULL;
    Type.insert("CNULL", &cf);
    cf = CTYPE;
    Type.insert("CTYPE", &cf);
    try {
        Type.commit(*this, "uni10Types/cflag");
    } catch(H5::DataTypeIException){};
    return Type;
}

const H5::CompType HDF5IO::initRCflagType(){
    H5::CompType Type( sizeof(uni10_rcflag_hdf5) );
    H5::EnumType BlockRflagEnumType = initBlockRflagEnumType();
    H5::EnumType BlockCflagEnumType = initBlockCflagEnumType();
    Type.insertMember("rflag", HOFFSET(uni10_rcflag_hdf5, r_flag), BlockRflagEnumType);
    Type.insertMember("cflag", HOFFSET(uni10_rcflag_hdf5, c_flag), BlockCflagEnumType);
    try {
        Type.commit(*this, "uni10Types/flags");
    } catch(H5::DataTypeIException){};
    return Type;
}

bool HDF5IO::fileExists(const std::string& FileName){
  try {
      H5::Exception::dontPrint();
      return isHdf5(FileName.c_str());
  } catch(H5::FileIException){
      return false;
  }
}

HDF5IO::HDF5IO (const std::string& fn) :
FileName(fn),
H5File(fn.c_str(), fileExists(fn) ? H5F_ACC_RDWR : H5F_ACC_TRUNC),
Uni10TypeGroup(getGroup("uni10Types")),
ComplexDataType(initCompexDataType()),
BondEnumType(initBondEnumType()),
RCflagType(initRCflagType()),
QnumCompType(initQnumCompType()){
    std::cout << "Open Hdf5 file - " << FileName << std::endl;
}

HDF5IO::HDF5IO (const std::string& fn, const bool force):
FileName(fn),
H5File(fn.c_str(), H5F_ACC_TRUNC),
Uni10TypeGroup(getGroup("uni10Types")),
ComplexDataType(initCompexDataType()),
BondEnumType(initBondEnumType()),
RCflagType(initRCflagType()),
QnumCompType(initQnumCompType()){
    std::cout << "Force opend HDF5 file - " << FileName << std::endl;
}

HDF5IO::~HDF5IO(){
    std::cout << "Close HDF5 file - " << this->getFileName() << std:: endl;
    this->close();
}

H5::Group HDF5IO::getGroup(const std::string &GroupName){
    H5::Group group;
    try{
        group = H5::Group( this->openGroup( GroupName.c_str() ) );
    } catch( H5::FileIException not_found_error ){
        group = H5::Group( this->createGroup( GroupName.c_str() ) );
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << std::endl;
        throw std::runtime_error("HDF5IO::getGroup");
    }
    return group;
}

/* Save/Load single value data */
void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    int x){
    H5::Group FG = getGroup( GroupName );
    try{
        H5::Exception::dontPrint();
        H5::DataSet dataset = FG.openDataSet( Name.c_str() );
        dataset.write(&x, H5::PredType::NATIVE_INT);
    } catch ( const H5::GroupIException not_found_error ){
        H5::DataSet dataset = FG.createDataSet( Name.c_str(),
                  H5::PredType::NATIVE_INT, H5::DataSpace());
        dataset.write(&x, H5::PredType::NATIVE_INT);
    }
    FG.close();
}

int HDF5IO::loadInt(const std::string& GroupName, const std::string& Name){
    try{
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet( Name.c_str());
        int x;
        DataSet.read(&x,H5::PredType::NATIVE_INT);
        FG.close();
        return x;
    }catch( H5::GroupIException not_found_error ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("No dataset found in loadInt. ");
    }
}

void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    unsigned long x){
    H5::Group FG = getGroup( GroupName );
    try{
        H5::Exception::dontPrint();
        H5::DataSet dataset = FG.openDataSet( Name.c_str() );
        dataset.write(&x, H5::PredType::NATIVE_ULONG);
    } catch ( const H5::GroupIException not_found_error ){
        H5::DataSet dataset = FG.createDataSet( Name.c_str(),
                H5::PredType::NATIVE_ULONG, H5::DataSpace());
        dataset.write(&x, H5::PredType::NATIVE_ULONG);
    }
    FG.close();
}

size_t HDF5IO::loadUlong(const std::string& GroupName, const std::string& Name){
  try{
      H5::Group FG = getGroup( GroupName );
      H5::DataSet DataSet = FG.openDataSet( Name.c_str() );
      size_t x;
      DataSet.read(&x, H5::PredType::NATIVE_ULONG);
      return x;
      FG.close();
  }catch( H5::GroupIException not_found_error ){
      std::cout << "In Group - " << GroupName << ", and Name is "
                << Name << std::endl;
      throw std::runtime_error("No dataset found in loadUlong. ");
  }
}

// Save/load double
void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    double x){
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dataset = FG.openDataSet(Name.c_str());
      dataset.write(&x,H5::PredType::NATIVE_DOUBLE);
    } catch ( const H5::GroupIException not_found_error ) {
      H5::DataSet dataset = FG.createDataSet(Name.c_str(),
                  H5::PredType::NATIVE_DOUBLE,H5::DataSpace());
      dataset.write(&x,H5::PredType::NATIVE_DOUBLE);
    }
    FG.close();
}

double HDF5IO::loadReal(const std::string& GroupName, const std::string& Name){
  try{
      H5::Group FG = getGroup( GroupName );
      H5::DataSet DataSet = FG.openDataSet(Name.c_str());
      double x;
      DataSet.read(&x,H5::PredType::NATIVE_DOUBLE);
      FG.close();
      return x;
  }catch( H5::GroupIException not_found_error ){
      std::cout << "In Group - " << GroupName << ", and Name is "
                << Name << std::endl;
      throw std::runtime_error("No dataset found in loadReal. ");
  }
}

void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    std::complex<double> C){
    // H5::CompType ComplexDataType = openCompType("uni10Types/complex");
    H5::Group FG = getGroup( GroupName );
    try{
        H5::Exception::dontPrint();
        H5::DataSet dataset = FG.openDataSet(Name.c_str());
        double RealImag[2] = {real(C),imag(C)};
        dataset.write(RealImag, ComplexDataType);
    } catch ( const H5::GroupIException not_found_error ){
        H5::DataSet dataset = FG.createDataSet(Name.c_str(), ComplexDataType,
                                 H5::DataSpace());
        double RealImag[2] = {real(C),imag(C)};
        dataset.write(RealImag, ComplexDataType);
    }
    FG.close();
}

std::complex<double> HDF5IO::loadComplex(const std::string& GroupName,
    const std::string& Name){
    try{
        // H5::CompType ComplexDataType = this->openCompType("uni10Types/complex");
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        std::complex<double> C;
        double RealImag[2];
        DataSet.read(RealImag, ComplexDataType);
        FG.close();
        return std::complex<double>(RealImag[0],RealImag[1]);
    }catch( H5::GroupIException not_found_error ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("No dataset found in loadComplex. ");
    }
}

/* Load/Save std::vector */
void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<int>& V){
    try{
      hsize_t Dim[1] = {hsize_t(V.size())};
      H5::DataSpace dspace(1,Dim);
      H5::Group FG = getGroup( GroupName );
      try{
        H5::Exception::dontPrint();
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        DataSet.write(V.data(),H5::PredType::NATIVE_INT, dspace);
      } catch ( const H5::GroupIException not_found_error ){
        H5::DataSet DataSet = FG.createDataSet(Name.c_str(),
          H5::PredType::NATIVE_INT, dspace);
        DataSet.write(V.data(),H5::PredType::NATIVE_INT);
      }
      FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveIntStdVector ");
    }
}

void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<size_t>& V){
    try{
        hsize_t Dim[1] = {hsize_t(V.size())};
        H5::DataSpace dspace(1,Dim);
        H5::Group FG = getGroup( GroupName );
        try{
            H5::Exception::dontPrint();
            H5::DataSet DataSet = FG.openDataSet(Name.c_str());
            DataSet.write(V.data(),H5::PredType::NATIVE_ULONG, dspace);
        } catch ( const H5::GroupIException not_found_error ){
            H5::DataSet DataSet = FG.createDataSet(Name.c_str(),
                                  H5::PredType::NATIVE_ULONG, dspace);
            DataSet.write(V.data(), H5::PredType::NATIVE_ULONG);
        }
        FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveUlongStdVector ");
    }
}

void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<double>& V){
    try{
        hsize_t Dim[1] = {hsize_t(V.size())};
        H5::DataSpace dataspace(1,Dim);
        H5::Group FG = getGroup( GroupName );
        try{
            H5::Exception::dontPrint();
            H5::DataSet dataset = FG.openDataSet(Name.c_str());
            dataset.write(V.data(),H5::PredType::NATIVE_DOUBLE, dataspace);
        } catch ( const H5::GroupIException not_found_error ){
            H5::DataSet dataset = FG.createDataSet(Name.c_str(),
                                  H5::PredType::NATIVE_DOUBLE,dataspace);
            dataset.write(V.data(),H5::PredType::NATIVE_DOUBLE);
        }
        FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveRealStdVector ");
    }
}

void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<std::complex<double> >& V){
    try{
        // H5::CompType ComplexDataType = this->openCompType("uni10Types/complex");
        hsize_t Dim[1] = {hsize_t(V.size())};
        H5::DataSpace dataspace(1,Dim);
        H5::Group FG = getGroup( GroupName.c_str() );
        try{
            H5::Exception::dontPrint();
            H5::DataSet dataset = FG.openDataSet(Name.c_str());
            dataset.write(V.data(), ComplexDataType, dataspace);
        } catch( const H5::GroupIException not_found_error ){
            H5::DataSet dataset = FG.createDataSet(Name.c_str(),
                                     ComplexDataType, dataspace);
            dataset.write(V.data(), ComplexDataType);
        } catch( const H5::FileIException error){
            error.printError();
        } catch( const H5::DataSetIException error){
            error.printError();
        }
        FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveComplexStdVector. ");
    }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<size_t>& V){
    try{
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        H5::DataSpace DataSpace = DataSet.getSpace();
        if(DataSpace.getSimpleExtentNdims() != 1)
            throw(H5::DataSpaceIException("HDF5IO::loadRealVector() ",
                          "Unexpected multidimentional dataspace. "));
        V.resize(DataSpace.getSimpleExtentNpoints());
        DataSet.read(V.data(),H5::PredType::NATIVE_ULONG);
        FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadUlongStdVector ");
    }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<int>& V){
    try{
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        H5::DataSpace DataSpace = DataSet.getSpace();
        if(DataSpace.getSimpleExtentNdims() != 1)
            throw(H5::DataSpaceIException("HDF5IO::loadRealVector() ",
                          "Unexpected multidimentional dataspace. "));
        V.resize(DataSpace.getSimpleExtentNpoints());
        DataSet.read(V.data(),H5::PredType::NATIVE_INT);
        FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadIntStdVector ");
    }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<double>& V){
    try{
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        H5::DataSpace DataSpace = DataSet.getSpace();
        if(DataSpace.getSimpleExtentNdims() != 1)
            throw(H5::DataSpaceIException("HDF5IO::loadRealVector()",
                          "Unexpected multidimentional dataspace."));
        V.resize(DataSpace.getSimpleExtentNpoints());
        DataSet.read(V.data(),H5::PredType::NATIVE_DOUBLE);
        FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadRealStdVector");
    }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<std::complex<double> >& V){
    try{
        // H5::CompType Type = this->openCompType("uni10Types/complex");
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        H5::DataSpace DataSpace = DataSet.getSpace();
        if(DataSpace.getSimpleExtentNdims() != 1)
          throw(H5::DataSpaceIException("HDF5IO::loadComplexVector()",
                            "Unexpected multidimentional dataspace."));
        V.resize(DataSpace.getSimpleExtentNpoints());
        DataSet.read(V.data(), ComplexDataType);
        FG.close();
    } catch( const H5::Exception err ){
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadComplexStdVector ");
    }
}

/* Load/Save C/C++ raw buffer */
void HDF5IO::saveRawBuffer(const std::string& GroupName, const std::string& Name,
    const size_t dim, const int* x){
    try{
        H5::Exception::dontPrint();
        hsize_t Dim[1] = {hsize_t(dim)};
        H5::DataSpace dspace(1,Dim);
        H5::Group FG = getGroup( GroupName );
        try{
            H5::DataSet dset = FG.openDataSet(Name.c_str());
            dset.write(x, H5::PredType::NATIVE_INT, dspace);
        } catch ( const H5::GroupIException not_found_error ){
            H5::DataSet dset = FG.createDataSet(Name.c_str(),
                           H5::PredType::NATIVE_DOUBLE, dspace);
            dset.write(x, H5::PredType::NATIVE_INT);
        }
        FG.close();
    }catch ( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO:: ");
    }
}

void HDF5IO::saveRawBuffer(const std::string& GroupName, const std::string& Name,
    const size_t dim, const double* x){
    try{
        H5::Exception::dontPrint();
        hsize_t Dim[1] = {hsize_t(dim)};
        H5::DataSpace dspace(1,Dim);
        H5::Group FG = getGroup( GroupName );
        try{
            H5::DataSet dset = FG.openDataSet(Name.c_str());
            dset.write(x, H5::PredType::NATIVE_DOUBLE, dspace);
        } catch ( const H5::GroupIException not_found_error ){
            H5::DataSet dset = FG.createDataSet(Name.c_str(),
                           H5::PredType::NATIVE_DOUBLE, dspace);
            dset.write(x, H5::PredType::NATIVE_DOUBLE);
        }
        FG.close();
    }catch ( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO:: ");
    }
}

void HDF5IO::saveRawBuffer(const std::string& GroupName, const std::string& Name,
    const size_t dim, const std::complex<double>* x){
    try{
        // H5::CompType Type = this->openCompType("uni10Types/complex");
        H5::Exception::dontPrint();
        hsize_t Dim[1] = {hsize_t(dim)};
        H5::DataSpace dspace(1,Dim);
        H5::Group FG = getGroup( GroupName );
        try{
            H5::DataSet dset = FG.openDataSet(Name.c_str());
            dset.write(x, ComplexDataType, dspace);
        } catch ( const H5::GroupIException not_found_error ){
            H5::DataSet dset = FG.createDataSet(Name.c_str(), ComplexDataType, dspace);
            dset.write(x, ComplexDataType);
        }
        FG.close();
    }catch ( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveRawBuffer ");
    }
}

void HDF5IO::loadRawBuffer(const std::string& GroupName, const std::string& Name,
    size_t& dim, int*& x){
    try{
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        H5::DataSpace DataSpace = DataSet.getSpace();
        if( DataSpace.getSimpleExtentNdims() != 1)
            throw(H5::DataSpaceIException("HDF5IO::loadRawBuffer",
                  "A dataspace must be precisely one-dimensional."));
        hsize_t Dims[1];
        DataSpace.getSimpleExtentDims(Dims);
        x = (int*)malloc( Dims[0] * sizeof(int));
        dim = Dims[0];
        DataSet.read(x, H5::PredType::NATIVE_INT);
        FG.close();
    }catch ( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadRawBuffer ");
    }
}

void HDF5IO::loadRawBuffer(const std::string& GroupName, const std::string& Name,
    size_t& dim, double*& x){
    try{
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        H5::DataSpace DataSpace = DataSet.getSpace();
        if( DataSpace.getSimpleExtentNdims() != 1)
            throw(H5::DataSpaceIException("HDF5IO::loadRawBuffer",
                  "A dataspace must be precisely one-dimensional."));
        hsize_t Dims[1];
        DataSpace.getSimpleExtentDims(Dims);
        x = (double*)malloc( Dims[0] * sizeof(double));
        dim = Dims[0];
        DataSet.read(x, H5::PredType::NATIVE_DOUBLE);
        FG.close();
    }catch ( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadRawBuffer ");
    }
}

void HDF5IO::loadRawBuffer(const std::string& GroupName, const std::string& Name,
    size_t& dim, std::complex<double>*& x){
    try{
        // H5::CompType Type = this->openCompType("uni10Types/complex");
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet(Name.c_str());
        H5::DataSpace DataSpace = DataSet.getSpace();
        if( DataSpace.getSimpleExtentNdims() != 1)
            throw(H5::DataSpaceIException("HDF5IO::loadRawBuffer",
                  "A dataspace must be precisely one-dimensional."));
        hsize_t Dims[1];
        DataSpace.getSimpleExtentDims(Dims);
        x = (std::complex<double>*)malloc(Dims[0] * sizeof(std::complex<double>));
        dim = Dims[0];
        DataSet.read(x, ComplexDataType);
        FG.close();
    }catch ( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadRawBuffer ");
    }
}

/* Load/Save uni10 properties */
void HDF5IO::saveBond(const std::string& GroupName, const std::string& Name, const bondType& _bt){
    // H5::EnumType Type = this->openEnumType("uni10Types/bondType");
    H5::Group FG = getGroup( GroupName );
    try{
        H5::Exception::dontPrint();
        H5::DataSet dset = FG.openDataSet(Name.c_str());
        dset.write(&_bt, BondEnumType);
    }catch( H5::GroupIException not_found_error ) {
        H5::DataSet dset = FG.createDataSet(Name.c_str(), BondEnumType, H5::DataSpace());
        dset.write(&_bt, BondEnumType);
    }catch( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveBond ");
    }
}
void HDF5IO::saveFlag(const std::string& GroupName, const std::string& Name,
    const rflag& _rf, const cflag& _cf){
    // H5::CompType Type = this->openCompType("uni10Types/flags");
    uni10_rcflag_hdf5 _input;
    _input.r_flag = _rf;
    _input.c_flag = _cf;
    H5::Group FG = getGroup( GroupName );
    try{
        H5::Exception::dontPrint();
        H5::DataSet dset = FG.openDataSet(Name.c_str());
        dset.write(&_input, RCflagType);
    }catch( H5::GroupIException not_found_error ) {
        H5::DataSet dset = FG.createDataSet(Name.c_str(), RCflagType, H5::DataSpace());
        dset.write(&_input, RCflagType);
    }catch( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveFlag ");
    }
}
void HDF5IO::saveQnum(const std::string& GroupName, const std::string& Name,
    const int &_U1, const parityType &_pt, const parityFType &_ptf) {
    uni10_qnum_hdf5 _input;
    _input.m_U1 = _U1;
    _input.m_prt = _pt;
    _input.m_prtF = _ptf;
    H5::Group FG = getGroup( GroupName );
    try{
        H5::Exception::dontPrint();
        H5::DataSet dset = FG.openDataSet(Name.c_str());
        dset.write(&_input, QnumCompType);
    }catch( H5::GroupIException not_found_error ) {
        H5::DataSet dset = FG.createDataSet(Name.c_str(), QnumCompType, H5::DataSpace());
        dset.write(&_input, QnumCompType);
    }catch( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::saveQnum ");
    }
}

void HDF5IO::loadBond(const std::string& GroupName, const std::string& Name, bondType& _bt){
    try{
        // H5::EnumType Type = this->openEnumType("uni10Types/bondType");
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet( Name.c_str() );
        DataSet.read(&_bt, BondEnumType);
        FG.close();
    }catch( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadBond ");
    }
}
void HDF5IO::loadFlag(const std::string& GroupName, const std::string& Name, rflag& _rf, cflag& _cf){
    try{
        // H5::CompType Type = this->openEnumType("uni10Types/flags");
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet( Name.c_str() );
        uni10_rcflag_hdf5 _output;
        DataSet.read(&_output, RCflagType);
        _rf = _output.r_flag;
        _cf = _output.c_flag;
        FG.close();
    }catch( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadFlag ");
    }
}
void HDF5IO::loadQnum(const std::string& GroupName, const std::string& Name,
    int &_U1, parityType &_pt, parityFType &_ptf){
    try{
        H5::Group FG = getGroup( GroupName );
        H5::DataSet DataSet = FG.openDataSet( Name.c_str() );
        uni10_qnum_hdf5 _output;
        DataSet.read(&_output, QnumCompType);
        _U1 = _output.m_U1;
        _pt = _output.m_prt;
        _ptf = _output.m_prtF;
        FG.close();
    }catch( const H5::Exception ) {
        std::cout << "In Group - " << GroupName << ", and Name is "
                  << Name << std::endl;
        throw std::runtime_error("HDF5IO::loadQnum ");
    }
}

}; /* end of namespace uni10 */

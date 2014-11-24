%module pyUni10
%{
  /* Put header files here or function declarations like below */
  #include <sstream>
  #include <uni10/datatype/Qnum.h>
  #include <uni10/data-structure/Bond.h>
  #include <uni10/tensor-network/Matrix.h>
  #include <uni10/tensor-network/UniTensor.h>
  #include <uni10/tensor-network/Network.h>
%}

%begin %{
#ifdef _MSC_VER
#define SWIG_PYTHON_INTERPRETER_NO_DEBUG
#endif
%}


%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "exception.i"
/*%include "typemaps.i"*/
namespace std{
  %template(int_arr) vector<int>;
  %template(double_arr) vector<double>;
  %template(Qnum_arr) vector<uni10::Qnum>;
  %template(Bond_arr) vector<uni10::Bond>;
  %template(Qnum2int) map<uni10::Qnum, int>;
  %template(Matrix_arr) vector<uni10::Matrix>;
  %template(Qnum2Matrix) std::map<uni10::Qnum, uni10::Matrix>;
  /*%template(Swap_arr)  std::vector<uni10::_Swap>;*/
}
%feature("autodoc");

%exception {
    try {
      $action
    } catch (const std::exception &e) {
      std::string s("\nException raised by pyUni10: "), s2(e.what());
        s = s + s2;
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "unknown exception");
    }
}

%inline{
  uni10::Qnum QnumF(uni10::parityFType _prtF, int _U1=0, uni10::parityType _prt=uni10::PRT_EVEN){
    return uni10::Qnum(_prtF, _U1, _prt);
  }
};

namespace uni10{
/* Qnum */
enum parityType{
  PRT_EVEN = 0,
  PRT_ODD = 1
};
enum parityFType{
  PRTF_EVEN = 0,
  PRTF_ODD = 1
};
typedef struct{
        int b1;
        int b2;
}_Swap;

class Qnum {
  public:
      Qnum(int _U1 = 0, parityType _prt = PRT_EVEN);
      /*Qnum(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);*/
      Qnum(const Qnum& _q);
      ~Qnum(){};
      int U1()const;
      parityType prt()const;
      parityFType prtF()const;
      void assign(int _U1 = 0, parityType _prt = PRT_EVEN);
      /*void assign(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);*/
      static bool isFermionic(){return Fermionic;}
      /*
         friend bool operator< (const Qnum& q1, const Qnum& q2);
         friend bool operator<= (const Qnum& q1, const Qnum& q2);
         friend bool operator== (const Qnum& q1, const Qnum& q2);
         friend Qnum operator* (const Qnum& q1, const Qnum& q2);
         friend Qnum operator- (const Qnum& q1);
         friend std::ostream& operator<< (std::ostream& os, const Qnum& q);
       */
      %extend {
          /*
          bool __eq__(const Qnum& q2){
              return (*self) == q2;
          }*/
          int __cmp__(const Qnum& q2){
            if((*self) < q2)
              return -1;
            else if((*self) == q2)
              return 0;
            else
              return 1;
          }
          long __hash__(){
            return (*self).hash();
          }
          Qnum __mul__(const Qnum& q2){
              return (*self) * q2;
          }
          Qnum __neg__(){
              return -(*self);
          }
          Qnum __copy__(){
              return (*self);
          }
          const std::string __repr__() {
              std::ostringstream oss(std::ostringstream::out);
              oss << (*self);
              return oss.str();
          }
          Qnum& assignF(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN){
              self->assign(_prtF, _U1, _prt);
              return *self;
          }
      }
      /* Make Qnum Class immutable */
      static const int U1_UPB = 1000;//Upper bound of U1
      static const int U1_LOB = -1000;//Lower bound of U1
};
/* End of Qnum */

/* Bond */
enum bondType{
  BD_IN = 1,
  BD_OUT = -1
};
class Bond {
  public:
      /*Bond(){};*/
      Bond(bondType _type, int dim);
      Bond(bondType, const std::vector<Qnum>& qnums);
      Bond(const Bond& _b);
      void assign(bondType, int dim);
      void assign(bondType, const std::vector<Qnum>& qnums);
      bondType type()const;
      int dim()const;
      std::map<Qnum, int> degeneracy()const;
      std::vector<Qnum> Qlist()const;
      void change(bondType tp);
      Bond& combine(const Bond bd);
      /*
      friend bool operator== (const Bond& b1, const Bond& b2);
      friend Bond combine(bondType tp, const std::vector<Bond>& bds);
      friend Bond combine(const std::vector<Bond>& bds);
      friend std::ostream& operator<< (std::ostream& os, const Bond& b);
      friend class UniTensor;
      friend class Node;
      */
      %extend {
          bool __eq__(const Bond& b2){
              return (*self) == b2;
          }
          Bond __copy__(){
              return (*self);
          }
          const std::string  __repr__() {
              std::ostringstream oss(std::ostringstream::out);
              oss << (*self);
              return oss.str();
          }
      }
      ~Bond();
};
extern Bond combine(bondType tp, const std::vector<Bond>& bds);
extern Bond combine(const std::vector<Bond>& bds);
/* End of Bond */

/* Matrix */
%apply int *OUTPUT { int *lanczos_iter};
class Matrix {
    public:
        Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
        Matrix(size_t _Rnum, size_t _Cnum, double* _elem, bool _diag=false, bool src_ongpu=false);
        Matrix(size_t _Rnum, size_t _Cnum, std::vector<double> _elem, bool _diag=false, bool src_ongpu=false);
		    /*Matrix(const Matrix& _m);*/
        Matrix();
        ~Matrix();
        /*Matrix& operator=(const Matrix& _m);*/
        int row()const;
        int col()const;
        bool isDiag()const{return diag;};
        bool isOngpu()const{return ongpu;};
        size_t elemNum()const;
		    /*double& operator[](size_t idx);*/
		    /*double& at(size_t i, size_t j);*/
        double* getElem()const;
        double* getHostElem();
		    void setElem(const double* elem, bool _ongpu = false);
    		void setElem(const std::vector<double>& elem, bool _ongpu = false);
        Matrix& resize(size_t row, size_t col);
        void save(const std::string& fname);
        void load(const std::string& fname);
        void identity();
        void set_zero();
        void randomize();
        void orthoRand();
        Matrix& transpose();
        std::vector<Matrix> eigh()const;
        std::vector<Matrix> svd()const;
        /*size_t lanczosEigh(double E0, Matrix& psi, size_t max_iter=200, double err_tol = 5E-15);*/
        double trace();
        double norm();
        double sum();
        Matrix& operator*= (double a);
		    Matrix& operator*= (const Matrix& Mb);
        Matrix& operator+= (const Matrix& Mb);
        /*
		    friend Matrix takeExp(double a, const Matrix& mat);
        friend Matrix operator* (const Matrix& Ma, const Matrix& Mb);
        friend Matrix operator*(const Matrix& Ma, double a);
        friend Matrix operator*(double a, const Matrix& Ma){return Ma * a;};
        friend Matrix operator+(const Matrix& Ma, const Matrix& Mb);
        friend bool operator== (const Matrix& m1, const Matrix& m2);
        friend std::ostream& operator<< (std::ostream& os, const Matrix& b);
        */
        %extend {
          bool __eq__(const Matrix& m2){
            return (*self) == m2;
          }
          Matrix __copy__(){
            return (*self);
          }
          const std::string __repr__() {
            std::ostringstream oss(std::ostringstream::out);
            oss << (*self);
            return oss.str();
          }
          Matrix __mul__(const Matrix& Ma){
            return (*self) * Ma;
          }
          Matrix __add__(const Matrix& Ma){
            return (*self) + Ma;
          }
          Matrix __mul__(double a){
            return a * (*self);
          }
          Matrix __rmul__(double a){
            return a * (*self);
          }
          double __getitem__(PyObject *parm) {
            if (PyTuple_Check(parm)){
              long r,c;
              r=PyInt_AsLong(PyTuple_GetItem(parm,0));
              c=PyInt_AsLong(PyTuple_GetItem(parm,1));
              if( (*self).isDiag() && r!=c) {
                return 0.0;
              }
              else {
                return (*self).at(r,c);
              }
            } else if (PyInt_Check(parm))
              return (*self)[PyInt_AsLong(parm)];
          }
          void __setitem__(PyObject *parm, double val){
            if (PyTuple_Check(parm)){
              long r,c;
              r=PyInt_AsLong(PyTuple_GetItem(parm,0));
              c=PyInt_AsLong(PyTuple_GetItem(parm,1));
              if((*self).isDiag()) {
                if (r==c) (*self)[r]=val;
              }
              else {
                (*self)[r*(*self).col()+c]=val;
              }
            }
            else
              if (PyInt_Check(parm)) (*self)[PyInt_AsLong(parm)]=val;
          }
          double lanczosEigh(Matrix& psi, int *lanczos_iter, size_t max_iter=200, double err_tol = 5E-15){
            double E0;
            *lanczos_iter = (*self).lanczosEigh(E0, psi, max_iter, err_tol);
            return E0;
          }
        }
};
Matrix takeExp(double a, const Matrix& mat);
Matrix otimes(const Matrix& Ta, const Matrix& Tb);
%clear int *lanczos_iter;

/* End of Matrix */


/*class UniTensor;*/
class UniTensor{
  public:
    UniTensor();
    UniTensor(double val);
    UniTensor(const std::string& fname);
    UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "");
    UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
    UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
    /*UniTensor(const UniTensor& UniT);*/
    ~UniTensor();
    /*UniTensor& operator=(const UniTensor& UniT);*/
    UniTensor& assign(const std::vector<Bond>& _bond);
    void setLabel(const std::vector<int>& newLabels);
    void setLabel(int* newLabels);
    std::vector<int> label()const;
    int label(size_t idx)const;
    size_t bondNum()const;
    int inBondNum()const;
    std::vector<Bond> bond()const;
    Bond bond(size_t idx)const;
    size_t elemNum()const;
    Matrix getRawElem()const;
    void setRawElem(const std::vector<double>& rawElem);
    void setRawElem(const double* rawElem);
    double at(const std::vector<int>& idxs)const;
    /*double at(const std::vector<size_t>& idxs)const;*/
    size_t blockNum()const;
    std::vector<Qnum> blockQnum()const;
    Qnum blockQnum(size_t idx)const;
    std::map<Qnum, Matrix> getBlocks()const;
		Matrix getBlock(bool diag = false)const;
    Matrix getBlock(const Qnum& qnum, bool diag = false)const;
		void putBlock(const Matrix& mat);
    void putBlock(const Qnum& qnum, const Matrix& mat);
    double* getElem();
    void setElem(const double* elem, bool _ongpu = false);
    void setElem(const std::vector<double>& elem, bool _ongpu = false);
    /*double operator[](size_t idx);*/
    void set_zero();
    void set_zero(const Qnum& qnum);
    void identity();
    void identity(const Qnum& qnum);
    void randomize();
    void orthoRand();
    void orthoRand(const Qnum& qnum);
    std::string getName();
    void setName(const std::string& _name);
    void save(const std::string& fname);
    UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
    /*UniTensor& permute(int* newLabels, int inBondNum);*/
    UniTensor& permute(int inBondNum);
    UniTensor& transpose();
    UniTensor& combineBond(const std::vector<int>& combined_labels);
    UniTensor& partialTrace(int la, int lb);
    double trace()const;
    std::vector<_Swap> exSwap(const UniTensor& Tb)const;
    void addGate(const std::vector<_Swap>& swaps);
    UniTensor& operator*= (double a);
    UniTensor& operator*= (const UniTensor& Tb);
    UniTensor& operator+= (const UniTensor& Tb);
    bool similar(const UniTensor& Tb)const;
    bool elemCmp(const UniTensor& UniT)const;
    void printRawElem()const;
    static void profile();

    /*
		friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
		friend UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
    friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
    friend UniTensor operator* (const UniTensor& Ta, double a);
    friend UniTensor operator* (double a, const UniTensor& Ta){return Ta * a;};
    friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
    friend std::ostream& operator<< (std::ostream& os, const UniTensor& UniT);
    */

    %extend {
      UniTensor __copy__(){
        return (*self);
      }
      const std::string __repr__() {
        std::ostringstream oss(std::ostringstream::out);
        oss << (*self);
        return oss.str();
      }
      UniTensor __mul__(const UniTensor& Ta){
        return (*self) * Ta;
      }
      UniTensor __add__(const UniTensor& Ta){
        return (*self) + Ta;
      }
      UniTensor __mul__(double a){
        return a * (*self);
      }
      UniTensor __rmul__(double a){
        return a * (*self);
      }
      double __getitem__(PyObject *parm) {
        return (*self)[PyInt_AsLong(parm)];
      }
    }
};
UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast=false);
UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
/* End of UniTensor */

/* class Network */
class Network {
  public:
    Network(const std::string& fname, const std::vector<UniTensor*>& tens);
    Network(const std::string& fname);
    ~Network();
    void putTensor(int idx, const UniTensor& UniT, bool force=true);
    void putTensor(const std::string& name, const UniTensor& UniT, bool force=true);
    void putTensorT(const std::string& nameT, const UniTensor& UniT, bool force=true);
    UniTensor launch(const std::string& name="");
    void profile();
    /*friend std::ostream& operator<< (std::ostream& os, Network& nd);*/
    %extend {
      Network __copy__(){
        return (*self);
      }
      const std::string __repr__() {
        std::ostringstream oss(std::ostringstream::out);
        oss << (*self);
        return oss.str();
      }
    }
};
/* End of Network */


};


%module pyUni10
%{
  /* Put header files here or function declarations like below */
  #include <sstream>
  #include <uni10/datatype/Qnum.h>
  #include <uni10/data-structure/Bond.h>
  #include <uni10/data-structure/Block.h>
  #include <uni10/tensor-network/Matrix.h>
  #include <uni10/tensor-network/UniTensor.h>
  #include <uni10/tensor-network/Network.h>
%}

%begin %{
#ifdef _MSC_VER
#define SWIG_PYTHON_INTERPRETER_NO_DEBUG
#endif
%}

%inline%{
  typedef double Real;
  typedef std::complex<double> Complex;
%}

%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "std_complex.i"
%include "exception.i"
/*%include "typemaps.i"*/
namespace std{
  %template(int_arr) vector<int>;
  %template(double_arr) vector<double>;
  %template(Qnum_arr) vector<uni10::Qnum>;
  %template(Bond_arr) vector<uni10::Bond>;
  %template(Qnum2int) map<uni10::Qnum, int>;
  %template(Matrix_arr) vector<uni10::Matrix>;
  %template(Qnum2Block) std::map<uni10::Qnum, uni10::Block>;
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

/* Block */

%apply int *OUTPUT { int *lanczos_iter};

enum rflag{
  RNULL = 0,
  RTYPE = 1
};

enum cflag{
  CNULL = 0,
  CTYPE = 2
};

class Matrix;
class Block{
  public:
    /******** develop ********/
    /* getElem();            */
    /*************************/
    Block();
    Block(size_t _Rnum, size_t _Cnum, bool _diag = false);
    Block(const Block& _b);
    virtual ~Block();
    size_t row()const;
    size_t col()const;
    bool isDiag()const;
    bool isOngpu()const;
    size_t elemNum()const;
    int typeID()const;
    void save(const std::string& fname)const;
    std::vector<Matrix> qr()const;
    std::vector<Matrix> rq()const;
    std::vector<Matrix> ql()const;
    std::vector<Matrix> lq()const;
    std::vector<Matrix> svd()const;
    std::vector<Matrix> eig()const;
    std::vector<Matrix> eigh()const;
    Matrix inverse()const;
    Real norm()const;
    Matrix getDiag()const;
    Real trace()const;
    Real sum()const;
    Real at(size_t i, size_t j)const;
    bool CelemIsNULL()const;
    bool RelemIsNULL()const;
    /* deafult is real, so rflag is ignored for conflict with cflag */
    Block(cflag _tp, size_t _Rnum, size_t _Cnum, bool _diag = false);
    void save(cflag _tp, const std::string& fname)const;
    std::vector<Matrix> qr(cflag _tp)const;
    std::vector<Matrix> rq(cflag _tp)const;
    std::vector<Matrix> ql(cflag _tp)const;
    std::vector<Matrix> lq(cflag _tp)const;
    std::vector<Matrix> svd(cflag _tp)const;
    std::vector<Matrix> eig(cflag _tp)const;
    std::vector<Matrix> eigh(cflag _tp)const;
    Matrix inverse(cflag _tp)const;
    Real norm(cflag _tp)const;
    Matrix getDiag(cflag _tp)const;
    Complex trace(cflag _tp)const;
    Complex sum(cflag _tp)const;
    /*Complex operator()(size_t idx)const;*/
    Complex at(cflag _tp, size_t i, size_t j)const;
    Complex* getElem(cflag _tp)const;
    /*
       friend Matrix operator*(const Block& Ma, const Block& Mb); //R*R C*C R*C C*R
       friend Matrix operator*(Real a, const Block& Ma);
       friend Matrix operator*(const Block& Ma, Real a);
       friend Matrix operator*(const Complex& a, const Block& Ma);
       friend Matrix operator*(const Block& Ma, const Complex& a);
       friend Matrix operator+(const Block& Ma, const Block& Mb);
       friend bool operator==(const Block& m1, const Block& m2);
       friend bool operator!=(const Block& m1, const Block& m2){return !(m1 == m2);};
       Block(rflag _tp, size_t _Rnum, size_t _Cnum, bool _diag = false);
       void save(rflag _tp, const std::string& fname)const;
       std::vector<Matrix> qr(rflag tp)const;
       std::vector<Matrix> rq(rflag tp)const;
       std::vector<Matrix> ql(rflag tp)const;
       std::vector<Matrix> lq(rflag tp)const;
       std::vector<Matrix> svd(rflag tp)const;
       std::vector<Matrix> eig(rflag tp)const;
       std::vector<Matrix> eigh(rflag tp)const;
       Matrix inverse(rflag tp)const;
       Real norm(rflag tp)const;
       Matrix getDiag(rflag tp)const;
       Real trace(rflag tp)const;
       Real sum(rflag tp)const;
       Real operator[](size_t idx)const;
       Real at(rflag tp, size_t i, size_t j)const;
       Real* getElem(rflag tp = RTYPE)const;
     */
    %extend {
      bool __eq__(const Block& b2){
        return (*self) == b2;
      }
      Block __copy__(){
        return (*self);
      }
      const std::string __repr__() {
        std::ostringstream oss(std::ostringstream::out);
        oss << (*self);
        return oss.str();
      }
      Matrix __mul__(const Block& Ma){
        return (*self) * Ma;
      }
      Matrix __add__(const Block& Ma){
        return (*self) + Ma;
      }
      Matrix __mul__(double a){
        return a * (*self);
      }
      Matrix __rmul__(double a){
        return a * (*self);
      }
      std::complex<double> __getitem__(PyObject *parm) {
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
        return 0;
      }
      /*
      double lanczosEigh(Matrix& psi, int *lanczos_iter, size_t max_iter=200, double err_tol = 5E-15){
        double E0;
        *lanczos_iter = (*self).lanczosEigh(E0, psi, max_iter, err_tol);
        return E0;
      }
      */
    }
};

/* End of Block */


/* Matrix */
class Matrix: public Block {
  public:
    double absMax(bool _ongpu=false);
    Matrix& maxNorm();
    Matrix& absMaxNorm();
    double* getHostElem();
    std::complex<double>* getHostElem(cflag _tp);
    Matrix();
    Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    Matrix(std::string tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    Matrix(const Matrix& _m);
    Matrix(const Block& _b);
    Matrix(const std::string& fname);
    ~Matrix();
    void identity();
    void set_zero();
    void randomize();
    void orthoRand();
    Matrix& normalize();
    Matrix& transpose();
    Matrix& cTranspose();
    Matrix& conj();
    Matrix& resize(size_t row, size_t col);
    double max(bool _ongpu=false);
    void assign(size_t _Rnum, size_t _Cnum);
    void load(const std::string& fname);
    bool toGPU();

    Matrix(size_t _Rnum, size_t _Cnum, const double* _elem, bool _diag=false, bool _ongpu=false, bool src_ongpu=false);
    Matrix(size_t _Rnum, size_t _Cnum, const std::vector<double>& _elem, bool _diag=false, bool _ongpu = false, bool src_ongpu=false);
    void setElem(const double* elem, bool src_ongpu = false);
    void setElem(const std::vector<double>& elem, bool src_ongpu = false);
    Matrix(size_t _Rnum, size_t _Cnum, const std::complex<double>* _elem, bool _diag=false, bool _ongpu=false, bool src_ongpu=false);
    Matrix(size_t _Rnum, size_t _Cnum, const std::vector< std::complex<double> >& _elem, bool _diag=false, bool _ongpu=false, bool src_ongpu=false);
    Matrix(cflag _tp, const std::string& fname);
    Matrix(cflag _tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    void setElem(const std::complex<double>* elem, bool src_ongpu = false);
    void setElem(const std::vector< std::complex<double> >& elem, bool src_ongpu = false);
    void identity(cflag _tp);
    void set_zero(cflag _tp);
    void randomize(cflag _tp);
    void orthoRand(cflag _tp);
    Matrix& normalize(cflag tp);
    Matrix& transpose(cflag _tp);
    Matrix& cTranspose(cflag _tp);
    Matrix& conj(cflag _tp);
    Matrix& resize(cflag _tp, size_t row, size_t col);
    std::complex<double>& operator()(size_t idx); //&
    std::complex<double>& at(cflag _tp, size_t i); //&
    void assign(cflag _tp, size_t _Rnum, size_t _Cnum);
    bool toGPU(cflag _tp);

    /*
       double absMax(rflag tp, bool _ongpu=false);
       Matrix& maxNorm(rflag tp);
       Matrix& absMaxNorm(rflag tp);
       double* getHostElem(rflag _tp);
       Matrix& operator=(const Matrix& _m);
       Matrix& operator=(const Block& _m);
       Matrix& operator*= (double a);
       Matrix& operator*= (std::complex<double> a);
       Matrix& operator*= (const Block& Mb);
       Matrix& operator+= (const Block& Mb);
       Matrix(rflag _tp, const std::string& fname);
       Matrix(rflag _tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
       void identity(rflag _tp);
       void set_zero(rflag _tp);
       void randomize(rflag _tp);
       void orthoRand(rflag _tp);
       Matrix& normalize(rflag tp);
       Matrix& transpose(rflag _tp);
       Matrix& cTranspose(rflag _tp);
       Matrix& conj(rflag _tp);
       Matrix& resize(rflag _tp, size_t row, size_t col);
       double max(rflag _tp, bool _ongpu=false);
       double& at(size_t i, size_t j); //&
       double& at(rflag _tp, size_t i); //&
       double& operator[](size_t idx); //&
       void assign(rflag _tp, size_t _Rnum, size_t _Cnum);
       bool toGPU(rflag _tp);
     */

        %extend {
          bool __eq__(const Matrix& b2){
            return (*self) == b2;
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
          std::complex<double> __getitem__(PyObject *parm) {
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
            return 0;
          }
          /*
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
          */
        }
};
Matrix takeExp(double a, const Block& mat);
Matrix otimes(const Block& Ta, const Block& Tb);
%clear int *lanczos_iter;
/* End of Matrix */


/*class UniTensor;*/
class UniTensor{
  public:

        double max() const;

        double absMax() const;

        UniTensor& maxNorm();

        UniTensor& absMaxNorm();

        double norm() const;
        double norm(cflag tp) const;
        /*
        UniTensor& normalize();
        double max(rflag tp) const;
        double absMax(rflag tp) const;
        UniTensor& maxNorm(rflag tp);
        UniTensor& absMaxNorm(rflag tp);
        double norm(rflag tp) const;
        UniTensor& normalize(rflag tp);
        UniTensor& normalize(cflag tp);
        */
        void printGraphy()const;

        std::vector<UniTensor> hosvd(int* group_labels, int* groups, size_t groupsSize, std::vector<Matrix>& Ls)const ;
        std::vector<UniTensor> hosvd(int* group_labels, int* groups, size_t groupsSize, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const ;
        std::vector<UniTensor> hosvd(std::vector<int>& group_labels, std::vector<int>& groups, std::vector<Matrix>& Ls)const ;
        std::vector<UniTensor> hosvd(std::vector<int>& group_labels, std::vector<int>& groups, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const ;

        std::vector<UniTensor> hosvd(cflag tp, int* group_labels, int* groups, size_t groupsSize, std::vector<Matrix>& Ls)const ;
        std::vector<UniTensor> hosvd(cflag tp, int* group_labels, int* groups, size_t groupsSize, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const ;
        std::vector<UniTensor> hosvd(cflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<Matrix>& Ls)const ;
        std::vector<UniTensor> hosvd(cflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const ;

        std::vector<UniTensor> hosvd(size_t modeNum, size_t fixedNum = 0)const;
        std::vector<UniTensor> hosvd(size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const;
        std::vector<UniTensor> hosvd(size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const;
        std::vector<UniTensor> hosvd(cflag tp, size_t modeNum, size_t fixedNum = 0)const;
        std::vector<UniTensor> hosvd(cflag tp, size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const;
        std::vector<UniTensor> hosvd(cflag tp, size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const;
        UniTensor& cTranspose();
        UniTensor& cTranspose(cflag tp);
        UniTensor();
        UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "");
        UniTensor(const std::string tp, const std::vector<Bond>& _bonds, const std::string& _name = "");
        UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
        UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
        UniTensor(const UniTensor& UniT);
        UniTensor(const std::string& fname);
        UniTensor(const std::string& fname, const bool hdf5);
        UniTensor(const Block& UniT);
        ~UniTensor();
        void setRawElem(const Block& blk);
        void putBlock(const Block& mat);
        void putBlock(const Qnum& qnum, const Block& mat);
        int typeID()const;
        void setLabel(const int newLabel, const size_t idx);
        void setLabel(const std::vector<int>& newLabels);
        void setLabel(int* newLabels);
        std::vector<int> label()const;
        int label(size_t idx)const;
        std::string getName() const;
        void setName(const std::string& _name);
        size_t bondNum()const;
        size_t inBondNum()const;
        std::vector<Bond> bond()const;
        Bond bond(size_t idx)const;
        size_t elemNum()const;
        size_t blockNum()const;
        std::vector<Qnum> blockQnum()const;
        Qnum blockQnum(size_t idx)const;
        const std::map<Qnum, Block>& const_getBlocks()const;
        const Block& const_getBlock()const;
        const Block& const_getBlock(const Qnum& qnum)const;
        std::map<Qnum, Matrix> getBlocks()const;
        Matrix getBlock(bool diag = false)const;
        Matrix getBlock(const Qnum& qnum, bool diag = false)const;
        void set_zero();
        void set_zero(const Qnum& qnum);
        void identity();
        void identity(const Qnum& qnum);
        void randomize();
        void orthoRand();
        void orthoRand(const Qnum& qnum);
        void save(const std::string& fname) const;
        void h5save(const std::string& fname);
        UniTensor& transpose();
        UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
        UniTensor& permute(int* newLabels, int inBondNum);
        UniTensor& permute(int inBondNum);
        UniTensor& combineBond(const std::vector<int>& combined_labels);
        static std::string profile(bool );
        std::vector<_Swap> exSwap(const UniTensor& Tb)const;
        void addGate(const std::vector<_Swap>& swaps);
        std::complex<double> trace()const;
        UniTensor& partialTrace(int la, int lb);
        Matrix getRawElem()const;
        UniTensor& assign(const std::vector<Bond>& _bond);
        bool CelemIsNULL();
        bool RelemIsNULL();
        bool similar(const UniTensor& Tb)const;
        bool elemCmp(const UniTensor& UniT)const;
        void clear();

        UniTensor(double val);
        void setRawElem(const std::vector<double>& rawElem);
        void setRawElem(const double* rawElem);
        void setElem(const double* elem, bool _ongpu = false);
        void setElem(const std::vector<double>& elem, bool _ongpu = false);

        UniTensor(std::complex<double> val);
        UniTensor(cflag tp, const std::vector<Bond>& _bonds, const std::string& _name = "");
        UniTensor(cflag tp, const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
        UniTensor(cflag tp, const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
        void setRawElem(const std::vector< std::complex<double> >& rawElem);
        void setRawElem(const std::complex<double>* rawElem);
        void setRawElem(cflag tp, const Block& blk);
        void putBlock(cflag tp, const Block& mat);
        void putBlock(cflag tp, const Qnum& qnum, const Block& mat);
        void setElem(const std::complex<double>* c_elem, bool _ongpu = false);
        void setElem(const std::vector< std::complex<double> >& c_elem, bool _ongpu = false);
        std::map<Qnum, Matrix> getBlocks(cflag tp)const;
        Matrix getBlock(cflag tp, bool diag = false)const;
        Matrix getBlock(cflag tp, const Qnum& qnum, bool diag = false)const;
        void set_zero(cflag tp);
        void set_zero(cflag tp, const Qnum& qnum);
        void identity(cflag tp);
        void identity(cflag tp, const Qnum& qnum);
        void randomize(cflag tp);
        void orthoRand(cflag tp);
        void orthoRand(cflag tp, const Qnum& qnum);
        UniTensor& transpose(cflag tp);
        UniTensor& permute(cflag tp, const std::vector<int>& newLabels, int inBondNum);
        UniTensor& permute(cflag tp, int* newLabels, int inBondNum);
        UniTensor& permute(cflag tp, int inBondNum);
        UniTensor& combineBond(cflag tp, const std::vector<int>& combined_labels);
        void addGate(cflag tp, const std::vector<_Swap>& swaps);
        std::complex<double>* getElem(cflag tp);
        Matrix getRawElem(cflag tp)const;
        std::complex<double> trace(cflag tp)const;
        UniTensor& partialTrace(cflag tp, int la, int lb);
        UniTensor& assign(cflag tp, const std::vector<Bond>& _bond);

        double at(size_t idx)const;
        double at(const std::vector<int>& idxs)const;
        double at(const std::vector<size_t>& idxs)const;
        std::complex<double> at(cflag tp, size_t idx)const;
        std::complex<double> at(cflag tp, const std::vector<int>& idxs)const;
        std::complex<double> at(cflag tp, const std::vector<size_t>& idxs)const;
        double* getElem();
        
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
      std::complex<double> __getitem__(PyObject *parm) {
        return (*self)[PyInt_AsLong(parm)];
      }
      static const std::string profile(){
        return uni10::UniTensor::profile(false);
      }
      const std::string printRawElem(){
        return (*self).printRawElem(false);
      }
    }
    /*
       std::vector<UniTensor> hosvd(rflag tp, int* group_labels, int* groups, size_t groupsSize, std::vector<Matrix>& Ls)const ;
       std::vector<UniTensor> hosvd(rflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<Matrix>& Ls)const ;
       std::vector<UniTensor> hosvd(rflag tp, int* group_labels, int* groups, size_t groupsSize, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const ;
       std::vector<UniTensor> hosvd(rflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const ;
       std::vector<UniTensor> hosvd(rflag tp, size_t modeNum, size_t fixedNum = 0)const;
       std::vector<UniTensor> hosvd(rflag tp, size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const;
       std::vector<UniTensor> hosvd(rflag tp, size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const;
       UniTensor& operator*= (double a);
       UniTensor& operator*= (std::complex<double> a);
       UniTensor& operator*= (const UniTensor& Tb);
       friend UniTensor operator*(const UniTensor& Ta, const std::complex<double>& a);
       friend UniTensor operator*(const std::complex<double>& a, const UniTensor& Ta);
       friend UniTensor operator*(const UniTensor& Ta, double a);
       friend UniTensor operator*(double a, const UniTensor& Ta);
       friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
       UniTensor& operator=(const UniTensor& UniT);
       UniTensor& operator+= (const UniTensor& Tb);
       friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
       friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
       friend UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
       std::string printRawElem(bool print=true)const;
       UniTensor(rflag tp, const std::vector<Bond>& _bonds, const std::string& _name = "");
       UniTensor(rflag tp, const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
       UniTensor(rflag tp, const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
       void setRawElem(rflag tp, const Block& blk);
       void putBlock(rflag tp, const Block& mat);
       void putBlock(rflag tp, const Qnum& qnum, const Block& mat);
       void set_zero(rflag tp);
       void set_zero(rflag tp, const Qnum& qnum);
       void identity(rflag tp);
       void identity(rflag tp, const Qnum& qnum);
       void randomize(rflag tp);
       void orthoRand(rflag tp);
       void orthoRand(rflag tp, const Qnum& qnum);
       std::map<Qnum, Matrix> getBlocks(rflag tp)const;
       Matrix getBlock(rflag tp, bool diag = false)const;
       Matrix getBlock(rflag tp, const Qnum& qnum, bool diag = false)const;
       UniTensor& transpose(rflag tp);
       UniTensor& permute(rflag tp, const std::vector<int>& newLabels, int inBondNum);
       UniTensor& permute(rflag tp, int* newLabels, int inBondNum);
       UniTensor& permute(rflag tp, int inBondNum);
       friend UniTensor contract(rflag tp, UniTensor& Ta, UniTensor& Tb, bool fast);
       friend UniTensor otimes(rflag tp, const UniTensor& Ta, const UniTensor& Tb);
       UniTensor& combineBond(rflag tp, const std::vector<int>& combined_labels);
       void addGate(rflag tp, const std::vector<_Swap>& swaps);
       double* getElem(rflag tp);
       Matrix getRawElem(rflag tp)const;
       UniTensor& partialTrace(rflag tp, int la, int lb);
       UniTensor& assign(rflag tp, const std::vector<Bond>& _bond);
       friend UniTensor contract(cflag tp, UniTensor& Ta, UniTensor& Tb, bool fast);
       friend UniTensor otimes(cflag tp, const UniTensor& Ta, const UniTensor& Tb);
       double at(rflag tp, size_t idx)const;
       double at(rflag tp, const std::vector<int>& idxs)const;
       double at(rflag tp, const std::vector<size_t>& idxs)const;
       double operator[](size_t idx) const;
       std::complex<double> operator()(size_t idx) const;

     */
};
UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast=false);
UniTensor contract(cflag tp, UniTensor& Ta, UniTensor& Tb, bool fast = false);
UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
UniTensor otimes(cflag tp, const UniTensor& Ta, const UniTensor& Tb);
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
    /*void profile();*/
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
      const std::string profile() {
        return (*self).profile(false);
      }
    }
};
/* End of Network */


};


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


%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "exception.i"
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
          bool __eq__(const Qnum& q2){
              return (*self) == q2;
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
      // Make Qnum Class immutable
      /* %pythoncode{
      def __setattr__(self, *args):
         raise TypeError("can't modify immutable instance")
      __delattr__ = __setattr__
      } */
      static const int U1_UPB = 100;//Upper bound of U1
      static const int U1_LOB = -100;//Lower bound of U1
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
      /*
      friend class UniTensor;
      friend class Node;
      friend std::ostream& operator<< (std::ostream& os, const Bond& b);
      friend bool operator== (const Bond& b1, const Bond& b2);
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
      void change(bondType tp);
      Bond& combine(const Bond bd);
      /*
      friend Bond combine(bondType tp, const std::vector<Bond>& bds);
      friend Bond combine(const std::vector<Bond>& bds);
      */
      std::map<Qnum, int> degeneracy()const;
      std::vector<Qnum> Qlist()const;
      ~Bond();
};
extern Bond combine(bondType tp, const std::vector<Bond>& bds);
extern Bond combine(const std::vector<Bond>& bds);
/* End of Bond */

/* Matrix */
class Matrix {
    public:
        Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
        Matrix(size_t _Rnum, size_t _Cnum, double* _elem, bool _diag=false, bool src_ongpu=false);
        Matrix(size_t _Rnum, size_t _Cnum, std::vector<double> _elem, bool _diag=false, bool src_ongpu=false);
        Matrix();
        ~Matrix();
        int row()const;
        int col()const;
        bool isDiag()const{return diag;};
        bool isOngpu()const{return ongpu;};
        size_t elemNum()const;
        /*Matrix& operator=(const Matrix& _m);*/
        /*Matrix& operator*= (const Matrix& Mb);*/
        std::vector<Matrix> eigh();
        std::vector<Matrix> svd();
        size_t lanczosEigh(double& E0, Matrix& psi, size_t max_iter=200, double err_tol = 5E-15);
        void setElem(double* elem, bool _ongpu = false);
        void setElem(std::vector<double> elem, bool _ongpu = false);
        double* getElem()const;
        double* getHostElem();
        void randomize();
        void orthoRand();
        void identity();
        void set_zero();
        Matrix& transpose();
        Matrix& resize(size_t row, size_t col);
        double trace();
        double norm();
        double sum();
        void save(const std::string& fname);
        void load(const std::string& fname);
        /*
        friend std::ostream& operator<< (std::ostream& os, const Matrix& b);
        friend Matrix operator* (const Matrix& Ma, const Matrix& Mb);
        friend Matrix operator*(const Matrix& Ma, double a);
        friend Matrix operator*(double a, const Matrix& Ma){return Ma * a;};
        friend bool operator== (const Matrix& m1, const Matrix& m2);
        friend Matrix operator+(const Matrix& Ma, const Matrix& Mb);
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
        }
        Matrix& operator*= (double a);
        Matrix& operator+= (const Matrix& Mb);
        bool toGPU();
};
Matrix takeExp(double a, const Matrix& mat);


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
    UniTensor(const UniTensor& UniT);
    /*UniTensor& operator=(const UniTensor& UniT);*/
    UniTensor& assign(const std::vector<Bond>& _bond);
    ~UniTensor();
    void setLabel(const std::vector<int>& newLabels);
    void setLabel(int* newLabels);
    void setRawElem(std::vector<double> rawElem);
    double at(std::vector<int>idxs)const;
    /*double& operator[](size_t idx);*/
    std::vector<Qnum> blockQnum()const;
    Qnum blockQnum(int idx)const;
    size_t blockNum()const;
    void save(const std::string& fname);
    std::vector<int> label()const;
    int label(int idx)const;
    std::vector<Bond> bond()const;
    Bond bond(int idx)const;
    void setName(const std::string& _name);
    std::string getName();
    size_t elemNum()const;
    size_t bondNum()const;
    int inBondNum()const;
    static void profile();
    UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
    UniTensor& permute(int* newLabels, int inBondNum);
    UniTensor& permute(int inBondNum);
    UniTensor& transpose();
    void randomize();
    /*friend std::ostream& operator<< (std::ostream& os, const UniTensor& UniT);*/

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
    }
    UniTensor& operator*= (const UniTensor& Tb);
    /*friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
      friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
     */
    UniTensor& operator+= (const UniTensor& Tb);
    /*friend UniTensor operator* (const UniTensor& Ta, double a);
      friend UniTensor operator* (double a, const UniTensor& Ta){return Ta * a;};*/
    %rename(__add__) UniTensor::operator+;
    %rename(__mul__) UniTensor::operator*;
    UniTensor& operator*= (double a);
    Matrix getBlock(const Qnum& qnum, bool diag = false)const;
    void putBlock(const Qnum& qnum, const Matrix& mat);
    std::map<Qnum, Matrix> getBlocks()const;
    Matrix getRawElem()const;
    void printRawElem()const;
    /*              friend class Node;
                    friend class Network;*/
    void orthoRand();
    void orthoRand(const Qnum& qnum);
    void identity();
    void identity(const Qnum& qnum);
    void set_zero(const Qnum& qnum);
    void set_zero();
    double* getElem();
    void setElem(double* elem, bool _ongpu = false);
    void setElem(std::vector<double>& elem, bool _ongpu = false);
    std::vector<_Swap> exSwap(const UniTensor& Tb)const;
    bool similar(const UniTensor& Tb)const;
    void addGate(std::vector<_Swap> swaps);
    bool elemCmp(const UniTensor& UniT)const;
    double trace()const;
    UniTensor& combineBond(const std::vector<int>& combined_labels);
    UniTensor& partialTrace(int la, int lb);
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
    /*
    void putTensor(int idx, const UniTensor* UniT, bool force=true);
    void putTensor(const std::string& name, const UniTensor* UniT, bool force=true);
    void putTensorT(const std::string& nameT, const UniTensor* UniT, bool force=true);
    */
    void putTensor(int idx, const UniTensor& UniT, bool force=true);
    void putTensor(const std::string& name, const UniTensor& UniT, bool force=true);
    void putTensorT(const std::string& nameT, const UniTensor& UniT, bool force=true);
    UniTensor launch(const std::string& name="");
    /*friend std::ostream& operator<< (std::ostream& os, Network& nd);*/
    //int rollcall();
    //size_t max_tensor_elemNum();
    //size_t sum_of_memory_usage();
    //size_t memory_requirement();
    void profile();
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
  /*
  private:
    void preprint(std::ostream& os, Node* nd, int layer);	//pre-order print
    std::vector<std::string> names;
    std::map<std::string, int> name2pos;
    std::vector< std::vector<int> > label_arr;
    std::vector< int > Rnums;
    std::vector<Node*> leafs;
    std::vector<UniTensor*> tensors;
    std::vector< std::vector<_Swap> > swaps_arr;
    std::vector<bool> swapflags;
    std::vector<int> conOrder;	//contraction order;
    std::vector<int> order;	//add order
    std::vector<int> brakets;	//add order
    Node* root;
    bool load;	//whether or not the network is ready for contraction, construct=> load=true, destruct=>load=false
    int times;	//construction times
    int tot_elem;	//total memory ussage
    int max_elem;	//maximum
    void construct();
    void destruct();
    void matching(Node* sbj, Node* tar);
    void branch(Node* sbj, Node* tar);
    UniTensor merge(Node* nd);
    void clean(Node* nd);
    void fromfile(const std::string& fname);
    void findConOrd(Node* nd);
    void addSwap();
    void _max_tensor_elemNum(Node* nd, size_t& max_num, Node& max_nd) const;
    size_t _sum_of_tensor_elem(Node* nd) const;
    size_t _elem_usage(Node* nd, size_t& usage, size_t& max_usage)const;
  */
};
/* End of Network */


};


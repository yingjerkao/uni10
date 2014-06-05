%module pyUni10
%{
  /* Put header files here or function declarations like below */
  #include <sstream>
  #include <uni10/datatype/Qnum.h>
  #include <uni10/data-structure/Bond.h>
  #include <uni10/tensor-network/Matrix.h>
  #include <uni10/tensor-network/UniTensor.h>
%}
%include "std_vector.i"
%include "std_map.i"
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
          Qnum cp(){
              return (*self);
          }
          const char* __str__() {
              std::ostringstream out;
              out << *self;
              return out.str().c_str();
          }
          Qnum& assignF(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN){
              self->assign(_prtF, _U1, _prt);
              return *self;
          }
      }
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
          Bond cp(){
              return (*self);
          }
          const char* __str__() {
              std::ostringstream out;
              out << *self;
              return out.str().c_str();
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
        Matrix(int _Rnum, int _Cnum, bool _diag=false);
        Matrix(int _Rnum, int _Cnum, double* _elem, bool _diag=false);
        Matrix(const Matrix& _m);
        ~Matrix();
        int row()const;
        int col()const;
        bool isDiag()const{return diag;};
        size_t elemNum()const;
        /*Matrix& operator=(const Matrix& _m);*/
        /*Matrix& operator*= (const Matrix& Mb);*/
        std::vector<Matrix> diagonalize();
        std::vector<Matrix> svd();
        void addElem(double* elem);
        void randomize();
        void orthoRand();
        void set_zero();
        void transpose();
        double trace();
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
          Matrix cp(){
              return (*self);
          }
          const char* __str__() {
              std::ostringstream out;
              out << *self;
              return out.str().c_str();
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
      }
        Matrix& operator*= (double a);
        Matrix& operator+= (const Matrix& Mb);
        /*double& operator[](size_t idx);*/
        double* elem()const;
        double& at(int i, int j);
};

/* End of Matrix */


/*class UniTensor;*/
class UniTensor{
        public:
                UniTensor(double val = 1.0);
                UniTensor(const std::string& fname);
                UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "");
                UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
                UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
                UniTensor(const UniTensor& UniT);
                /*UniTensor& operator=(const UniTensor& UniT);*/
                UniTensor& assign(const std::vector<Bond>& _bond);
                ~UniTensor();
                void addLabel(const std::vector<int>& newLabels);
                void addLabel(int* newLabels);
                void addRawElem(double* rawElem);
                void elemSet(const Qnum& qnum, double* _elem);
                void elemSet(double* _elem);
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
                static void check();
                UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
                UniTensor& permute(int* newLabels, int inBondNum);
                UniTensor& permute(int inBondNum);
                UniTensor& transpose();
                void randomize();
                /*friend std::ostream& operator<< (std::ostream& os, const UniTensor& UniT);*/

      %extend {
          UniTensor copy__(){
              return (*self);
          }
          UniTensor cp(){
              return (*self);
          }
          const char* __str__() {
              std::ostringstream out;
              out << *self;
              return out.str().c_str();
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
                Matrix rawElem()const;
                void printRawElem()const;
/*                friend class Node;
                friend class Network;*/
                void orthoRand();
                void orthoRand(const Qnum& qnum);
                void eye();
                void eye(const Qnum& qnum);
                void set_zero(const Qnum& qnum);
                void set_zero();
                std::vector<_Swap> exSwap(const UniTensor& Tb)const;
                bool similar(const UniTensor& Tb)const;
                void addGate(std::vector<_Swap> swaps);
                bool elemCmp(const UniTensor& UniT)const;
                double trace()const;
                UniTensor& combineBond(const std::vector<int>& combined_labels);
                UniTensor& partialTrace(int la, int lb);
};
UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast); 
UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
};

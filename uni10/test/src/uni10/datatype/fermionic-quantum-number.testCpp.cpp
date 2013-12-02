/*****************************************************************************
*
* Universal Tensor Network Library
*
* Copyright (C) 2013-2014 
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
* Please direct any enquiry to <development@uni10.org>
*
*****************************************************************************/

#include <iostream>
#include <uni10.hpp>

int main(int agrc, char** argv)
{
  // Usual matter
  uni10::print_copyright(std::cout);

  std::cout << "Test -- fermionic-quantum-number\n"
            << "================================\n"
            << "\n";

  // Test 1
  if (true)
  {
    typedef uni10::datatype::fermionic_quantum_number<int,char> Qnum_t;

    // Constructors
    Qnum_t q1, q2(2), q3(2,1), q4(2,1,1);

    std::cout << "Constructors:\n\n"
              << "q1 : " << q1 << " , q1.U1() : " << q1.U1() << " , q1.prt() : " << q1.prt() << "\n"
              << "q2 : " << q2 << " , q2.U1() : " << q2.U1() << " , q2.prt() : " << q2.prt() << "\n" 
              << "q3 : " << q3 << " , q3.U1() : " << q3.U1() << " , q3.prt() : " << q3.prt() << "\n"
              << "q4 : " << q4 << " , q4.U1() : " << q4.U1() << " , q4.prt() : " << q4.prt() << "\n"
              << "\n";

    // Assignment operator
    Qnum_t q4_copy = q4;

    std::cout << "Assignment operator:\n\n"
              << "q4_copy = q4;\n"
              << "q4_copy : " << q4_copy << "\n"
              << "\n";

    // Set 
    q1.set(1,1);
    q2.set(1,1,1);
    q3.set_U1(3);
    q3.set_prt(0);
    q3.set_prtF(1);

    std::cout << "Set:\n\n"
              << "q1.set(1,1);\n"
              << "q1 : " << q1 << "\n"
              << "q2.set(1,1,1);\n"
              << "q2 : " << q2 << "\n"
              << "q3.set_U1(3);\n"
              << "q3.set_prt(0);\n"
              << "q3.set_prtF(1);\n"
              << "q3 : " << q3 << "\n" 
              << "\n";


    // Comparison operators
/*
    std::cout << "Comparison operators:\n\n"
              << "Qnum_t(1,1) <  Qnum_t(2,1,1) : " << (Qnum_t(1,1) <  Qnum_t(2,1,1)) << "\n"
              << "Qnum_t(1,1) <= Qnum_t(2,1,1) : " << (Qnum_t(1,1) <= Qnum_t(2,1,1)) << "\n"
              << "Qnum_t(1,1) == Qnum_t(2,1,1) : " << (Qnum_t(1,1) == Qnum_t(2,1,1)) << "\n"
              << "\n"; 

    // Arithmetic operator
    std::cout << "Arithmetic operators:\n\n"
              << "-Qnum_t(1)               : " << (-Qnum_t(1)) << "\n"
              << "Qnum_t(1) * Qnum_t(2,1)  : " << (Qnum_t(1) * Qnum_t(2,1)) << "\n"
              << "\n";
*/

  }

  std::cout << "\n"
            << "End Test\n"
            << "\n";

  return 0;
}

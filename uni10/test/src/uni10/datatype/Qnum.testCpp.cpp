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

#include <uni10.hpp>

int main(int agrc, char** argv)
{
  std::cout << uni10::print_copyright() << "\n\n"
            << "test -- Qnum.h\n"
            << "==============\n"
            << "\n\n";

  if (true)
  {
    typedef uni10::datatype::Qnum<int,char> Qnum_t;

    Qnum_t q1, q2(2), q3(2,4);

    std::cout << "q1 = " << q1 << "\n"
              << "q2 = " << q2 << "\n" 
              << "q3 = " << q3 << "\n"
    ;


  }

  return 0;
}

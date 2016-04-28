/****************************************************************************
*  @file uni10_tools.cpp
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
*  @brief Implementation file for useful functions
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <uni10/tools/uni10_tools.h>
#include <string.h>
namespace uni10 {

size_t MEM_USAGE = 0;
size_t GPU_MEM_USAGE = 0;

std::vector<_Swap> recSwap(std::vector<int>& _ord) { //Given the reshape order out to in.
    //int ordF[n];
    int n = _ord.size();
    std::vector<int> ordF(n);
    for(int i = 0; i < n; i++)
        ordF[i] = i;
    return recSwap(_ord, ordF);
}
std::vector<_Swap> recSwap(std::vector<int>& _ord, std::vector<int>& ordF) { //Given the reshape order out to in.
    int n = _ord.size();
    std::vector<int> ord = _ord;
    std::vector<_Swap> swaps;
    _Swap sg;
    int tmp;
    for(int i = 0; i < n - 1; i++)
        for(int j = 0; j < n - i - 1; j++)
            if(ord[j] > ord[j + 1]) {
                sg.b1 = ordF[ord[j + 1]];
                sg.b2 = ordF[ord[j]];
                tmp = ord[j];
                ord[j] = ord[j + 1];
                ord[j + 1] = tmp;
                swaps.push_back(sg);
            }
    return swaps;
}

void propogate_exception(const std::exception& e, const std::string& msg) {
    std::string except_str("\n");
    except_str.append(msg);
    except_str.append(e.what());
    throw std::logic_error(except_str);
}

std::string exception_msg(const std::string& msg) {
    std::string except_str("\n  ");
    except_str.append(msg);
    except_str.append("\n");
    return except_str;
}

};  /* namespace uni10 */

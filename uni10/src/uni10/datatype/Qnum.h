/****************************************************************************
*  @file Qnum.h
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
*  @brief Header file for Qnum (Quantum Number) class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef QNUM_H
#define QNUM_H
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <stdexcept>
#include <sstream>
#include <exception>

namespace uni10 {

//! Parity/Z2 types
enum parityType {
    PRT_EVEN = 0, ///< Parity/Z2 even
    PRT_ODD = 1   ///< Parity/Z2 odd
};
//! Fermion parity types
enum parityFType {
    PRTF_EVEN = 0, ///< Fermion parity even
    PRTF_ODD = 1   ///< Fermion parity odd
};



/// @class Qnum
/// @brief The Qnum class defines the quantum number
///
///  A Qnum is comprised of one or more eigenvalues of following three symmetry operators,
///
/// \f$ U(1)\f$ : \c U(1) eigenvalues; \ref U1 in Qnum.
///
/// \f$ Z_2 \f$ (particle number parity): conservation of the  parity of total
/// particle number or total up/down spins.
///
///  \f$Z_2^F \f$ (fermionic number parity): conservation of fermionic number parity.
///  In a fermionic system, this symmetry is used to keep track the fermionic signs in operations.
class Qnum {
public:

    /// @brief Default constructor
    ///
    /// Creates a Qnum
    /// @param _U1 U1 quantum number, defaults to 0
    /// @param _prt Particle parity/Z2, defaults to \c PRT_EVEN
    Qnum(int _U1=0, parityType _prt=PRT_EVEN);

    /// @brief Creates a fermionic quantum number
    ///
    /// Creates a Qnum with fermionic parity
    /// @param _prtF Fermionic number parity
    /// @param _U1 U1 quantum number, defaults to 0
    /// @param _prt Particle parity/Z2, defaults to \c PRT_EVEN
    Qnum(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);

    /// @brief Copy constructor
    ///
    Qnum(const Qnum& _q);

    /// @brief Default destructor
    ///
    ~Qnum() {};

    /// @brief Access U1
    ///
    /// Returns the value of U1 quantum number, bounded by (Qnum::U1_UPB, Qnum::U1_LOB)
    /// @return U1 quantum number
    int U1()const;

    /// @brief Access paritcle parity/Z2
    ///
    /// Returns the value of parity/Z2 quantum number.
    /// @return Parity/Z2 quantum number
    parityType prt()const;

    /// @brief Access fermionic number parity
    ///
    /// Returns the value of fermionic number parity.
    /// @return Fermionic number parity
    parityFType prtF()const;

    /// @brief Assign Qnum content
    ///
    /// Assigns new content to Qnum, replacing its current content.
    /// @param _U1 U1 quantum number, defaults to 0
    /// @param _prt parity/Z2, defaults to \c PRT_EVEN
    void assign(int _U1 = 0, parityType _prt = PRT_EVEN);

    /// @brief Assign fermionic content
    ///
    /// Assigns new fermionic quantum number to Qnum, replacing its current content.
    /// @param _prtF fermionic number parity
    /// @param _U1 U1 quantum number, defaults to 0
    /// @param _prt parity/Z2 quantum number, defaults to \c PRT_EVEN
    ///
    /// @warning Notice the default values for _U1 and _prt will be assigned if no values are given
    void assign(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);

    /// @brief Test whether the system is fermionic
    ///
    /// Tests whether fermionic parity \c PRTF_ODD exists
    /// @return \c True if the fermionic odd parity exists; \c False otherwise.
    static bool isFermionic() {
        return Fermionic;
    }
    long int hash()const;

    /// @brief Define less than operator
    ///
    /// Defines \c q1 < \c q2
    /// @return \c True if \c q1 < \c q2; \c False otherwise
    friend bool operator< (const Qnum& q1, const Qnum& q2);

    /// @brief Define less than or equal to operator
    ///
    /// Defines \c q1 <= \c q2
    /// @return \c True if \c q1 <= \c q2; \c False otherwise
    friend bool operator<= (const Qnum& q1, const Qnum& q2);

    /// @brief Define equal to operator
    ///
    /// Defines \c q1 == \c q2
    /// @return \c True if \c q1 == \c q2; \c False otherwise
    friend bool operator== (const Qnum& q1, const Qnum& q2);

    /// @brief Define negation operator
    ///
    /// Negation is defined as the operation on the quantum number when an incoming bond is permuted to an
    /// outcoming bond, or vice versa.
    /// For U(1) quantum number \c U1,  negation give \c -U1. No effects on the parity/fermionic parity.
    friend Qnum operator- (const Qnum& q1);

    /// @brief Define multiplication operator
    ///
    /// Defines the fusion rules for quantum numbers.
    ///
    /// For U1 symmetry, <tt> q1*q2 </tt> corresponds to <tt> q1.U1()+q2.%U1()</tt>.
    ///
    /// For parity, <tt> q1*q2  </tt> corresponds to <tt> q1.prt()^q2.%prt()</tt>.
    ///(<tt> q1.prtF()^q2.%prtF() </tt> in the fermionic case)
    /// @return <tt>q1 * q2</tt>
    friend Qnum operator* (const Qnum& q1, const Qnum& q2);
    /// @brief Print out Qnum
    ///
    /// Prints out a quantum number \c q as,
    /// \verbatim (U1 = 1, P = 0, 1) \endverbatim
    /// Where \c q.U1() is 1, \c q.prt() is 0(\c PRT_EVEN) and q.prtF() is 1(\c PRTF_ODD)
    friend std::ostream& operator<< (std::ostream& os, const Qnum& q);
    static const int U1_UPB = 1000; ///<Upper bound of U1 quantum number
    static const int U1_LOB = -1000;///<Lower bound of U1 quantum number
private:
    static bool Fermionic;
    int m_U1;
    parityType m_prt;
    parityFType m_prtF;
};

/// @example egQ1.cpp
/// @example egQ2.cpp

};
#endif /* QNUM_H */

/****************************************************************************
*  @file Bond.h
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
*  @brief Declaration file for Bond Class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef BOND_H
#define BOND_H

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
#include <uni10/datatype.hpp>
#include <stdexcept>
#include <sstream>

namespace uni10 {

//!  Bond types
enum bondType {
    BD_IN = 1, ///<Defines an incoming Bond
    BD_OUT = -1  ///<Defines an outgoing Bond
};

class UniTensor;

/// @brief The Bond class holds the information of a bond.
///
/// A bond is defined from the quantum states specified by quantum numbers defined by the Qnum class.
/// For example, a bond representing the states of a spin-1 particle with local basis states of \f$|-1\rangle,
/// |0\rangle, |1\rangle\f$. In this case, the dimension of the bond is three. Each state of the bond carries
///    a Qnum, which is the eigenvalue of the symmetry operators.
///
///  For example, For a spin system with \f$U(1)\f$ symmetry which conserves the total \f$S_z\f$, we can
/// define a bond with three distinct Qnum with \f$S_z= -1, 0, 1\f$.
///
/// @see \ref bondType, Qnum, UniTensor
class Bond {
public:


    /// @brief Default constuctor
    ///
    Bond() {};

    ///
    /// @brief Create a Bond of type \c tp and dimension \c dim
    /// @param tp  Type of bond
    /// @param dim  Dimension of bond
    Bond(bondType tp, size_t dim);

    /// @brief Create a Bond of type \c tp and dimension \c dim
    /// @param tp  Type of bond
    /// @param qnums Vector of quantum numbers
    Bond(bondType tp, const std::vector<Qnum>& qnums);

    /// @brief Copy constructor
    /// @param bd Reference to a second Bond

    Bond(const Bond& bd);

    /// @brief Destructor
    ///
    ~Bond();

    /*!
     @brief Assign bond content

     Assigns type tp and dimension dim to Bond, replacing the current content.
     @param tp Type of bond,  either \ref BD_IN or \ref BD_OUT
     @param dim Dimension
     */
    void assign(bondType tp, size_t dim);
    /*!
     @brief Assign bond content

     Assigns type \c tp with list of Qnum \c qnums to Bond, replacing the current content.
     @param tp Type of bond
     @param qnums  Vector of Qnum
     */
    void assign(bondType tp, const std::vector<Qnum>& qnums);

    /// @brief Access bond type
    ///
    /// Returns the bond type of Bond
    /// @return Type of bond
    bondType type()const;

    /// @brief Access bond dimension
    ///
    /// Returns the dimension of Bond
    /// @return Dimension of Bond
    int dim()const;

    /// @brief Returns the degeneracies for quantum numbers
    /// @return  Map of Qnum's to their degeneracies
    std::map<Qnum, int> degeneracy()const;

    /// @brief Access quantum numbers
    ///
    /// Returns a vector of Qnum's for states in Bond. The size of the vector is the same
    /// as the dimension of Bond.
    /// @return Vector of Qnum
    std::vector<Qnum> Qlist()const;

    /// @brief Change the type of Bond
    ///
    /// Changes the type of Bond and the Qnum's of the bond when necssary.
    /// If bond is changed from incoming to outgoing or vice versa,
    /// the Qnum's are changed to - Qnum's
    /// @param tp %Bond type to change to
    Bond& change(bondType tp);

    /// @brief Dummy Change type of Bond
    /// @param tp %Bond type to change to
    Bond& dummy_change(bondType tp);
    /// @brief Combine bonds
    ///
    /// Combines Bond with another bond \c bd,  and expands the bond dimension by the direct product
    /// of Qnum's of two bonds. The resulting bond type is unchanged.
    ///
    /// @param bd Bond to be combined
    /// @return \c *this
    Bond& combine(Bond bd);

    /// @brief Compare two bonds
    ///
    /// The equality condition is such that
    /// \c b1.type() == \c b2.type() && \c b1.Qlist() == \c b2.Qlist()
    /// @param b1, b2 the bonds to compare
    /// @return \c True if two bonds are equal; \c False otherwise
    friend bool operator== (const Bond& b1, const Bond& b2);

    /// @brief Combines a list of bonds
    ///
    /// Combines a list of bonds by successively combining the bonds in the order in the given list \c bds.
    ///
    /// @return A bond with type \c tp
    friend Bond combine(bondType tp, const std::vector<Bond>& bds);

    /// @brief Combines a list of bonds
    /// @return A bond with bond type of the first bond in \c bds
    friend Bond combine(const std::vector<Bond>& bds);

    /// @brief Print out Qnum
    ///
    /// Prints out a bond  as:
    ///
    /// @code
    /// IN : (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|2, (U1 = -1, P = 0, 0)|1, Dim = 4
    /// @endcode
    /// The bond is an incoming bond with three Qnum's: one \c U1 =1, two \c U1 =0 and one \c U1= -1.
    /// The dimension of the bond is 4.
    friend std::ostream& operator<< (std::ostream& os, const Bond& b);

    friend class UniTensor;
    friend class CUniTensor;
    friend class Node;
    friend class CNode;
private:
    void setting(const std::vector<Qnum>& qnums);
    bondType m_type;
    int m_dim;
    std::vector<Qnum>Qnums; //Quantum numbers
    std::vector<int>Qdegs;  //Degeneracy in each quantum sector
    std::vector<int>offsets;
};

Bond combine(bondType tp, const std::vector<Bond>& bds);

Bond combine(const std::vector<Bond>& bds);


/// @example egB1.cpp

};


#endif /* BOND_H */

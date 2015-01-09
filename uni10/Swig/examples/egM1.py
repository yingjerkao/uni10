import sys
import pyUni10 as uni10
import copy

# Spin 1/2 Heisenberg hamiltonian
elem = [1.0/4,      0,      0,     0,\
				    0, -1.0/4,  1.0/2,     0,\
				    0,  1.0/2, -1.0/4,     0,\
				    0,      0,      0, 1.0/4]

H = uni10.Matrix(4, 4, elem)

# Diagonlize H
print H
results = H.eigh()
print "The eigen values:\n\n", results[0]
print "The eigen vectors:\n\n", results[1]

# Access element in a diagonal matrix
D = results[0];
print "D[2] =", D[2]

# Assign element
print "\nAssign D[3] = 7\n"
D[3] = 7
print D

# Access element
print "H[5] =", H[5]


#Make a pure density matrix from ground state
U = results[1];
# Generate ground state by taking the first H.rol() elements from U.
GS = uni10.Matrix(1, H.col(), U.getElem());
# Transpose GS
GST = copy.copy(GS)
GST.transpose()

Rho = GST * GS
print "\nPure density matrix of the ground state:\n"
print Rho

print "\nThe ground state energy: \n", (Rho * H).trace()

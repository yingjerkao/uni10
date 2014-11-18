import sys
import pyUni10 as uni10

# Construct spin 1 Heisenberg model by reading in the tensor which is written out in example egU1
H_U1 = uni10.UniTensor("egU1_H_U1")

# Get the block of quantum number q0 as a matrix "block0"
q0 = uni10.Qnum(0)
block0 = H_U1.getBlock(q0)
print block0

# Permute bonds by its labels, the default labels are [0, 1, 2, 3]
#Permute bonds to which with labels [1, 2, 3, 0] and leaving 1 bond as in-coming bonds.
H_U1.permute([1, 2, 3, 0], 1)


# combine the two bonds with label 2 and 3
H_U1.combineBond([2, 3])
print H_U1

import sys
import pyUni10 as uni10
import copy

# Construct spin 1 Heisenberg model by reading in the tensor which is written out in example egU1
H_U1 = uni10.UniTensor("egU1_H_U1")

# Randomly create an isometry tensor
q0 = uni10.Qnum(0)
q1 = uni10.Qnum(1)
q_1 = uni10.Qnum(-1)
q2 = uni10.Qnum(2)
q_2 = uni10.Qnum(-2)

bdi = uni10.Bond(uni10.BD_IN, [q2, q1, q0, q0, q_1, q_2]) # Do truncation here
bdo = uni10.Bond(uni10.BD_OUT, [q1, q0, q_1]) # Physical dimension

# Create isometry tensor W and transposed WT
W = uni10.UniTensor([bdi, bdo, bdo], "W");
W.orthoRand();
WT = copy.copy(W)
WT.transpose()

# Operate W and WT on H_U1, see the contraction labels in the documentation.
H_U1.setLabel([1, 2, 3, 4])
W.setLabel([-1, 1, 2])
WT.setLabel([3, 4, -2])
print W * H_U1 * WT;

# Write the tensors W and WT out to file
W.save("egU3_W");
WT.save("egU3_WT");

# Check the memory usage
uni10.UniTensor.profile()

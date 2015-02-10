import sys
import pyUni10 as uni10

M = uni10.Matrix(4, 5)
M.randomize()
print M
# carry out SVD
rets = M.svd();

# write matrice out to file
rets[0].save("mat_U")
rets[1].save("mat_Sigma")
rets[2].save("mat_VT")

# read in the matrice we've just written out
U = uni10.Matrix(rets[0].row(), rets[0].col(), rets[0].isDiag());
S = uni10.Matrix(rets[1].row(), rets[1].col(), rets[1].isDiag());
VT = uni10.Matrix(rets[2].row(), rets[2].col(), rets[2].isDiag());

U.load("mat_U");
S.load("mat_Sigma");
VT.load("mat_VT");
print U, S, VT
print U * S * VT

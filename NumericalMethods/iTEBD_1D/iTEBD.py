import pyUni10 as uni10
import sys
import numpy as np
if("../Hamiltonians" not in sys.path):
  sys.path.append("../Hamiltonians")
import hamiltonian as ham

def bondcat(T, L, bidx):
	labels = T.label();
	per_labels = list(T.label())
	per_labels.insert(0, per_labels.pop(bidx))
	inBondNum = T.inBondNum();
	T.permute(per_labels, 1)
	T.putBlock(L * T.getBlock())
	T.permute(labels, inBondNum);

def bondrm(T, L, bidx):
	invL = uni10.Matrix(L.row(), L.col(), True)
	for i in xrange(L.elemNum()):
		invL[i] = 0 if L[i] == 0 else 1.0 / L[i]
	bondcat(T, invL, bidx)

chi = 30
delta = 0.02
N = 10000
H = ham.Heisenberg()

bdi_chi = uni10.Bond(uni10.BD_IN, chi);
bdo_chi = uni10.Bond(uni10.BD_OUT, chi);
Gs = []
Gs.append(uni10.UniTensor([bdi_chi, bdo_chi, H.bond(2)], "Ga"))
Gs.append(uni10.UniTensor([bdi_chi, bdo_chi, H.bond(2)], "Gb"))
Gs[0].randomize(), Gs[1].randomize()

Ls = []
Ls.append(uni10.Matrix(chi, chi, True))  # Diagonal matrix
Ls.append(uni10.Matrix(chi, chi, True))  # Diagonal matrix
Ls[0].randomize(), Ls[1].randomize()

U = uni10.UniTensor(H.bond(), "U");
U.putBlock(uni10.takeExp(-delta, H.getBlock()))

for step in range(N):
	# Construct theta
	A = step % 2
	B = (step + 1) % 2
	bondcat(Gs[A], Ls[A], 1);
	bondcat(Gs[A], Ls[B], 0);
	bondcat(Gs[B], Ls[B], 1);
	Gs[A].setLabel([-1, 3, 1]);
	Gs[B].setLabel([3, -3, 2]);
	U.setLabel([1, 2, -2, -4]);
	theta = uni10.contract(Gs[A], Gs[B], True) # Gs[A], Gs[B] is permuted atfer the execution
	theta *= U;
	theta.permute([-1, -2, -3, -4], 2);

	# SVD
	svd = theta.getBlock().svd()

	# Truncation
	sv = svd[1]
	norm = sv.resize(chi, chi).norm()
	sv *= 1.0 / norm;
	Ls[A] = sv
	Gs[A].putBlock(svd[0].resize(svd[0].row(), chi));
	Gs[B].putBlock(svd[2].resize(chi, svd[2].col()));
	Gs[A].permute([-1, 3, 1], 1);
	bondrm(Gs[A], Ls[B], 0);
	bondrm(Gs[B], Ls[B], 1);

val = (theta * theta)[0]
print "E =", -np.log(val) / delta / 2

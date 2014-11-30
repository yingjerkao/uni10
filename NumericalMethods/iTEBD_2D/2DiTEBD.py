import pyUni10 as uni10
import sys
import numpy as np
import copy
if("../Hamiltonians" not in sys.path):
  sys.path.append("../Hamiltonians")
import hamiltonian as ham

def bondcat(T, L, bidx):
	labels = T.label()
	per_labels = list(T.label())
	per_labels.insert(0, per_labels.pop(bidx))
	inBondNum = T.inBondNum()
	T.permute(per_labels, 1)
	T.putBlock(L * T.getBlock())
	T.permute(labels, inBondNum)
	return T

def bondrm(T, L, bidx):
	invL = uni10.Matrix(L.row(), L.col(), True)
	for i in xrange(L.elemNum()):
		invL[i] = 0 if L[i] == 0 else 1.0 / L[i]
	return bondcat(T, invL, bidx)

def sqrt_bondcat(T, L, bidx):
	sqL = uni10.Matrix(L.row(), L.col(), True)
	for i in xrange(L.elemNum()):
		sqL[i] = np.sqrt(L[i])
	return bondcat(T, sqL, bidx);

def updateU(AL, BL, LU, LR, LD, LL, U, iTEBD, updateA):
	BL = bondcat(BL, LR, 4);
	iTEBD.putTensor("A", AL)
	iTEBD.putTensor("B", BL);
	iTEBD.putTensor("expH", U);
	C = iTEBD.launch();
	Theta = copy.copy(C);
	Theta = bondcat(Theta, LD, 2)
	Theta = bondcat(Theta, LL, 3)

	svds = Theta.getBlock().svd();
	dim = LU.row()
	LU = svds[1];
	norm = LU.resize(dim, dim).norm()
	LU = LU * (1.0 / norm)
	BL.assign(BL.bond());
	BL.permute([3, 0, 4, 1, 2], 1);
	BL.putBlock(svds[2].resize(dim, svds[2].col()));
	BL.permute(range(5), 1);

	updateA.putTensor("B", BL);
	updateA.putTensor("C", C);
	AL = updateA.launch() * (1 / norm);
	AL.setLabel(range(5));
	BL = bondrm(BL, LR, 4);
	return AL, BL, LU;

def updateR(AL, BL, LU, LR, LD, LL, U, iTEBD, updateA):
	BL = bondcat(BL, LU, 3);
	iTEBD.putTensor("A", AL)
	iTEBD.putTensor("B", BL);
	iTEBD.putTensor("expH", U);
	C = iTEBD.launch();
	Theta = copy.copy(C);
	Theta = bondcat(Theta, LD, 1)
	Theta = bondcat(Theta, LL, 2)

	svds = Theta.getBlock().svd();
	dim = LR.row()
	LR = svds[1]
	norm = LR.resize(dim, dim).norm();
	LR = LR * (1 / norm);
	BL.assign(BL.bond());
	BL.permute([4, 0, 1, 2, 3], 1);
	BL.putBlock(svds[2].resize(dim, svds[2].col()));
	BL.permute(range(5), 1);

	updateA.putTensor("B", BL);
	updateA.putTensor("C", C);
	AL = updateA.launch();
	AL *= (1 / norm);
	AL.setLabel(range(5));
	BL = bondrm(BL, LU, 3);
	return AL, BL, LR

def updateD(AL, BL, LU, LR, LD, LL, U, iTEBD, updateA):
	BL, AL, LD = updateU(BL, AL, LD, LL, LU, LR, U, iTEBD, updateA);
	return AL, BL, LD


def updateL(AL, BL, LU, LR, LD, LL, U, iTEBD, updateA):
	BL, AL, LL = updateR(BL, AL, LD, LL, LU, LR, U, iTEBD, updateA);
	return AL, BL, LL

def bond_expectation(ALL, BLL, Ob, measure_nt, norm_net):
	measure_net.putTensor("ALL", ALL);
	measure_net.putTensor("BLL", BLL);
	measure_net.putTensorT("ALLT", ALL);
	measure_net.putTensorT("BLLT", BLL);
	measure_net.putTensorT("Ob", Ob);
	val = measure_net.launch()

	norm_net.putTensor("ALL", ALL);
	norm_net.putTensor("BLL", BLL);
	norm_net.putTensorT("ALLT", ALL);
	norm_net.putTensorT("BLLT", BLL);
	norm = norm_net.launch();
	return val[0] / norm[0];


def measure(AL, BL, LU, LR, LD, LL, Ob, measure_net, norm_net):
	ALL = copy.copy(AL)
	BLL = copy.copy(BL)
	bondcat(ALL, LD, 3);
	bondcat(ALL, LL, 4);
	bondcat(BLL, LU, 3);
	bondcat(BLL, LR, 4);
	Ls = [LU, LR, LD, LL]

	val = 0.0
	for i in xrange(4):
		BLL1 = copy.copy(BLL);
		rotate_label = [x if x == 0 else ((x-1) + i) % 4 + 1 for x in xrange(5)]
		ALL.permute(rotate_label, 1);
		BLL1.permute(rotate_label, 1);
		BLL1 = bondrm(BLL1, Ls[i], 3);
		val += bond_expectation(ALL, BLL1, Ob, measure_net, norm_net);
	return val / 4.0;


chi = 2
N = 1000
tau = 0.01
H = ham.transverseIsing(0.7)
print H

bdi_chi = uni10.Bond(uni10.BD_IN, chi);
bdo_chi = uni10.Bond(uni10.BD_OUT, chi);

GaL = uni10.UniTensor([H.bond(0), bdo_chi, bdo_chi, bdo_chi, bdo_chi], "GaL"); #gammaA * lambda
GbL = uni10.UniTensor([H.bond(0), bdo_chi, bdo_chi, bdo_chi, bdo_chi], "GbL"); #gammaA * lambda
GaL.randomize();
GbL.randomize();

I_chi = uni10.Matrix(chi, chi, True);
I_chi.randomize();
LU = copy.copy(I_chi) # up
LR = copy.copy(I_chi) # right
LD = copy.copy(I_chi) # down
LL = copy.copy(I_chi) #left

GaL = bondcat(GaL, LU, 1);
GaL = bondcat(GaL, LR, 2);
GbL = bondcat(GbL, LD, 1);
GbL = bondcat(GbL, LL, 2);

U = uni10.UniTensor(H.bond(), "U")
U.putBlock(uni10.takeExp(-tau, H.getBlock()));

iTEBD_V_net = uni10.Network("2DiTEBD_V.net");
updateA_V_net = uni10.Network("updateA_V.net");
iTEBD_H_net = uni10.Network("2DiTEBD_H.net")
updateA_H_net = uni10.Network("updateA_H.net");
measure_net = uni10.Network("measure.net");
norm_net = uni10.Network("norm.net");

for i in xrange(N):
	GaL, GbL, LU = updateU(GaL, GbL, LU, LR, LD, LL, U, iTEBD_V_net, updateA_V_net);
	GaL, GbL, LR = updateR(GaL, GbL, LU, LR, LD, LL, U, iTEBD_H_net, updateA_H_net);
	GaL, GbL, LD = updateD(GaL, GbL, LU, LR, LD, LL, U, iTEBD_V_net, updateA_V_net);
	GaL, GbL, LL = updateL(GaL, GbL, LU, LR, LD, LL, U, iTEBD_H_net, updateA_H_net);
	E = measure(GaL, GbL, LU, LR, LD, LL, H, measure_net, norm_net)
	print "E =", E

Mz = uni10.UniTensor(H.bond())
Mx = uni10.UniTensor(H.bond())
I = uni10.Matrix(H.bond(0).dim(), H.bond(0).dim());
I.identity();
Mx.putBlock(uni10.otimes(2*ham.matSx(), I));
Mz.putBlock(uni10.otimes(2*ham.matSz(), I));

mx = abs(measure(GaL, GbL, LU, LR, LD, LL, Mx, measure_net, norm_net))
mz = abs(measure(GaL, GbL, LU, LR, LD, LL, Mz, measure_net, norm_net))
print mx, mz

GaL = bondrm(GaL, LU, 1);
GaL = bondrm(GaL, LR, 2);
GbL = bondrm(GbL, LD, 1);
GbL = bondrm(GbL, LL, 2);

GaL = sqrt_bondcat(GaL, LU, 1);
GaL = sqrt_bondcat(GaL, LR, 2);
GaL = sqrt_bondcat(GaL, LD, 3);
GaL = sqrt_bondcat(GaL, LL, 4);

GbL = sqrt_bondcat(GbL, LD, 1);
GbL = sqrt_bondcat(GbL, LL, 2);
GbL = sqrt_bondcat(GbL, LU, 3);
GbL = sqrt_bondcat(GbL, LR, 4);

GaL.save("Data/GaL.ten");
GbL.save("Data/GbL.ten");

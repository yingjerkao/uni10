import pyUni10 as uni10
import os
import sys
import numpy as np
import copy

def makeT(GL):
	GLT = copy.copy(GL);
	GLT.transpose();
	GL.setLabel([0, 1, 2, 3, 4]);
	GLT.setLabel([-1, -2, -3, -4, 0]);
	T = GLT * GL
	for i in xrange(1, 5):
		T.combineBond([-i, i]);
	return T

def takeSVD(T, chi):
		bdi_chi = uni10.Bond(uni10.BD_IN, chi)
		bdo_chi = uni10.Bond(uni10.BD_OUT, chi)
		svd = T.getBlock().svd()
		S0 = uni10.UniTensor([T.bond(0), T.bond(1), bdo_chi])
		S1 = uni10.UniTensor([bdi_chi, T.bond(2), T.bond(3)])
		chi = chi if chi < svd[1].row() else svd[1].row()
		svd[1].resize(chi, chi);
		for i in xrange(chi):
			svd[1][i] = np.sqrt(svd[1][i])
		S0.putBlock(svd[0].resize(svd[0].row(), chi) * svd[1]);
		S1.putBlock(svd[1] * svd[2].resize(chi, svd[2].col()));
		return S0, S1

def TRG_contract(Ss, TRG_net):
	for i in xrange(4):
		TRG_net.putTensor(i, Ss[i])
	return TRG_net.launch()

def updateT(Ts, chi, TRG_net):
	Ss = ()
	for step in xrange(2): # 2 types of directions for svd
		T = Ts[step] if step < len(Ts) else Ts[step % len(Ts)]
		if step == 0:
			T.permute([-4, -3, -1, -2], 2);
		else:
				T.permute([-1, -4, -2, -3], 2);
		Ss += takeSVD(T, chi)
	return [TRG_contract(Ss, TRG_net)]

def updateAll(Ts, Imps, TRG_net):
	pass


def takeExpectation(Ts, exp_net):
	for step in xrange(4):
		T = Ts[step] if step < len(Ts) else Ts[step % len(Ts)]
		exp_net.putTensor(step, T);
	return exp_net.launch()[0]

root_dir = "../iTEBD_2D/Data"
GaL = uni10.UniTensor(os.path.join(root_dir, "GaL.ten"));
GbL = uni10.UniTensor(os.path.join(root_dir, "GbL.ten"));
TRG_net = uni10.Network("TRG.net");
exp_net = uni10.Network("expectation.net");

Ts = [makeT(GaL), makeT(GbL)]
chi = Ts[0].bond(0).dim()
print "E:", takeExpectation(Ts, exp_net);
for i in xrange(30):
	oriT = copy.copy(Ts[0]);
	Ts = updateT(Ts, 4, TRG_net)
	diff = (oriT.getBlock() + (-1 * Ts[0].getBlock())).norm()
	print "E:", takeExpectation(Ts, exp_net);


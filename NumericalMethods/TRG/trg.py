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

def makeImpurity(GL, Ob):
	GLT = copy.copy(GL);
	GLT.transpose();
	Ob.setLabel([0, 5])
	GL.setLabel([5, 1, 2, 3, 4])
	GLT.setLabel([-1, -2, -3, -4, 0]);
	I = Ob * GL;
	I = GLT * I
	for i in xrange(1, 5):
		I.combineBond([-i, i]);
	return I

def trgSVD(wch, T, chi):
		if wch % 2 == 0:
			T.permute([-4, -3, -1, -2], 2);
		else:
			T.permute([-1, -4, -2, -3], 2);
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

def trgContract(Ss, TRG_net):
	for i in xrange(4):
		TRG_net.putTensor(i, Ss[i])
	return TRG_net.launch()

def updateT(Ts, chi, TRG_net):
	S = []
	S.append(trgSVD(0, Ts[0], chi))
	S.append(trgSVD(1, Ts[1%len(Ts)], chi))
	a, b = 0, 1
	return [trgContract([S[a][0], S[b][0], S[a][1], S[b][1]], TRG_net)]

def updateAll(Ts, Imps, chi, TRG_net):
	S = []
	I = []
	for r in xrange(2):
		S.append(trgSVD(r, Ts[r%len(Ts)], chi)); # get [(Ta1, Ta2), (Tb1, Tb2)]
	for r in xrange(4):
		I.append(trgSVD(r, Imps[r], chi)); # get [(Ia1, Ia2), (Ib1, Ib2), (Ic1, Ic2), (Id1, Id2)]
	a, b, c, d = 0, 1, 2, 3
	Ts = [trgContract([S[a][0], S[b][0], S[a][1], S[b][1]], TRG_net)]
	Imps = [
			trgContract([S[a][0], I[b][0], I[a][1], S[b][1]], TRG_net),\
			trgContract([S[a][0], S[b][0], I[c][1], I[b][1]], TRG_net),\
			trgContract([I[c][0], S[b][0], S[a][1], I[d][1]], TRG_net),\
			trgContract([I[a][0], I[d][0], S[a][1], S[b][1]], TRG_net)]
	return Ts, Imps

def trgExpectation(Ts, exp_net):
	for step in xrange(4):
		exp_net.putTensor(step, Ts[step % len(Ts)]);
	return exp_net.launch()[0]

root_dir = "../iTEBD_2D/Data"
GaL = uni10.UniTensor(os.path.join(root_dir, "GaL.ten"));
GbL = uni10.UniTensor(os.path.join(root_dir, "GbL.ten"));
TRG_net = uni10.Network("TRG.net");
exp_net = uni10.Network("expectation.net");

bdi = uni10.Bond(uni10.BD_IN, 2);
bdo = uni10.Bond(uni10.BD_OUT, 2);
Mz = uni10.UniTensor([bdi, bdo], "Mz");
Mz.setElem([1, 0, 0, -1])

Ts = [makeT(GaL), makeT(GbL)]
Imps = [makeT(GaL), makeT(GbL), makeImpurity(GaL, Mz), makeT(GbL)]
chi = Ts[0].bond(0).dim()

print "<O>:", trgExpectation(Imps, exp_net) / trgExpectation(Ts, exp_net),
print "norm:", trgExpectation(Ts, exp_net)
for i in xrange(25):
	Ts, Imps = updateAll(Ts, Imps, chi, TRG_net);
	print "<O>:", trgExpectation(Imps, exp_net) / trgExpectation(Ts, exp_net),
	print "norm:", trgExpectation(Ts, exp_net)

import pyUni10 as uni10

def matSp():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [0, 1, 0, 0]);

def matSm():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [0, 0, 1, 0]);

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [0.5, 0, 0, -0.5]);

def Heisenberg():
	spin = 0.5
	sp = matSp();
	sm = matSm();
	sz = matSz();
	ham = uni10.otimes(sz, sz);
	ham += 0.5 * (uni10.otimes(sp, sm) + uni10.otimes(sm, sp));
	dim = int(spin * 2 + 1)
	bdi = uni10.Bond(uni10.BD_IN, dim);
	bdo = uni10.Bond(uni10.BD_OUT, dim);
	H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg");
	H.putBlock(ham)
	return H

def Heisenberg_U1():
	spin = 0.5
	q1 = uni10.Qnum(1);
	bdi = uni10.Bond(uni10.BD_IN, [q1, -q1]);
	bdo = uni10.Bond(uni10.BD_OUT, [q1, -q1]);
	H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
	print Heisenberg().getBlock()
	H.setRawElem(Heisenberg().getBlock().getElem());
	return H

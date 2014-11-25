import pyUni10 as uni10
import sys
import numpy as np
if("../Hamiltonians" not in sys.path):
  sys.path.append("../Hamiltonians")
import hamiltonian as ham
import copy

def bondcat(T, L, bidx):
	labels = T.label()
	inBondNum = T.inBondNum();
	if bidx == 0:
		# -- L -- T --
		#         |
		L.setLabel([labels[0], 77])
		T.setLabel([77, labels[1], labels[2]])
		T = L * T;
	elif bidx == 1:
		# -- T -- L --
		#    |
		L.setLabel([77, labels[1]])
		T.setLabel([labels[0], 77, labels[2]])
		T *= L;
	else:
		raise Exception("In bondcat(): Invalid Usage.")
	T.permute(labels, inBondNum);
	return T;

def bondrm(T, L, bidx):
  invL = copy.copy(L);
  qnums = invL.blockQnum();
  for qnum in qnums:
    blk = invL.getBlock(qnum, True); # get diagonal block
    for i in xrange(blk.elemNum()):
      blk[i] = 0 if blk[i] == 0 else 1.0 / blk[i];
    invL.putBlock(qnum, blk);
  return bondcat(T, invL, bidx);

def sv_merge(svs, bidxs, bidx, sv_mat, chi):
	if(len(svs)):
		length = len(svs) + sv_mat.elemNum();
		length = length if length < chi else chi
		ori_svs = svs
		ori_bidxs = bidxs
		svs = [0] * length
		bidxs = [0] * length
		cnt = 0;
		cur1 = 0;
		cur2 = 0
		while cnt < length:
			if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
				if ori_svs[cur1] >= sv_mat[cur2]:
					svs[cnt] = ori_svs[cur1]
					bidxs[cnt] = ori_bidxs[cur1];
					cur1 += 1
				else:
					svs[cnt] = sv_mat[cur2]
					bidxs[cnt] = bidx;
					cur2 += 1
			elif cur2 < sv_mat.elemNum():
				for i in xrange(cnt, len(svs)):
					svs[i] = sv_mat[cur2];
					bidxs[i] = bidx;
					cur2 += 1
				break
			else:
				svs[cnt:] = ori_svs[cur1 : (cur1 + len(svs) - cnt)];
				bidxs[cnt:] = ori_bidxs[cur1 : (cur1 + len(svs) - cnt)];
				break;
			cnt += 1;
	else:
		bidxs = [bidx] * sv_mat.elemNum()
		svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())];
	return svs, bidxs

def setTruncation(theta, GA, GB, LA, chi):
	svds = {};
	blk_qnums = theta.blockQnum();
	for qnum in blk_qnums:
		svds[qnum] = theta.getBlock(qnum).svd();
	svs = []
	bidxs = []
	for bidx in xrange(len(blk_qnums)):
		svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi)
	dims = [0] * len(blk_qnums);
	for bidx in bidxs:
		dims[bidx] += 1;
	qnums = []
	for bidx in xrange(len(blk_qnums)):
		qnums += [blk_qnums[bidx]] * dims[bidx];
	bdi_mid = uni10.Bond(uni10.BD_IN, qnums);
	bdo_mid = uni10.Bond(uni10.BD_OUT, qnums);
	labelGa = GA.label();
	GA.assign([GA.bond(0), GA.bond(1), bdo_mid]);
	GB.assign([bdi_mid, GB.bond(1), GB.bond(2)]);
	LA.assign([bdi_mid, bdo_mid]);
	GA.setLabel(labelGa);
	degs = bdi_mid.degeneracy();
	sv_mat = uni10.Matrix(bdi_mid.dim(), bdo_mid.dim(), svs, True);
	norm = sv_mat.norm();
	for qnum, dim in degs.iteritems():
		if qnum not in svds:
			raise Exception("In setTruncaton(): Fatal error.")
		svd = svds[qnum]
		GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim));
		GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()));
		LA.putBlock(qnum, svd[1].resize(dim, dim) * (1/norm));


chi = 40
delta = 0.1
N = 2000
H = ham.Heisenberg_U1()

Gs = []
bdi_mid = uni10.combine([H.bond(0), H.bond(0)]);
bdo_mid = uni10.combine([H.bond(2), H.bond(2)]);
Gs.append(uni10.UniTensor([H.bond(0), bdo_mid, H.bond(2)], "GA"))
Gs.append(uni10.UniTensor([bdi_mid, H.bond(2), H.bond(2)], "GB"))
Gs[0].randomize(), Gs[1].randomize();

Ls = []
Ls.append(uni10.UniTensor([bdi_mid, bdo_mid], "LA"))
Ls.append(uni10.UniTensor([H.bond(0), H.bond(2)], "LB"))
Ls[0].randomize(), Ls[1].randomize()

U = uni10.UniTensor(H.bond(), "U");
blk_qnums = H.blockQnum();
for qnum in blk_qnums:
	U.putBlock(qnum, uni10.takeExp(-delta, H.getBlock(qnum)));

for step in xrange(N):
	A = step % 2;
	B = (step + 1) % 2;
	Gs[A] = bondcat(Gs[A], Ls[A], 1);
	Gs[A] = bondcat(Gs[A], Ls[B], 0);
	Gs[B] = bondcat(Gs[B], Ls[B], 1);
	Gs[A].setLabel([-1, 3, 1]);
	Gs[B].setLabel([3, -3, 2]);
	U.setLabel([1, 2, -2, -4]);
	theta = uni10.contract(Gs[A], Gs[B], True) # Gs[A], Gs[B] is permuted atfer the execution
	theta *= U;
	theta.permute([-1, -2, -3, -4], 2);
  # update Gs[A], Gs[B], Ls[A]
	setTruncation(theta, Gs[A], Gs[B], Ls[A], chi);
	Gs[A].permute([-1, 3, 1], 1);
	Gs[A] = bondrm(Gs[A], Ls[B], 0);
	Gs[B] = bondrm(Gs[B], Ls[B], 1);

val = (theta * theta)[0]
print "E =", -np.log(val) / delta / 2

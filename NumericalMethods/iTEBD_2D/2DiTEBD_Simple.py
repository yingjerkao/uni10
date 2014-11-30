import pyUni10 as uni10
import sys
import numpy as np
import copy

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


def transverseIsing(h):
    spin = 0.5
    sx = 0.5*(matSp()+matSm())
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    ham =uni10.otimes(2*sz,2*sz)*(-1) +0.5*float(h)*(uni10.otimes(iden,2*sx)+uni10.otimes(2*sx,iden))
    dim = int(spin * 2 + 1)
    bdi = uni10.Bond(uni10.BD_IN, dim);
    bdo = uni10.Bond(uni10.BD_OUT, dim);
    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "TFIM");
    H.putBlock(ham)
    return H





def bondcat(T, L, bidx):
    labels = T.label();
    per_labels = list(T.label())
    per_labels.insert(0, per_labels.pop(bidx))
    inBondNum = T.inBondNum();
    T.permute(per_labels, 1)
    T.putBlock(L * T.getBlock())
    T.permute(labels, inBondNum);
    return T

def bondrm(T, L, bidx):
    invL = uni10.Matrix(L.row(), L.col(), True)
    for i in xrange(L.elemNum()):
        invL[i] = 0 if L[i] == 0 else 1.0 / L[i]
    return bondcat(T, invL, bidx)

chi = 30
delta = 0.001
Ns = 1000
H = transverseIsing(0.8)
#H=Heisenberg()

z=2

bdi_chi = uni10.Bond(uni10.BD_IN, chi);
bdo_chi = uni10.Bond(uni10.BD_OUT, chi);
bond_list=[]
for i in range(z):
    bond_list.append(bdi_chi)
bond_list.append(H.bond(2))

G = []
G.append(uni10.UniTensor(bond_list, "Ga"))
G.append(uni10.UniTensor(bond_list, "Gb"))
G[0].randomize(), G[1].randomize()


I_chi = uni10.Matrix(chi, chi, True)
I_chi.randomize()
L = []
for i in range(z):
    L.append(uni10.Matrix(chi, chi, True))  # Diagonal matrix
    L[i]=copy.copy(I_chi)




def update(Gs, Ls, z, delta, H0,N):
    """
    Performs two-site updates
    """
    U = uni10.UniTensor(H.bond(), "U");
    U.putBlock(uni10.takeExp(-delta, H0.getBlock()))
    for step in range(N):
        E=0.0
        for bond in range(z):
            Gs_A_l=copy.copy(Gs[0])
            Gs_A_l.setName("Gs_Al")

            for b in range(1,z):
                Gs_A_l=bondcat(Gs_A_l, Ls[(bond+b) % z], ((bond+b) % z))

            Gs_A_l=bondcat(Gs_A_l, Ls[bond], bond)

            per_labels = list(Gs_A_l.label())
            per_labels.insert(0, per_labels.pop(bond))
            inBondNum = Gs_A_l.inBondNum();
            Gs_A_l.permute(per_labels, inBondNum)

            A_merged_bonds=[ (bond+i)%z for i in range(1,z)]

            if len(A_merged_bonds) > 1 :
                Gs_A_l.combineBond(A_merged_bonds)


            Gs_B_l=copy.copy(Gs[1])
            Gs_B_l.setName("Gs_Bl")


            for b in range(1,z):
                Gs_B_l=bondcat(Gs_B_l, Ls[(bond+b) % z], (bond+b) % z)


            per_labels = list(Gs_B_l.label())
            per_labels.insert(0, per_labels.pop(bond))
            inBondNum = Gs_B_l.inBondNum();
            Gs_B_l.permute(per_labels, inBondNum)

            B_merged_bonds=[ (bond+i)%z for i in range(1,z)]

            if len(B_merged_bonds) > 1:
                Gs_B_l.combineBond(B_merged_bonds)

            Gs_A_l.setLabel([3, -1, 1]);
            Gs_B_l.setLabel([3, -3 , 2]);
            U.setLabel([1, 2, -2, -4]);

            theta = uni10.contract(Gs_A_l, Gs_B_l, True) # Gs[0], Gs[1] is permuted atfer the execution
            Ntheta=theta
            theta *= U;
            theta.permute([-1, -2, -3, -4], 2);



            # SVD
            svd = theta.getBlock().svd()

            # Truncation
            sv = svd[1]
            norm = sv.resize(chi, chi).norm()
            sv = sv * (1.0 / norm);
            Ls[bond] = sv



            Gs_A_l.putBlock(svd[0].resize(svd[0].row(),chi));
            Gs_B_l.putBlock(svd[2].resize(chi,svd[2].col()));



            Gs_A_l.permute([3,-1,1],2)
            Gs_B_l.permute([3,-3,2],2)

            labels=list(Gs[0].label())

            per_labels=[bond]+A_merged_bonds+[z]

            Gs[0].permute(per_labels,z)
            Gs[0].putBlock(Gs_A_l.getBlock())

            Gs[0].permute(labels,z)

            labels=list(Gs[1].label())
            per_labels=[bond]+B_merged_bonds+[z]


            Gs[1].permute(per_labels,z)
            Gs[1].putBlock(Gs_B_l.getBlock())
            Gs[1].permute(labels,z)

            for j in A_merged_bonds:
                Gs[0] = bondrm(Gs[0], Ls[j], j);
            for j in B_merged_bonds:
                Gs[1] = bondrm(Gs[1], Ls[j], j);

            Ebond = (theta*theta)[0]
            Nbond = (Ntheta * Ntheta)[0]
            E_b=-np.log(Ebond)/delta/2/z
            E+=E_b/Nbond
        print "E=",E


update(G,L,delta,H,z,Ns)

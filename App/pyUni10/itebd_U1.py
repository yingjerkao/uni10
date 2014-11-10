#!/usr/bin/env python

# import modules
import math
import copy
import pyUni10 as uni10
from copy import *
# todo
# create dictionary
# U1 = dict([('q0', uni10.Qnum(0))])
# for i in range(1, 3):
#     qp = 'q+' + str(i)
#     qm = 'q-' + str(i)
#     U1[qp] = uni10.Qnum(+i)
#     U1[qm] = uni10.Qnum(-i)
#
#
# print U1
# print U1['q+1']
# exit()

# qlist = list of uni10.Qnum()
# q1 = uni10.Qnum(1)
# qlist=[q1, q2, q3]

# bond_in= uni10.Bond(uni10.BD_IN, qlist)
# bond_out= uni10.Bond(uni10.BD_OUT, qlist)
# bond_list = [bond_in, bond_out, bond_in]
# tensor = uni10.UniTensor(bond_list)

# set up quantum number for single-site operators
q0 = uni10.Qnum(0)
q1 = uni10.Qnum(1)
q_1 = uni10.Qnum(-1)
q2 = uni10.Qnum(2)
q_2 = uni10.Qnum(-2)

def QL_2_BL(Q_list):
    q_list = []
    for Q in Q_list:
        for i in range(Q[1]):
            q_list.append(Q[0])
    return q_list

# init a rank-2 tensor
Q_list = [(q0, 2), (q1, 3), (q_1, 3)]

def rank2tensor(Q_list):
    q_list = QL_2_BL(Q_list)
    bdi = uni10.Bond(uni10.BD_IN, q_list)
    bdo = uni10.Bond(uni10.BD_OUT, q_list)
    tensor = uni10.UniTensor([bdi, bdo])
    tensor.set_zero()
    return tensor

T = rank2tensor(Q_list)
print T

# S=1/2 XXZ model
# q=2S
# S=+1/2 => q=+1
# S=-1/2 => q=-1

# set up operators and Hamiltonian
Sp_raw = [0, 1,\
          0, 0]

Sm_raw = [0, 0,\
          1, 0]

Sz_raw = [+0.5, 0,\
          0, -0.5]

Id_raw = [1, 0,\
          0, 1]

# bonds without U(1) symmetry
bdi_nosym = uni10.Bond(uni10.BD_IN, [q0, q0])
bdo_nosym = uni10.Bond(uni10.BD_OUT, [q0, q0])

# bonds with U(1) symmetry
bdi = uni10.Bond(uni10.BD_IN, [q1, q_1])
bdo = uni10.Bond(uni10.BD_OUT, [q1, q_1])

# bond that represent the U(1) charge of the operator
bd_op1 = uni10.Bond(uni10.BD_OUT, [q2])
bd_op_1 = uni10.Bond(uni10.BD_OUT, [q_2])

# Single-site operators
# Sp_nosym
Sp_nosym = uni10.UniTensor([bdi_nosym, bdo_nosym])
Sp_nosym.setName('Sp_nosym')
Sp_nosym.setRawElem(Sp_raw)
print Sp_nosym

Sp = uni10.UniTensor([bdi, bd_op1, bdo])
Sp.setName('Sp')
Sp.setRawElem(Sp_raw)
print Sp

# Sm
Sm_nosym = uni10.UniTensor([bdi_nosym, bdo_nosym])
Sm_nosym.setName('Sm_nosym')
Sm_nosym.setRawElem(Sm_raw)
print Sm_nosym

Sm = uni10.UniTensor([bdi, bd_op_1, bdo])
Sm.setName('Sm')
Sm.setRawElem(Sm_raw)
print Sm

# Sz
Sz = uni10.UniTensor([bdi, bdo])
Sz.setName('Sz')
Sz.setRawElem(Sz_raw)
print Sz

#  Id
Id = uni10.UniTensor([bdi, bdo])
Id.setName('Id')
Id.setRawElem(Id_raw)
print Id

# two-sites operators
# SpSm_nosym
SpSm_nosym = uni10.otimes(Sp_nosym, Sm_nosym)
SpSm_nosym.setName('SpSm_nosym')
print SpSm_nosym

# SmSp_nosym
SmSp_nosym = uni10.otimes(Sm_nosym, Sp_nosym)
SmSp_nosym.setName('SmSp_nosym')
print SmSp_nosym

print SpSm_nosym + SmSp_nosym

# SpSm
SpSm = uni10.otimes(Sp, Sm)
SpSm.setName('SpSm')
# print SpSm
SpSm.combineBond([2, 3, 4])
SpSm.setLabel([-1, -2, 1, 2])
print SpSm

# SmSp
SmSp = uni10.otimes(Sm, Sp)
SmSp.setName('SmSp')
# print SmSp
SmSp.combineBond([2, 3, 4])
SmSp.setLabel([-1, -2, 1, 2])
print SmSp

print SpSm + SmSp

# SzSz
SzSz = uni10.otimes(Sz, Sz)
SzSz.setName('SzSz')
SzSz.setLabel([-1, -2, 1, 2])
print SzSz

# IdId
IdId = uni10.otimes(Id, Id)
IdId.setName('IdId')
IdId.setLabel([-1, -2, 1, 2])
print IdId

# set up Hamiltonian
# Ham = SpSm+SmSp+SzSz
Ham = 0.5 * (SpSm + SmSp) + SzSz
Ham.setName('Ham')
print Ham

####################

# set up Exp(-Ham*dt)
def Ham_Exp(Ham, dt):
    Ham_Exp_dt = copy(Ham)

    QList = Ham_Exp_dt.getBlocks()
    for q in QList:
        print 'q=', q
        M = Ham_Exp_dt.getBlock(q)
        (D, U) = M.eigh()
        UT=copy(U)
        UT.transpose()
        for i in range(D.row()):
            D[i]=math.exp(-D[i] * dt)
        Ham_Exp_dt.putBlock(q, UT * D * U)
    return Ham_Exp_dt


# initialize Gamma and Lambda
# virtual bonds with U(1) symmetry
# bdi_v = uni10.Bond(uni10.BD_IN, [q1, q0, q_1])
# bdo_v = uni10.Bond(uni10.BD_OUT, [q1, q0, q_1])
bdi_v = uni10.Bond(uni10.BD_IN, [q1, q1, q0, q0, q_1, q_1])
bdo_v = uni10.Bond(uni10.BD_OUT, [q1, q1, q0, q0, q_1, q_1])

# Gamma
Gamma_A = uni10.UniTensor([bdi_v, bdi, bdo_v])
Gamma_A.setName('Gamma_A')

Gamma_B = uni10.UniTensor([bdi_v, bdi, bdo_v])
Gamma_B.setName('Gamma_B')

Gamma_A.randomize()
Gamma_B.randomize()

print "Gamma A:", Gamma_A
print "Gamma B:", Gamma_B

# Lambda
Lambda_A = uni10.UniTensor([bdi_v, bdo_v])
Lambda_B = uni10.UniTensor([bdi_v, bdo_v])

# uniform distribution of lambda
Lambda_A.identity()
Lambda_B.identity()

# fix the normalization
Lambda_A = (1.0/3.0) * Lambda_A
Lambda_B = (1.0/3.0) * Lambda_B
Lambda_A.setName('Lambda_A')
Lambda_B.setName('Lambda_B')

print Lambda_A
print Lambda_B


# find inverse
def tensor_inv(tensor):
    tensor_inv = copy(tensor)
    for block in range(tensor_inv.blockNum()):
        q = tensor_inv.blockQnum(block)
        M = tensor_inv.getBlock(q)
        # print M.col(), M.row()
        for i in range(M.col()):
            M[i, i] = 1.0/M[i, i]
        tensor_inv.putBlock(q, M)
    return tensor_inv

Lambda_A_inv = tensor_inv(Lambda_A)
Lambda_A_inv.setName('Lambda_A_inv')
Lambda_B_inv = tensor_inv(Lambda_B)
Lambda_B_inv.setName('Lambda_B_inv')
print Lambda_B_inv

# setup dt and exp(-H*dt)
dt = 0.1
Ham_Exp_dt = Ham_Exp(Ham, dt)
Ham_Exp_dt.setName('Ham_Exp_dt')
print Ham_Exp_dt

########################################
# iTEBD update
# todo
# def iTEBD_update(dt):
#     Ham_Exp_dt = Ham_Exp(Ham, dt)
#     print Ham_Exp_dt
#     pass

# return (Gamma_L, Lambda_C, Gamma_R)
# (Gamma_L, Lambda_C, Gamma_R) = iTEBD_update(....)
def iTEBD_update(Ham_Expt_dt, Lambda_L, Gamma_L, Lambda_C, Gamma_R, Lambda_R):
    print ">>>>> iTEBD_update >>>>>"
    # construct S via network
    net = uni10.Network("t.net")
    net.putTensor("Ham_Exp_dt", Ham_Exp_dt)
    net.putTensor("Lambda_L", Lambda_L)
    net.putTensor("Gamma_L", Gamma_L)
    net.putTensor("Lambda_C", Lambda_C)
    net.putTensor("Gamma_R", Gamma_R)
    net.putTensor("Lambda_R", Lambda_R)
    S = net.launch()
    S.setName('S')
    print S

    # all eigenvalues, to be sorted
    Lambda_all = []
    # use dictionary to store the Lambda, U, VT for each block
    Lambda = {}
    U = {}
    VT = {}
    # use Q_list to keep track of (quantum numbers, degeneracy)
    Q_list =[]

    # SVD on all blocks
    for block in range(S.blockNum()):
        q = S.blockQnum(block)
        M_q = S.getBlock(q)
        Q = (q, M_q.col())
        Q_list.append(Q)
        (U_q, Lambda_q, VT_q) = M_q.svd()
        U[repr(q)] = U_q
        VT[repr(q)] = VT_q
        Lambda[repr(q)] = Lambda_q
        for i in range(M_q.col()):
            Lambda_all.append(Lambda_q[i])
    print '>>>>> Lambda=\n', Lambda
    print '>>>>> U=\n', U
    print '>>>>> VT=\n', VT
    print '>>>>> Q_list=\n', Q_list
    Lambda_all.sort(reverse=True)
    print '>>>>> Lambda_all=\n', Lambda_all

    # identify the smallest lambda to keep
    D_cut = 12
    if D_cut >= len(Lambda_all):
        Lambda_smallest = Lambda_all[-1]
    else:
        Lambda_smallest = Lambda_all[D_cut-1]
    print '>>>>> Lambda_smallest = \n', Lambda_smallest

    # copy the in and out bond of the S separately
    bdi_old_list = []
    bdo_old_list = []
    for block in range(S.bondNum()):
        bond = S.bond(block)
        if bond.type() == uni10.BD_IN :
            bdi_old_list.append( bond )
        else:
            bdo_old_list.append( bond )

    # perform the truncation
    print ">>>>> perform the truncation >>>>> "
    Q_list_cut = []
    Lambda_M_new = {}
    U_cut = {}
    VT_cut = {}

    for block in range(S.blockNum()):
        q = S.blockQnum(block)
        M_q = Lambda[repr(q)]
        Lambda_q_cut = []
        for i in range(M_q.col()):
            # print M_q[i], Lambda_cut, M_q[i] >= Lambda_cut
            if M_q[i] >= Lambda_smallest:
                Lambda_q_cut.append(M_q[i])
        print '>>>>> Lambda_q_cut=\n', Lambda_q_cut
        M_q_cut = uni10.Matrix(len(Lambda_q_cut), len(Lambda_q_cut), Lambda_q_cut, True)
        Lambda_M_new[repr(q)] = M_q_cut
        # U
        U_q = U[repr(q)]
        # print U_q
        U_q.transpose()
        # print U_q
        U_q_cut = uni10.Matrix(len(Lambda_q_cut), U_q.row(), U_q.getElem())
        U_q_cut.transpose()
        # print U_q_cut
        U_cut[repr(q)] = U_q_cut
        # VT
        VT_q = VT[repr(q)]
        VT_q_cut = uni10.Matrix(len(Lambda_q_cut), VT_q.col(), VT_q.getElem())
        print VT_q
        print VT_q_cut
        VT_cut[repr(q)] = VT_q_cut
        # Q
        Q = (q, len(Lambda_q_cut))
        Q_list_cut.append(Q)

    # print U_cut
    # print VT_cut
    # print Q_list_cut
    # print Lambda_M_new

    Lambda_new = rank2tensor(Q_list_cut)
    Lambda_new.setName('Lambda_new')
    for block in range(Lambda_new.blockNum()):
        q = Lambda_new.blockQnum(block)
        # print q
        # print Lambda_M_new[repr(q)]
        Lambda_new.putBlock(q, Lambda_M_new[repr(q)])
    print Lambda_new

    # new IN and OUT bound after the truncation
    q_new_list = QL_2_BL(Q_list_cut)
    bdi_new = uni10.Bond(uni10.BD_IN, q_new_list)
    bdo_new = uni10.Bond(uni10.BD_OUT, q_new_list)

    Lambda_Gamma_new = uni10.UniTensor(bdi_old_list + [bdo_new])
    Lambda_Gamma_new.setName('Lambda_Gamma_new')
    Lambda_Gamma_new.set_zero()
    for block in range(Lambda_Gamma_new.blockNum()):
        q = Lambda_Gamma_new.blockQnum(block)
        # print 'block:', block, q
        # print U_cut[repr(q)]
        Lambda_Gamma_new.putBlock(q, U_cut[repr(q)])#
    print Lambda_Gamma_new

    Gamma_Lambda_new = uni10.UniTensor([bdi_new] + bdo_old_list)
    Gamma_Lambda_new.setName('Gamma_Lambda_new')
    Gamma_Lambda_new.set_zero()
    for block in range(Gamma_Lambda_new.blockNum()):
        q = Gamma_Lambda_new.blockQnum(block)
        Gamma_Lambda_new.putBlock(q, VT_cut[repr(q)])#
    print Gamma_Lambda_new



    pass


iTEBD_update(Ham_Exp_dt, Lambda_B, Gamma_A, Lambda_A, Gamma_B, Lambda_B)

exit()

exit()

net = uni10.Network("LambdaInv_GammaLambda.net")
net.putTensor("Lambda_inv", Lambda_A_inv)
net.putTensor("Gamma_Lambda", Gamma_Lambda_new)
Lambda_A_new = net.launch()
print Y

exit()

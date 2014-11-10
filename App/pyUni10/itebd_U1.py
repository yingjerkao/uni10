#!/usr/bin/env python

# import modules
import math
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


# S=1/2 XXZ model
# q=2S
# S=+1/2 => q=+1
# S=-1/2 => q=-1

# set up quantum number for single-site operators
q0 = uni10.Qnum(0)
q1 = uni10.Qnum(1)
q_1 = uni10.Qnum(-1)
q2 = uni10.Qnum(2)
q_2 = uni10.Qnum(-2)


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


# iTEBD update
# todo
# def iTEBD_update(dt):
#     Ham_Exp_dt = Ham_Exp(Ham, dt)
#     print Ham_Exp_dt
#     pass


# setup dt and exp(-H*dt)
dt = 0.1
Ham_Exp_dt = Ham_Exp(Ham, dt)
Ham_Exp_dt.setName('Ham_Exp_dt')
print Ham_Exp_dt

# construct S via network
net = uni10.Network("t.net")
net.putTensor("Lambda_L", Lambda_B)
net.putTensor("Gamma_L", Gamma_A)
net.putTensor("Lambda_C", Lambda_A)
net.putTensor("Gamma_R", Gamma_B)
net.putTensor("Lambda_R", Lambda_B)
net.putTensor("Ham_Exp_dt", Ham_Exp_dt)
S = net.launch()
S.setName('S')
print S


# all eigenvalues, to be sorted
Lambda_all = []
# use dictionary to store the Lambda, U, VT for each block
Lambda = {}
U = {}
VT = {}

# SVD on all blocks
for block in range(S.blockNum()):
    q = S.blockQnum(block)
    M = S.getBlock(q)
    (U_q, Lambda_q, VT_q) = M.svd()
    U[repr(q)] = U_q
    VT[repr(q)] = VT_q
    Lambda[repr(q)] = Lambda_q
    for i in range(Lambda_q.col()):
        Lambda_all.append( Lambda_q[i] )

for block in range(S.blockNum()):
    q = S.blockQnum(block)
    print Lambda[repr(q)]

print 'Lambda=\n', Lambda
# print 'U=\n', U
# print 'VT=\n', VT
# print 'Lambda_all=\n', Lambda_all
Lambda_all.sort(reverse=True)
# print 'Lambda_all=\n', Lambda_all

# identify the smallest lambda to keep
D_cut = 4
if D_cut >= len(Lambda_all):
    Lambda_cut = Lambda_all[-1]
else:
    Lambda_cut = Lambda_all[D_cut-1]
print 'Lambda_cut=\n', Lambda_cut

# copy the in and out bond of the S separately
bdi_old_list=[]
bdo_old_list=[]
for b in range( S.bondNum() ):
    # print 'b=', b
    bond = S.bond(b)
    # print type(bond)
    if bond.type() == uni10.BD_IN :
        bdi_old_list.append( bond )
    else:
        bdo_old_list.append( bond )
# find the D for each q after truncation
# find Lambda_q for each block

# Lambda_cut
print '>>>>> Lambda_cut >>>>>'
L_cut = {}
for block in range(S.blockNum()):
    q = S.blockQnum(block)
    print q
    M = Lambda[repr(q)]
    print M
    lam = []
    for i in range(M.col()):
        print M[i], Lambda_cut
        if M[i] >= Lambda_cut:
            lam.append(M[i])
        else:
            break
    M_cut = uni10.Matrix(len(lam), len(lam), lam, True)
    print M_cut
    L_cut[repr(q)] = M_cut
print L_cut
exit()

Lambda_new = []
for n in range( S.blockNum() ):
    # print 'n= ', n, Lambda[n]
    Lambda_q = []
    for lam in Lambda[n]:
        if lam < Lambda_cut:
            print 'cut'
        else:
            Lambda_q.append(lam)
    # print Lambda_q
    Lambda_new.append( Lambda_q )

# print Lambda
print 'Lambda_new=\n', Lambda_new
for i in range(len(Lambda_new)):
    print len(Lambda_new[i])


# list of quantum numbers after the truncation
q_new_list = []
for n in range( S.blockNum() ):
    # print len( Lambda_new[n] )
    q = S.blockQnum( n )
    # print q
    for i in range( len(Lambda_new[n] ) ):
        q_new_list.append( S.blockQnum( n ) )

# new IN and OUT bound after the truncation
bdi_new = uni10.Bond(uni10.BD_IN, q_new_list)
bdo_new = uni10.Bond(uni10.BD_OUT, q_new_list)

# print bdi_old_list + bd_new_list
Gamma_Lambda_new = uni10.UniTensor(bdi_old_list + [bdo_new])
Gamma_Lambda_new.set_zero()
print Gamma_Lambda_new
for block in range(Gamma_Lambda_new.blockNum()):
    q = Gamma_Lambda_new.blockQnum(block)
    print 'block:', block, q
    print 'D_c', len(Lambda_new[block])
    M = Gamma_Lambda_new.getBlock(q)
    print M
    U = U_all[block]
    print U
    # Gamma_Lambda_new.putBlock(q, U)

# print Gamma_Lambda_new

# print Lambda_A_inv
exit()

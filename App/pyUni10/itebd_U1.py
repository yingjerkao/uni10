#!/usr/bin/env python

# import modules
import math
import copy
import pyUni10 as uni10

def debug_print(text, info):
    if debug == True:
        print ">>>>> debug >>>>>", text
        print info
debug = True
debug = False

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

########################################
# utility functions
# list of Qnums into list of bonds
def QL_2_BL(Q_list):
    q_list = []
    for Q in Q_list:
        for i in range(Q[1]):
            q_list.append(Q[0])
    return q_list

# init a rank-2 tensor
# Q_list = [(q0, 2), (q1, 3), (q_1, 3)]

# initialize a rank-2 tensor (1-in, 1-out) from a lsit of Qnums
def rank2tensor(Q_list):
    q_list = QL_2_BL(Q_list)
    bdi = uni10.Bond(uni10.BD_IN, q_list)
    bdo = uni10.Bond(uni10.BD_OUT, q_list)
    tensor = uni10.UniTensor([bdi, bdo])
    tensor.set_zero()
    return tensor

# T = rank2tensor(Q_list)
# debug_print("T", T)

# find inverse
# todo: do pseudo-inverse
def tensor_inv(tensor):
    tensor_inv = copy.copy(tensor)
    for block in range(tensor_inv.blockNum()):
        q = tensor_inv.blockQnum(block)
        M = tensor_inv.getBlock(q)
        # print M.col(), M.row()
        for i in range(M.col()):
            # print M[i, i]
            M[i, i] = 1.0/M[i, i]
        tensor_inv.putBlock(q, M)
    return tensor_inv

########################################
# iTEBD related functions
# set up Exp(-Ham*dt)
def Ham_Exp(Ham, dt):
    Ham_Exp_dt = copy.copy(Ham)

    QList = Ham_Exp_dt.getBlocks()
    for q in QList:
        debug_print('q=', q)
        M = Ham_Exp_dt.getBlock(q)
        (D, U) = M.eigh()
        UT=copy.copy(U)
        UT.transpose()
        for i in range(D.row()):
            D[i]=math.exp(-D[i] * dt)
        Ham_Exp_dt.putBlock(q, UT * D * U)
    return Ham_Exp_dt

# iTEBD update
def iTEBD_update(Ham_Expt_dt, Lambda_L, Gamma_L, Lambda_C, Gamma_R, Lambda_R):
    debug_print("iTEBD_update start", "")

    Lambda_L_Inv = tensor_inv(Lambda_L)
    Lambda_L_Inv.setName('Lambda_L_Inv')
    Lambda_R_Inv = tensor_inv(Lambda_L)
    Lambda_R_Inv.setName('Lambda_R_Inv')

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
    debug_print("S", S)

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
    debug_print("Lambda", Lambda)
    debug_print("U", U)
    debug_print("VT", VT)
    debug_print("Q_list", Q_list)

    Lambda_all.sort(reverse=True)
    debug_print("Lambda_all", Lambda_all)
    Lambda_all_M=uni10.Matrix(len(Lambda_all), len(Lambda_all), True)
    Lambda_all_M.setElem(Lambda_all)
    debug_print("Lambda_all_M", Lambda_all_M)
    M_norm = Lambda_all_M.norm()
    debug_print("Norm", M_norm)

    # identify the smallest lambda to keep
    if D_cut >= len(Lambda_all):
        Lambda_smallest = Lambda_all[-1]
    else:
        Lambda_smallest = Lambda_all[D_cut-1]

    debug_print("Lambda_smallest", Lambda_smallest)
    # print Lambda_all
    # print Lambda_smallest

    if D_cut >= len(Lambda_all):
        Lambda_all_M_cut=uni10.Matrix(len(Lambda_all), len(Lambda_all), True)
    else:
        Lambda_all_M_cut=uni10.Matrix(D_cut, D_cut, True)
    Lambda_all_M_cut.setElem(Lambda_all)
    debug_print("Lambda_all_M_cut", Lambda_all_M_cut)

    M_cut_norm = Lambda_all_M_cut.norm()
    debug_print("M_cut_norm", M_cut_norm)
    debug_print("M_norm-M_cut_norm", M_norm-M_cut_norm)

    # copy.copy the in and out bond of the S separately
    bdi_old_list = []
    bdo_old_list = []
    for block in range(S.bondNum()):
        bond = S.bond(block)
        if bond.type() == uni10.BD_IN :
            bdi_old_list.append( bond )
        else:
            bdo_old_list.append( bond )

    # perform the truncation
    debug_print("Perform the truncation", "")
    Q_list_cut = []
    Lambda_M_new = {}
    U_cut = {}
    VT_cut = {}

    for block in range(S.blockNum()):
        q = S.blockQnum(block)

        # M_q => M_q_cut
        M_q = Lambda[repr(q)]
        Lambda_q_cut = []
        for i in range(M_q.col()):
            if M_q[i] >= Lambda_smallest:
                Lambda_q_cut.append(M_q[i]/Lambda_all_M_cut.norm())
        debug_print("Lambda_q_cut", Lambda_q_cut)
        # init to zero??
        M_q_cut = uni10.Matrix(len(Lambda_q_cut), len(Lambda_q_cut), Lambda_q_cut, True)
        Lambda_M_new[repr(q)] = M_q_cut

        # U => U_cut
        U_q = U[repr(q)]
        U_q.transpose()
        # print U_q
        U_q_cut = uni10.Matrix(len(Lambda_q_cut), U_q.row(), U_q.getElem())
        U_q_cut.transpose()
        # print U_q_cut
        U_cut[repr(q)] = U_q_cut

        # VT => VT_cut
        VT_q = VT[repr(q)]
        VT_q_cut = uni10.Matrix(len(Lambda_q_cut), VT_q.col(), VT_q.getElem())
        # print VT_q
        # print VT_q_cut
        VT_cut[repr(q)] = VT_q_cut
        # Q
        Q = (q, len(Lambda_q_cut))
        Q_list_cut.append(Q)

    #
    # print U_cut
    # print VT_cut
    # print Q_list_cut
    # print Lambda_M_new

    Lambda_C_new = rank2tensor(Q_list_cut)
    Lambda_C_new.setName('Lambda_C_new')
    for block in range(Lambda_C_new.blockNum()):
        q = Lambda_C_new.blockQnum(block)
        Lambda_C_new.putBlock(q, Lambda_M_new[repr(q)])
    debug_print("Lambda_C_new", Lambda_C_new)

    # new IN and OUT bound after the truncation
    q_new_list = QL_2_BL(Q_list_cut)
    bdi_new = uni10.Bond(uni10.BD_IN, q_new_list)
    bdo_new = uni10.Bond(uni10.BD_OUT, q_new_list)

    Lambda_L_Gamma_L_new = uni10.UniTensor(bdi_old_list + [bdo_new])
    Lambda_L_Gamma_L_new.setName('Lambda_L_Gamma_L_new')
    Lambda_L_Gamma_L_new.set_zero()
    for block in range(Lambda_L_Gamma_L_new.blockNum()):
        q = Lambda_L_Gamma_L_new.blockQnum(block)
        Lambda_L_Gamma_L_new.putBlock(q, U_cut[repr(q)])
    debug_print("Lambda_L_Gamma_L_new", Lambda_L_Gamma_L_new)

    net = uni10.Network("LambdaInv_LambdaGamma.net")
    net.putTensor("Lambda_Inv", Lambda_L_Inv)
    net.putTensor("Lambda_Gamma", Lambda_L_Gamma_L_new)
    Gamma_L_new = net.launch()
    debug_print("Gamma_L_new", Gamma_L_new)

    Gamma_R_Lambda_R_new = uni10.UniTensor([bdi_new] + bdo_old_list)
    Gamma_R_Lambda_R_new.setName('Gamma_R_Lambda_R_new')
    Gamma_R_Lambda_R_new.set_zero()
    for block in range(Gamma_R_Lambda_R_new.blockNum()):
        q = Gamma_R_Lambda_R_new.blockQnum(block)
        Gamma_R_Lambda_R_new.putBlock(q, VT_cut[repr(q)])
    debug_print("Gamma_R_Lambda_R_new", Gamma_R_Lambda_R_new)

    net = uni10.Network("GammaLambda_LambdaInv.net")
    net.putTensor("Lambda_Inv", Lambda_R_Inv)
    net.putTensor("Gamma_Lambda", Gamma_R_Lambda_R_new)
    Gamma_R_new = net.launch()
    debug_print("Gamma_R_new", Gamma_R_new)

    debug_print("iTEBD_update end", "")

    return Gamma_L_new, Lambda_C_new, Gamma_R_new



# perform measurement
def measurement(Operator, Lambda_L, Gamma_L, Lambda_C, Gamma_R, Lambda_R):
    # transpose
    Lambda_L_t = copy.copy(Lambda_L)
    Lambda_L_t.transpose()
    Gamma_L_t = copy.copy(Gamma_L)
    Gamma_L_t.transpose()
    Lambda_C_t = copy.copy(Lambda_C)
    Lambda_C_t.transpose()
    Gamma_R_t = copy.copy(Gamma_R)
    Gamma_R_t.transpose()
    Lambda_R_t = copy.copy(Lambda_R)
    Lambda_R_t.transpose()

    # measurement network
    net = uni10.Network("measure.net")
    net.putTensor("Lambda_L", Lambda_L)
    net.putTensor("Gamma_L", Gamma_L)
    net.putTensor("Lambda_C", Lambda_C)
    net.putTensor("Gamma_R", Gamma_R)
    net.putTensor("Lambda_R", Lambda_R)
    net.putTensor("Operator", Operator)
    net.putTensor("Lambda_L_t", Lambda_L_t)
    net.putTensor("Gamma_L_t", Gamma_L_t)
    net.putTensor("Lambda_C_t", Lambda_C_t)
    net.putTensor("Gamma_R_t", Gamma_R_t)
    net.putTensor("Lambda_R_t", Lambda_R_t)

    result_T = net.launch()
    q = result_T.blockQnum(0)
    result_M = result_T.getBlock(q)
    result = result_M[0, 0]
    return result

########################################
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
debug_print("Sp_nosym", Sp_nosym)

Sp = uni10.UniTensor([bdi, bd_op1, bdo])
Sp.setName('Sp')
Sp.setRawElem(Sp_raw)
debug_print("Sp", Sp)

# Sm
Sm_nosym = uni10.UniTensor([bdi_nosym, bdo_nosym])
Sm_nosym.setName('Sm_nosym')
Sm_nosym.setRawElem(Sm_raw)
debug_print("Sm_nosym", Sm_nosym)

Sm = uni10.UniTensor([bdi, bd_op_1, bdo])
Sm.setName('Sm')
Sm.setRawElem(Sm_raw)
debug_print("Sm", Sm)

# Sz
Sz = uni10.UniTensor([bdi, bdo])
Sz.setName('Sz')
Sz.setRawElem(Sz_raw)
debug_print("Sz", Sz)

#  Id
Id = uni10.UniTensor([bdi, bdo])
Id.setName('Id')
Id.setRawElem(Id_raw)
debug_print("Id", Id)

# two-sites operators
# SpSm_nosym
SpSm_nosym = uni10.otimes(Sp_nosym, Sm_nosym)
SpSm_nosym.setName('SpSm_nosym')
debug_print("SpSm_nosym", SpSm_nosym)

# SmSp_nosym
SmSp_nosym = uni10.otimes(Sm_nosym, Sp_nosym)
SmSp_nosym.setName('SmSp_nosym')
debug_print("SmSp_nosym", SmSp_nosym)

# SpSm
SpSm = uni10.otimes(Sp, Sm)
SpSm.setName('SpSm')
SpSm.combineBond([2, 3, 4])
SpSm.setLabel([-1, -2, 1, 2])
debug_print("SpSm", SpSm)

# SmSp
SmSp = uni10.otimes(Sm, Sp)
SmSp.setName('SmSp')
SmSp.combineBond([2, 3, 4])
SmSp.setLabel([-1, -2, 1, 2])
debug_print("SmSp", SmSp)

debug_print("SpSm + SmSp", SpSm + SmSp)

# SzSz
SzSz = uni10.otimes(Sz, Sz)
SzSz.setName('SzSz')
SzSz.setLabel([-1, -2, 1, 2])
debug_print("SzSz", SzSz)

# IdId
IdId = uni10.otimes(Id, Id)
IdId.setName('IdId')
IdId.setLabel([-1, -2, 1, 2])
debug_print("IdId", IdId)

# set up Hamiltonian
# Ham = SpSm+SmSp+SzSz
Ham = 0.5 * (SpSm + SmSp) + SzSz
Ham.setName('Ham')
debug_print("Ham", Ham)


########################################
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

debug_print("Gamma A:", Gamma_A)
debug_print("Gamma B:", Gamma_B)

# Lambda
Lambda_A = uni10.UniTensor([bdi_v, bdo_v])
Lambda_B = uni10.UniTensor([bdi_v, bdo_v])

# uniform distribution of lambda
Lambda_A.identity()
Lambda_B.identity()

# fix the normalization
Lambda_A = math.sqrt(1.0/3.0) * Lambda_A
Lambda_B = math.sqrt(1.0/3.0) * Lambda_B
Lambda_A.setName('Lambda_A')
Lambda_B.setName('Lambda_B')

debug_print("Lambda_A", Lambda_A)
debug_print("Lambda_B", Lambda_B)


########################################
# main loop for imaginary evolution
print ">>>>> main loop starts >>>>>"
D_cut = 100

# run_sequence = [(0.1, 100), (0.25, 100)]
run_sequence = [(0.1, 100), (0.01, 100)]
for step in run_sequence:
    print step
    dt = step[0]
    n_iter = step[1]
    print dt, n_iter

    # setup dt and exp(-H*dt)
    # dt = 0.1
    Ham_Exp_dt = Ham_Exp(Ham, dt)
    Ham_Exp_dt.setName('Ham_Exp_dt')
    debug_print("Ham_Exp_dt", Ham_Exp_dt)

    for i in range(n_iter):
        (Gamma_A, Lambda_A, Gamma_B) = iTEBD_update(Ham_Exp_dt, Lambda_B, Gamma_A, Lambda_A, Gamma_B, Lambda_B)
        Gamma_A.setName('Gamma_A')
        Lambda_A.setName('Lambda_A')
        Gamma_B.setName('Gamma_B')
        # print Lambda_A
        # print Lambda_B

        (Gamma_B, Lambda_B, Gamma_A) = iTEBD_update(Ham_Exp_dt, Lambda_A, Gamma_B, Lambda_B, Gamma_A, Lambda_A)
        Gamma_B.setName('Gamma_B')
        Lambda_B.setName('Lambda_B')
        Gamma_A.setName('Gamma_A')
        # print Lambda_A
        # print Lambda_B

        debug_print("Lambda_A", Lambda_A)
        debug_print("Lambda_B", Lambda_B)

    E_a = measurement(Ham, Lambda_B, Gamma_A, Lambda_A, Gamma_B, Lambda_B)
    E_b = measurement(Ham, Lambda_A, Gamma_B, Lambda_B, Gamma_A, Lambda_A)
    Norm_a = measurement(IdId, Lambda_B, Gamma_A, Lambda_A, Gamma_B, Lambda_B)
    Norm_b = measurement(IdId, Lambda_A, Gamma_B, Lambda_B, Gamma_A, Lambda_A)
    print "Norm_a", Norm_a
    print "Norm_b", Norm_b
    print "E_a", E_a / Norm_a
    print "E_b", E_b / Norm_b
    E_exact = 1./4.-math.log(2.)
    print "E_exact", E_exact
    print "Err", (E_a+E_b)/2.0-E_exact

print ">>>>> main loop ends >>>>>"


print "HAPPY"
exit()

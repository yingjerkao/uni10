import sys
import pyUni10 as uni10

# U1 = 1, parity even.
q10 = uni10.Qnum(1, uni10.PRT_EVEN)
# U1 = -1, parity odd.
q_11 = uni10.Qnum(-1, uni10.PRT_ODD)

print "q10:", q10
print "q10:", q_11

print "q_11: U1 =", q_11.U1(), q_11.prt()

q_11.assign(-2, uni10.PRT_EVEN)
print "q_11(after assign):", q_11
# check the for fermionic
print "isFermioinc:", uni10.Qnum.isFermionic()

# Fermionic system
print "----- Fermionic -----"
# fermionic parity even, U1 = 1, parity even.
f0_q10 = uni10.QnumF(uni10.PRTF_EVEN, 1, uni10.PRT_EVEN);	# !!!
# fermionic parity odd, U1 = 1, parity even.
f1_q10 = uni10.QnumF(uni10.PRTF_ODD, 1, uni10.PRT_EVEN);	# !!!

print "f0_q10:", f0_q10
print "f1_q10:", f1_q10
print "f1_q10: fermionic parity =", f1_q10.prtF()
print "isFermioinc:", uni10.Qnum.isFermionic()



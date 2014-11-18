import sys
import pyUni10 as uni10

q1 = uni10.Qnum(1)
q0 = uni10.Qnum(0)
q_1 = -q1

#Constrcut Bond with Qnum array
bd = uni10.Bond(uni10.BD_IN, [q1, q1, q0, q0, q0, q_1]);
# Print out a Bond
print "Bond bd:",  bd

#---------
qnums = bd.Qlist();
dges = bd.degeneracy();
for key in dges:
	print key, dges[key]

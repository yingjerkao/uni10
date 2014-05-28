import pyUni10 as uni10

elem = uni10.double_arr([1.0/4,  0,      0,     0,\
				         0, -1.0/4,  1.0/2,     0,\
				         0,  1.0/2, -1.0/4,     0,\
				         0,      0,      0, 1.0/4])

m = uni10.Matrix(4, 5)
m.set_zero();
m.orthoRand();
for e in elem:
	print e

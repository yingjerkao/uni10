import sys
sys.path.append('..')
import pyUni10 as uni10

elem = [1.0/4,      0,      0,     0,\
				    0, -1.0/4,  1.0/2,     0,\
				    0,  1.0/2, -1.0/4,     0,\
				    0,      0,      0, 1.0/4]

H = uni10.Matrix(4, 4, elem)

rets = H.svd()
for ret in rets:
	print ret


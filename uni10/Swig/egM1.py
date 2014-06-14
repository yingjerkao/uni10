import pyUni10 as uni10

elem = uni10.double_arr([\
	1.0/4,  0,      0,     0,\
	0, -1.0/4,  1.0/2,     0,\
	0,  1.0/2, -1.0/4,     0,\
	0,      0,      0, 1.0/4])

m = uni10.Matrix(4, 4);
m.randomize();
m2 = uni10.Matrix(5, 5);
m2.randomize();
sub = uni10.Matrix(3, 5, m2.elem());

print m
print m2.trace()
print sub * 10

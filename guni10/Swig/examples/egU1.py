import sys
import pyUni10 as uni10

# Construct spin 1 Heisenberg model
# Raw element
heisenberg_s1 = [\
		 1, 0, 0, 0, 0, 0, 0, 0, 0,\
		 0, 0, 0, 1, 0, 0, 0, 0, 0,\
		 0, 0,-1, 0, 1, 0, 0, 0, 0,\
		 0, 1, 0, 0, 0, 0, 0, 0, 0,\
		 0, 0, 1, 0, 0, 0, 1, 0, 0,\
		 0, 0, 0, 0, 0, 0, 0, 1, 0,\
		 0, 0, 0, 0, 1, 0,-1, 0, 0,\
		 0, 0, 0, 0, 0, 1, 0, 0, 0,\
		 0, 0, 0, 0, 0, 0, 0, 0, 1]

# Create in-coming and out-going bonds, without any symmetry.
bdi = uni10.Bond(uni10.BD_IN, 3 );
bdo = uni10.Bond(uni10.BD_OUT, 3 );

# Create tensor from the bonds and name it "H".
H = uni10.UniTensor([bdi, bdi, bdo, bdo]);
H.setRawElem(heisenberg_s1);
print H

# Since it has U1 symmetry(total Sz conserved)
# Add U1 quantum number to the states of bonds.
q0 = uni10.Qnum(0);
q1 = uni10.Qnum(1);
q_1 = uni10.Qnum(-1);
# Create in-coming and out-going bonds
bdi = uni10.Bond(uni10.BD_IN, [q1, q0, q_1] );
bdo = uni10.Bond(uni10.BD_OUT, [q1, q0, q_1] );

# Create tensor from the in-coming and out-going bonds.
H_U1 = uni10.UniTensor([bdi, bdi, bdo, bdo], "H_U1");
H_U1.setRawElem(heisenberg_s1);
print H_U1

# See the details of the blocks
print "The number of the blocks =", H_U1.blockNum(), "\n"
blocks = H_U1.getBlocks();
for q in blocks:
	print q, blocks[q]

H_U1.save("egU1_H_U1");

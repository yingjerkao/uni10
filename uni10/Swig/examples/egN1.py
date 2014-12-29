import sys
import pyUni10 as uni10

# Read in the tensor H_U1 which is written out in example egU1 and W, WT in example egU3
H_U1 = uni10.UniTensor("egU1_H_U1");
W = uni10.UniTensor("egU3_W");
WT = uni10.UniTensor("egU3_WT");

# Create network by reading in network file "egN1_network"
net = uni10.Network("egN1_network");
# Put tensors to the Network net
net.putTensor("H", H_U1);
net.putTensor("W", W);
net.putTensor("WT", WT);

# Perform contractions inside the tensor network
print net.launch();
# Print out the network
print net;

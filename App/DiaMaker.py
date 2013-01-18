 #!/usr/bin/python
import sys
import os

#Ordering policy:
#1. The contraction of the two tensors has a resulting tensor with least bonds.
#2. If the resulting tensor has the same number of bonds, choose the one that the 
#   original two tensors have the maximun number of bonds.
#3. If two pair of tensors, which are to be contracted, are on par with respect to 
#   the two conditions above, by default is contracting from the begining of the label list.

#labels is a list of lable
def ordering(inlabels):
	labels = []
	for num in inlabels:
		labels.append(num)
	order = []
	while len(labels) > 1:
		smallest_out_bonds = 100
		largest_total_bonds = -1
		pair = [-1, -1]
		for i in range(len(labels) - 1):
			for j in range(i+1, len(labels)):
				total_bonds = len(labels[i]) + len(labels[j])
				#find intersection
				out_bonds = total_bonds - 2 * len([val for val in labels[i] if val in labels[j]])
				flag = False
				if total_bonds > out_bonds:
					if 	out_bonds < smallest_out_bonds:
						smallest_out_bonds = out_bonds
						largest_total_bonds = total_bonds
						flag = True
					elif out_bonds == smallest_out_bonds:
						if total_bonds > largest_total_bonds:
							smallest_out_bonds = out_bonds
							largest_total_bonds = total_bonds
							flag = True
				if flag:
					pair[0] = i
					pair[1] = j
		#Some error
		if pair == [-1, -1]:
			return [], []
		tar1 = labels.pop(pair[0])
		tar2 = labels.pop(pair[1] - 1)
		intersect = [val for val in tar1 if val in tar2]
		new_label = [val for val in tar1 + tar2 if val not in intersect]
		labels.append(new_label)
		order += pair
	return order, labels[0]


if len(sys.argv) < 2:
	print "Give a diagram file"
	sys.exit(0)
f = open(sys.argv[1])
lines = f.readlines()
f.close()
tensors = []
labels = []
i = 0
Did = ''
Tout = False
func = ''
while i < len(lines):
	line = lines[i]
	if line.find('Function') >= 0 and line.find(':') >= 0:
		func = line.split(':')[-1].strip()
	elif line.find('Did') >= 0 and line.find(':') >= 0:
		Did = 'DIAG_' + line.split(':')[-1].strip()
	elif line.find('Tensors') >= 0 and line.find(':') >= 0:
		i += 1
		while i < len(lines):
			tmp = lines[i].split(':')
			tensors.append(tmp[0])
			label = tmp[-1].split()
			label = [int(l.strip()) for l in label]
			labels.append(label)
			i += 1
			if lines[i].find('Tout') >= 0:
				Tout = True
				break
		if not Tout:
			print "No 'Tout' in diagram file, Error!"
			sys.exit(0)
		outLabel = [int(l.strip()) for l in  lines[i].split(':')[-1].split()]
	i += 1
if not func:
	print "No 'Function Name' in diagram file, Error!"
	sys.exit(0)
if not Did:
	print "No 'Did' in diagram file, Error!"
	sys.exit(0)
if len(tensors) == 0:
	print "No tensor in diagram file, Error!"
	sys.exit(0)

od, outL = ordering(labels)
reshape = []
for l in outL:
	for i in range(len(outLabel)):
		if l == outLabel[i]:
			reshape.append(i)
			break
#sanity check
if len(outLabel) != len(outL):
	print 'There is some problem in your diagram, output bond number is incorrect!!'
	sys.exit(0)
for l in outLabel:
	if l not in outL:
		print '"%d" is a incorrect output bond label!! Check your "Tout"' % l
		sys.exit(0)
for l in outL:
	if l not in outLabel:
		print '"%d" is not in "Tout, Error!!!"' % l
		sys.exit(0)
for i in range(len(outL)):
	if i not in reshape:
		print 'Error in your labeling!!!'
		sys.exit(0)

myDiasfn = 'myDias.c'
if os.access(myDiasfn, os.F_OK):
	f = open(myDiasfn, 'r')
	lines = f.readlines()
	for line in lines:
		if 'void' in line:
			tmp = line.split()
			if tmp[1].split('(')[0].strip() == func:
				print 'Function: "%s" has already existed!!! Error' % (func)
				sys.exit(0)
else:
	f = open(myDiasfn, 'w')
	f.write('/*Function API*/\n')
	f.write('/*Function Definition*/\n')
	f.close()

#Add Did
didfn = 'network.h'
num = 0
dids = []
if os.access(didfn, os.F_OK):
	f = open(didfn, 'r')
	lines = f.readlines()
	f.close()
	flag = True
	begin = False
	for line in lines:
		if line.find('/*Did definition*/') >= 0:
			begin = True
		if line.find('#define') >= 0 and begin:
			if line.find(Did) >= 0:
				flag = False
			dids.append(line.split()[-2].strip())
			num += 1
	#check diasNet
	inusefn = 'diasNet'
	count = 0
	correct = True
	if os.access(inusefn, os.F_OK):
		f = open(inusefn)
		lines = f.readlines()
		for line in lines:
			if line.find('Did:') >= 0: 
				if len(dids) > count and dids[count] != line.split(':')[-1].strip():
					correct = False
					break
				count += 1
		if len(dids) < count:
			correct = False
			print 'The number of diagrams in "diasNet" is more than the number of did defined in "network.h"'
		if len(dids) > count:
			correct = False
			print 'There\'s some diagrams in "diasNet" that is not defined in "network.h"'
	elif len(dids) > 0:
		correct = False	
	if not correct:
		print 'Fatal Error!!!  Diagrams\' Dids in network.h is not match for that in "diasNet" file'
		print 'Please reconstruct the network.h.'
		ans = raw_input('Do you want to remove the incorrect "network.h" and "diasNet" (Y or N): ')
		if ans == 'Y':
			os.remove('diasNet')
			os.remove('network.h')
		sys.exit(0)
	#end of checking diasNet
	if not begin:
		f = open(didfn, 'r')
		lines = f.readlines()
		f.close()
		f = open(didfn, 'w')
		lines.insert(0, '/*Did definition*/\n')
		lines.insert(1, '#define %s %d\n' %(Did, num))
		lines.insert(2, '/*End Did definition*/\n')
		lines.insert(0, '#define DIAG_NUM %d\n' % (num + 1))
		lines.append('#include "myDias.c"\n')
		f.write(''.join(lines))
		f.close()
	elif flag:
		f = open(didfn, 'r')
		lines = f.readlines()
		f.close()
		f = open(didfn, 'w')
		lines = lines[lines.index('/*Did definition*/\n') :]
		lines.insert(0, '#define DIAG_NUM %d\n' % (num + 1))
		lines.insert(lines.index('/*End Did definition*/\n'), '#define %s %d\n' %(Did, num))
		if '#include "myDias.c"\n' not in lines:
			lines.append('#include "myDias.c"\n')
		f.write(''.join(lines))
		f.close()
	else:
		print 'Did: "%s" has existed!!!' % Did
		sys.exit(0)
else:
	f = open(didfn, 'w')
	f.write('#define DIAG_NUM 1\n')
	f.write('/*Did definition*/\n')
	f.write('#define %s %d\n' %(Did, num))
	f.write('/*End Did definition*/\n')
	f.write('#ifdef GPU\n')
	f.write('# include "Tensor.cu"\n')
	f.write('#else\n')
	f.write('# include "Tensor.cpp"\n')
	f.write('#endif\n')
	f.write('#include "netBasic.cpp"\n')
	f.write('#include "myDias.c"\n')
	f.close()

#Add diagram to diasNet	
inusefn = 'diasNet'
f = open(inusefn, 'a')
f.write('Did: %s\n'% (Did))
f.write('%d\n' %(len(tensors)))
for label in labels:
	f.write('%d ' % len(label))
f.write('%d ' % (len(outLabel)))
f.write('\n')
for label in labels:
	for i in label:
		f.write('%d ' % (i))
	f.write('\n')
for r in outLabel:
	f.write('%d ' % (r))
f.write('\n')		
for o in od:
	f.write('%d ' % (o))
f.write('\n')		
f.close()

#Create a function in myDias.c
myDiasfn = 'myDias.c'
api = 'void %s(' %(func)
for ten in tensors:
	api += 'Tensor* %s, '%ten
api += 'Tensor* Tout' + ');\n'
cnt = api[:-2] + '{\n' + "\tvector<Tensor*>Tlist;\n"
for ten in tensors:
	cnt += "\tTlist.push_back(%s);\n" % (ten)
cnt += "\toperate(%s, Tlist, Tout);\n}\n" % (Did)
f = open(myDiasfn, 'r')
lines = f.readlines()
f.close()
f = open(myDiasfn, 'w')
lines.insert(lines.index('/*Function API*/\n') + 1, api)
lines.append(cnt)
f.write(''.join(lines))
f.close()

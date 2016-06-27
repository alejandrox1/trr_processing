#!/usr/bin/python

import sys

file_handle = str(sys.argv[1])
new_top = str(sys.argv[2])
itp = str(sys.argv[3]) 

def top_line(file_handle, itp, new_top):
	f = open(file_handle, 'r')
	f.readline()
	x1 = []
	x2 = []	

	for num,line in enumerate(f):
		if "moleculetype" in line:
			x1.append(num)
		if "Include" in line:
			#if num > x1[0] :
			x2.append(num)
	f.close()	
	return (x1[0],x2[2])

def top_cut(file_handle, itp, new_top):
	f = open(file_handle, 'r')
        f.readline()

	t = open(new_top, 'a')
	n = open(itp, 'a')
	x1, x2 = top_line(file_handle, itp, new_top)
	
	for num,line in enumerate(f):
		if ( num <= x1 ):
                        t.write(line)
		if ( num > x1 ) and ( num < x2 ):
			n.write(line)
		if ( num >= x2 ):
			t.write(line)

if __name__=="__main__":
	
	top_cut(file_handle,itp, new_top)
		

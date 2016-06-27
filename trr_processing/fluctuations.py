#!/usr/bin/python

import matplotlib as mplt
import matplotlib.pyplot as plt
import sys
import numpy
from matplotlib.backends.backend_pdf import PdfPages
from random import randrange
import math 

file_name = str(sys.argv[1])

def tofloat(l):
        l = [float(x) for x in l]
        return l


def plotxvg(filename):
	datax=[]
	datay=[]
	infile=open(filename, 'r')
	for line in infile:
		if (line.find('#')==-1) and (line.find('@')==-1) and (line.find('&')==-1):
			if (line!='0'): 
				x,y = line.split()
				x = float(x); y = float(y)
				datax.append(x)
				datay.append(y)
	infile.close()
	
	#datax = tofloat(datax)
        #datay = tofloat(datay)

	return datax, datay

def fluctuations( n , out_handle ):
	datax, datay = plotxvg(file_name)

	energy_aver = sum(datay)/len(datay)	

	e_square = []      ######## LIST        
        for i in datay:
        	e_square.append( (i)**2 )

	# Fluctuatiuon in n ns blocks
	num = n * 1000
        interval = len(datay)/num
        blocks = []
        block_sums = []
        for i in range( len(datay)/interval ):
                #block_sums = []
                block_sums.append( sum(e_square[(interval * i):(interval * (i + 1))]) / len(e_square[(interval * i):(interval * (i + 1))]) )
        for j in block_sums:
                blocks.append( j - energy_aver**2 )
        
	px = range(len(datay)/interval)
        pdf_pages = PdfPages( out_handle )
        plt.rc('text', usetex=False)
        fig = plt.figure(figsize=(10,8))
        plt.plot(px,blocks, '-o')
        plt.title("Fluctuations")
#       plt.xlabel("")
        plt.ylabel("(Delta E)^2")
        pdf_pages.savefig(fig)
        pdf_pages.close()

	return out_handle
		

if __name__=='__main__':
#	plotxvg(file_name)
	
	# Average Energy	
#	energy_aver = sum(datay)/len(datay)   ##### VALUE
#       print "average energy",energy_aver
#	print "average energy squared",energy_aver**2
#	e_square = []      ######## LIST	
#	for i in datay:
#		e_square.append( (i)**2 )
#	energy_square_aver = sum(e_square)/len(e_square)
#	print "average of the square energy",energy_square_aver
#	print "delta square", (energy_square_aver - energy_aver**2)
#	k = 1.3806488e-23
#	na = 6.022141e23
#	c_non = (energy_square_aver - energy_aver**2) / ( (300**2) * k)
#	cv = (c_non/na) * 1000
#	print "The heat capacity is",cv,"KJ/K"
#	print

	fluctuations( 1, 'py_fluctuations.pdf' ) 	

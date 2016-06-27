#!/usr/bin/python

import os
import glob
import time


# current slurm output file
newest = max( glob.iglob('*.[Oo][Uu][Tt]'), key=os.path.getctime )


# Date (d:m:y) 
q = time.strftime("%d:%m:%y")
q = q.split(':')
queue = [ int(i) for i in q ]


# Time (H:M:S)
t = time.strftime("%H:%M:%S")
t = t.split(':')
time = [ int(i) for i in t ]


def slurm_parser(newest):
	if newest == "squeue.out":
            slurm = open("squeue.out", 'r')
            info_all = slurm.readlines()
            info = []
            for line in info_all:
                line = line.split()
                info.append(line)

            if ( info[-1][4] != 'R' ):         
                day_str = time.strftime("%d")
                day= int(day_str)
                info=[day, '1:00:00', 2014]

        else:
            slurm = open(newest, 'r')
	    info_all = slurm.readlines()
	    info = []
	    for line in info_all:
	            line = line.split()
		    info.append(line)
	
        # Return the expected finish time (line before last)
	return info[-3],info[-2] # day, hh:mm:ss



def format(newest, queue, time):
	# Set up Slurm parser output in proper format
	day, hour = slurm_parser(newest)
	hour = hour.split(':')
        day = int(day)

	# Day and time differences
	d_diff = abs(day - queue[0])
	t_diff = [ abs( hour[i] - time[i] ) for i in range(len(time)) ]


def day_secs(d_diff):
	return ( d_diff *24 * 60 *60 )


def time_insecs(t_diff):
	s = 0				# Counter
	s += t_diff[0] * 60 * 60	# Hours to seconds
	s += t_diff[1] * 60		# Minutes to seconds
	s += t_diff[2] 			
	return s


if __name__=="__main":

	s1 = day_secs(d_diff); s2 = time_insecs(t_diff)
	seconds = s1 + s2i + 600
	time.sleep(seconds)	# Sleep until simulation finishes


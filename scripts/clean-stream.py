#!/usr/bin/python
# coding=utf-8

# clean-stream.py
#
# Remove non-indexed frames from a stream
#
# Copyright © 2013 Fedor Chervinskii
#
# Authors:
#   2013 Fedor Chervinskii <fedor.chervinskii@gmail.com>
#
# This file is part of CrystFEL.
#
# CrystFEL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CrystFEL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
from itertools import islice
import sys
import re


infile_1 = open (sys.argv[1],'r')
infile_2 = open (sys.argv[1],'r')
outfile = open (sys.argv[2],'w')
Nfile = open ('N.dat','w')

suited = False
counter1 = -1
counter2 = 0
start_index = 0
n_patt = 0
num_peaks = 0
indexed = False
num_suited = 0

line = infile_1.readline()

reg0 = re.compile('----- Begin chunk')
reg1 = re.compile('----- End chunk')
reg6 = re.compile('^indexed_by')
reg7 = re.compile('none')

while  (line != ''):

	counter1 += 1

	if (reg6.match(line)):
	        if reg7.search(line) :
	          indexed = False
	        else:
	          indexed = True

	if (reg0.match(line)):
		if (start_index == 0) :
			while (counter2 < counter1) :
				outline = infile_2.readline()
				outfile.write(outline)
				counter2 += 1
		start_index = counter1

	if (reg1.match(line)):
		n_patt += 1
		if indexed :
			suited = True		
		if suited :
			num_suited += 1 
		while (counter2 <= counter1) :
			outline = infile_2.readline()
			if suited :	
				outfile.write(outline)
			counter2 += 1
		suited = False 

	line = infile_1.readline()	
	
print '%d suited of %d patterns have been extracted and saved as %s' % (num_suited, n_patt, sys.argv[2])
Nfile.write('%d' % num_suited)

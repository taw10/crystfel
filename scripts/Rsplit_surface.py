#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Plot Rsplit as a contour map
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
from array import array
from numpy import *

import matplotlib

"""there could be problems with dependencies, in this case, """
#matplotlib.use('PS')

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data_name = sys.argv[1]

Nfile = open('N.dat', 'r')
line = array(Nfile.readlines())[0]
N0 = eval(line.split()[0])
steps = arange(1, 10)

verts = []

for i in steps :

	filename='rsplit'+str(i)+'.dat'

	file = open(filename, 'r')

	lines = array(file.readlines())
	N = lines.size

	for j in range(1, N) :
		numbers = lines[j].split()
		if 0 < eval(numbers[1]) <= 100 : verts.append((eval(numbers[0]),i,eval(numbers[1])))

x, y, z = zip(*verts)

fig = plt.figure()

ax = host_subplot(111, axes_class=AA.Axes)

X = arange(1, 5, 0.05)
Y = arange(0, 10, 0.1)
X, Y = meshgrid(X, Y)
Z = griddata(x,y,z,X,Y)

cm = ax.pcolormesh(X, Y, Z, cmap = plt.get_cmap('spectral'), vmin=0, vmax=100)

plt.title(r"$R_{split}$ surface for %s" % data_name, fontsize=16)

ax.set_xlabel(r"Resolution (1/d) / $nm^{-1}$", fontsize=16)
ax.set_ylabel(r"Number of patterns", fontsize=16)
ypoints = [50000,40000,30000,20000,10000,9000,8000,7000,6000,5000,4000,3000,2000,1000,900,800,700,600,500,400,300,200,100]
temp = []
for ypoint in ypoints : temp.append(N0/ypoint)
ax.set_yticks(log2(temp)+1)
labels = []
for ypoint in ypoints :
	if ypoint in [50000,10000,5000,1000,500] : labels.append(str(ypoint))
	else : labels.append("")

ax.set_yticklabels(labels)

ax.set_ylim(log2(temp[0])+1,log2(temp[-1])+1)

cb = fig.colorbar(cm, shrink=0.5, aspect=10)
cb.set_label("Rsplit / \%")

plt.savefig(data_name[:-7] + '.png')
#plt.savefig(data_name[:-7] + '.ps')
plt.show()


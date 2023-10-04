#!/usr/bin/env python

import os
import numpy as np

bias = np.loadtxt("COLVAR",usecols=3)
bmax = max(bias)
kbt = 2.494339

colvar = np.loadtxt("COLVAR")
dw = []
for line in colvar:
    a = line[1], line[2], np.exp((line[3]-bmax)/kbt)
    dw.append(a)

np.savetxt('dist.s.weight',dw,fmt='%1.9f')

# calculate the free-energy landscape and the associated errors for different block-sizes (i) 
# then calculate the average error in each bin
# python3 do_block_fes.py CV-bias-file No-CVs Min-val Max-val No-bins Kbt block-size

errors = []

for i in range(1,1000,10):
    os.system('python3 do_block_fes.py dist.s.weight 2 0.3 1.5 51 -2 1 51 2.494339 {i}'.format(i=i))
    os.system('grep -v "Infinity" fes.{i}.dat > fes.dat'.format(i=i))
    os.system('sed -i "/^$/d" fes.dat')
    file = "fes.dat"
    err = np.genfromtxt(file,usecols=3)
    average = i, sum(err)/len(err)
    errors.append(average)

np.savetxt('errors.block',errors,fmt='%1.9f')


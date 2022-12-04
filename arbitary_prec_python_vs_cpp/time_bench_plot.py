# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 00:16:35 2022

@author: elekp
"""

from sys import argv
import numpy as np
import matplotlib.pyplot as plt

num=0
if (len(argv)>=2):
    num=argv[1]
else:
    print("nincsenek fileok")

mat_size=np.array([])
time1=np.array([])
time2=np.array([])
time3=np.array([])

for i in range(1,num):
    name="tmp/"+str(i)+".txt"
    file=open(name,"r")
    inputs0=file.readline()
    inputs1=file.readline()
    inputs2=file.readline()
    inputs3=file.readline()
    mat_size=np.append(mat_size,int(inputs0))
    time1=np.append(time1,float(inputs1))
    time2=np.append(time2,float(inputs2))
    time3=np.append(time3,float(inputs3))


fig, ax = plt.subplots()
ax.set_title('run time')

l1=ax.scatter(mat_size,time1, color='b')
l2=ax.scatter(mat_size,time2)
l3=ax.scatter(mat_size,time3)

#ax.legend([l1, l3], ['Python', 'C++ -O3 -flto'])
ax.legend([l1, l2, l3], ['Python', 'g++ -O1' , 'g++ -O3 -flto'])


#ax.set_yscale('log')
ax.set_xlabel('num of points (p)')
ax.set_ylabel('Run time (s)')

fig.savefig("kep.svg")

"""
fig1, ax1 = plt.subplots()
ax1.set_title('python/C run time')
ax1.set_xlabel('Mátrix mérete (nxn-es mátrix esetén n)')

#ax1.scatter(mat_size,time_py/time_c)
#ax1.scatter(mat_size,time_numpy_py/time_c)
ax1.scatter(mat_size,time1/time2)
#ax1.set_yscale('log')

fig2, ax2 = plt.subplots()
ax2.set_title('python/python+numpy run time')
ax2.set_xlabel('Mátrix mérete (nxn-es mátrix esetén n)')

ax2.scatter(mat_size,time_py/time_numpy_py)
#ax1.set_yscale('log')


#"""

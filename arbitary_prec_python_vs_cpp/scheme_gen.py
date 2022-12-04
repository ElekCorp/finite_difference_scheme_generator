# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 00:40:05 2022

@author: elekp
"""

from sympy import Matrix, zeros, factorial, Integer
from sympy import init_printing, pprint
init_printing(use_unicode=True, wrap_line=False)

from sys import argv

def matrix_coef(i,j,p,k):
    return Integer(j-k)**i
    #return (j-k)**i
    

def init_mat(p,k):
    num=2*p+1
    A=zeros(num,num)
    for i in range(num):
        for j in range(num):
            A[i,j]=matrix_coef(i,j,p,k)
    return A

def init_vec(m,p):#m-th derivative
    num=2*p+1

    b=zeros(num,1)
    b[m,0]=factorial(m)
    return b


def get_coeff_vector(p,m,k):
    A=init_mat(p,k)
    b=init_vec(m,p)
    
    x=A.LUsolve(b)
    return x

p = 1
if(len(argv) >= 2):
    p = int(argv[1])
if(p <= 0):
    p = 1

k=p
m=2

x=get_coeff_vector(p,m,k)
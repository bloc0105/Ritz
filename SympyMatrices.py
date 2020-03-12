'''
Created on Feb 17, 2020

@author: Andy Block
'''
import sympy as sym
import numpy as num


a = num.linspace(0,1,5,endpoint=True)

# print(a)

q =  num.meshgrid(a,a)

print(q)
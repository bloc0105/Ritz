import numpy
from matplotlib import pyplot as plotting

'''
Created on Dec 19, 2019

@author: andy
'''
x = numpy.linspace(0,10,10)
print(x)
y = [zerp**2 for zerp in x]

plotting.plot(x,y)
plotting.ylabel("Stuff")
plotting.show()
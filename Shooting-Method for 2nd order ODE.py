# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 00:10:47 2021

@author:CoderNova-asr

"""

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
#the eq to be solves by shooting:u"-(1-x/5)u=x
print('genral form of differential eq is taken as:\n')
print("y'' =p(x)y' +q(x)y +f(x)")
print('\n.........\n')
X=[]
Y=[]
def p(x):
    return 0    #p(x)
def q(x):
    q=1-x/5     #q(x)
    return q
def f(x):
    return x    #f(x)

xi=float(input('enter the initial vlaue of x'))
xf=float(input('enter the final vlaue of x'))
yi=float(input('enter the initial vlaue of y'))
yf=float(input('enter the final vlaue of y'))




"""
taking z=dy/dx we gwt two differential eq.s

dy/dx=z
dz/dx=pz+qy+x

Hence

[y'][z']=[y][z]*[0,1][p,(qy+x)/z]
we take

i=[y][z]
hence

           [y'][z']= i * [0,1][p,(qy+x)/z]
"""

#making the differential equation matrix to use
eq=lambda x,i:np.dot(np.array([[0,1],[p(x),(q(x)*i[0]+x)/i[1]]]),i)

x=np.linspace(xi,xf,10)
#a function to solve IVPs and return an array of solns for y(x)

def sol(m):
    sol= solve_ivp(eq, [xi,xf],[yi, m],t_eval= x)
    y=sol.y[0]
    return y
def plot(n):
    soln=solve_ivp(eq, [xi,xf],[yi, n], t_eval = x)
    return soln
#guesses
g1=-1.5
g2=-3

#solutions for the guesses
y1,y2=sol(g1),sol(g2)

#interpolated guess
g3= g2+ ((yf)-y2[-1])*(g1-g2)/(y1[-1]-y2[-1])

#output
print("\n**********\n")
print("Initial Guess: y(x) values-\n",sol(g1))

print("Secondary Guess: y(x) values-\n",sol(g2))

print("Interpolated value: y(x)s\n",sol(g3))

print("The interpolated y'(x)",g3)

#plotting
plt.plot(x,y1,label='First guess')    #plotting 1st guessed u'(x)
plt.plot(x,y2, label='Second guess')  #plotting 2nd guessed u'(x)
plt.plot(x,sol(g3),label= 'Sol by interpolation') #" interpolated "
plt.plot()
plt.xlabel('x')                                #labelling    
plt.ylabel('y(x)')
plt.legend(loc='best')


# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 13:24:55 2021

@author:CoderNova-asr
"""
#importing required packages
import matplotlib.pyplot as plt               #plotting
import numpy as np            #in-built mathematical functions 


# solve 2nd order boundary value problem using finite difference
""""Formula for finding y using finite difference method:
(2+p*h)*y(i-1)-(2*q*h^2+4)*y(i)+(2-p*h)*y(i+1)=2*(h^2)*F(x)"""

print("The general expression of a 2nd order differential eq.")
print("\ny''=p(x)*y'+q(x)*y+F(x)")
p=float(input('enter p(x)'))
q=float(input('enter q(x)'))

n=int(input('enter the no. of steps to be formed'))

#boundary conditions
a=float(input('enter a'))
b=float(input('enter b'))
h=(b-a)/n
print('stepsize found',h)
y1=float(input('Enter the value of y at a '))
y2=float(input('Enter the value of y at b '))
#coefficients of y
A=-2*q*h**2-4         
B=h*p+2
C=-h*p+2
x=a                        #initial limit of x
    
X=[]
Y=[]
Z=[]    
for i in range(n-1):
    z=C
    Z.append(C)
    
for i in range(n):
    z=A
    Y.append(z)
    
for i in range(n+1):
   z=B
   X.append(z)

#print(X,Y,Z)  
#print(np.diag(X,0),np.diag(Y,1),np.diag(Z,2)) 
#making the tri diagonal matrix
M=np.diag(X,0)+np.diag(Y,1)+np.diag(Z,2)
print(M)

#deleting unnecessary rows
M=np.delete(M,(n-1,n),axis=0) #axis=0 signifies rows are to be deleted
M=np.delete(M,(0,n),axis=1)#axis=1 for columns

#the required square matrix for the solution to be found  
print('tri-diagonal matrix consisting of the coefficients of y\n',M)


#making arrays for thomas method
A=[]
B=[]
C=[]
for i in range (n-2):                 #upper and lower diagonals
    a=M[i+1][i]
    A.append(a)
    c=M[i][i+1]
    C.append(c)
for i in range(n-1):
    b=M[i][i]
    B.append(b)                          #central diagonal

print('\n........\n')          
print(A,B,C)                    #matrices formed using the elements
#finding values of y at the discrete values of x

f=x*(x-4)-(h*p+2)*y1 #F(X) in the given eq
X=[x]
D=[f]

for i in range(n-1):
    x=x+h
    if i==n-3:
        f=x*(x-4)-2*y2+h*p*y2
        D.append(f)
    else:
        f=x*(x-4)
        D.append(2*f*h**2)
    X.append(x)
    
D=np.delete(D,(0,n))       #values not required to find the values of y

print('\n........\n')
print('constant matrix',D)  
W=np.linalg.solve(M,D)

#finding the values of y


"""thomas algorithm"""
D[0]=D[0]/B[0]
C[0]=C[0]/B[0]
for i in range(1,n-2):
     C[i]=C[i]/(B[i]-A[i-1]*C[i-1])
for i in range(1,n-1):
     D[i]=(D[i]-A[i-1]*D[i-1])/(B[i]-A[i-1]*C[i-1])
     
y=D[n-2]
Yin=[]
Y=[y1]
Yin.append(y)
for i in range(0,n-2):
    y=D[n-3-i]-y*C[n-3-i]
    Yin.append(y)
#print(Yin)                  matrix is formed with elements in inverse order
for i in range(n-1):
    y=Yin[n-i-2]
    Y.append(y)


Y.append(y2)    

print('\n........\n')
print('values of x at which y is found',X)
print('\n........\n')
print('discrete values of y found using thomas algorithm')
print(Y)
print('\n........\n')

print('solution using in-built function',W)
plt.plot(X,Y,label="sol. by thomas algorithm")

for i in range(5):
    a=(n//4)*i
    plt.plot(X[a],Y[a],'*')

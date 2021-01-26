# -*- coding: utf-8 -*-
"""Untitled0.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1usZgCbCfF2S_SJy45PIyr0rO0KWHANb5
"""

#import library
import numpy as np
import matplotlib.pyplot as plt
import math

#load input
x=[]
y=[]
check=True
with open ('interpolation_input.txt','r') as f:
  for line in f:
    #print(line)
    try:
      (_x,_y)=line.split()
      x.append(float(_x))
      y.append(float(_y))
    except ValueError:
      check=False
      print("Oops!Check Input,plz...")
#print(x)
#print(y)
#check input
if check:
  if len(x)==len(y):
    print("Input Satisfied")
  else:
    print("Check Input")

#plot input (x,y)
#plt.plot(x,y,'*')
plt.scatter(x,y)

#define function
def u0(x):
  return 1
def u1(x):
  return x
def u2(x):
  return x**2
def u3(x):
  return x**3
def u4(x):
  return x**4
def u5(x):
  return math.sin(x)
def u6(x):
  return math.cos(x)
def u7(x):
  return math.sin(2*x)
def u8(x):
  return math.cos(2*x)
def u9(x):
  return math.exp(x)
def u10(x):
  return math.exp(-x)

#tra ve ma tran sau khi thay cac x vao function u  
def pack1(u,x):
    result=[]
    for i in range(len(u)):
        temp=list(map(u[i],x))
        result.append(temp)
    return np.array(result).T

def pack4(x) :
    change=[]
    for i in range(len(x)):
        temp1=math.log10(x[i])
        change.append(temp1)
    return np.array(change)

#pack2 : nhân ma trận nghich dao voi ma tran theta 
def pack2(theta):
    return theta.T@theta

#Dùng viền quanh để tính nghịch đảo của M    
# np.linalg.pinv(matrix) tinh nghich dao su dung thu vien co san
def vienquanh_inverse(A):
    n,_=A.shape
    if n==1:
        return 1/A
    elif n>1:
        start=1/A[0,0]
        for i in range(n-1):
            alpha11=start
            alpha12=A[:(i+1),i+1].reshape(i+1,1)
            alpha21=A[i+1,:i+1]
            alpha22=A[i+1,i+1]
            if i==0:
                X=alpha11*alpha12
            else :
                X=alpha11@alpha12
            if i==0:
                Y=alpha21*alpha11
            else :
                Y=alpha21@alpha11
                Y=Y.reshape(1,-1)
            if i==0:
                theta=alpha22-Y*alpha12
            else :
                theta=alpha22-Y@alpha12           
            if i==0:
                beta11=alpha11+(1/theta)*(X*Y)
            else :
                beta11=alpha11+(1/theta)*(X@Y)            
            beta12=-(1/theta)*X
            beta21=-(1/theta)*Y
            beta22=1/theta
            tempt_result=np.vstack((np.hstack((beta11,beta12)),np.hstack((beta21,beta22))))
            start=tempt_result
            np.savetxt('output.txt',tempt_result)    
        return tempt_result

def pack3(theta,M,y):
    return vienquanh_inverse(M)@theta.T@y
def find_y(x,u,a):
  y=0
  for i in range (0,len(u)):
    y=y+a[i]*u[i](x)
  return y

def graph(x,y,u,a):
  x_test=np.linspace(min(x),max(x),100000)
  y_test=[]
  for i in range (0,len(x_test)):
    y_test.append(find_y(x_test[i],u,a))
  plt.scatter(x,y,s=30,cmap='palete')
  plt.plot(x_test,y_test,'r')



u=[u0,u1]

theta=pack1(u,x)
print(theta)

M=pack2(theta)
print(M)

a=pack3(theta,M,y)
print('Ma tran he so la: ')
print(a)



graph(x,y,u,a)
plt.show()





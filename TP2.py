# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 08:18:14 2021

@author: basti
"""
from math import pi, cos, sin
import numpy as np
import matplotlib.pyplot as plt


#################################### Oscillations du pendule ####################################

g=9.81
L=1
wo=(g/L)**0.5
f=wo/(2*pi)
T=1/f

print ('La pulsation est de',wo,'rad/s')
print ('La fréquence est de',f,'Hz')
print ('La période est de',T,'s')

t = np.arange(0,4,0.04)   
A=pi/12 
Y0 = np.array([[A],[0]])    
h=0.04 
N=len(t)   

def teta_t():
    t=np.arange(0,4,0.04)
    teta=[]
    for i in range (len(t)):
        tet=(A)*cos(wo*t[i])
        teta.append(tet)
    return (teta)
#     plt.grid(True)
#     plt.plot(t,teta,label='Equation linéaire')
#     plt.title("Evolution de teta au cours du temps")
#     plt.xlabel('temps (s)')
#     plt.ylabel('teta (rad)')
#     plt.show()
    
   
    
def pendule (Y,t):
    Yprime=np.array([Y[1],-(g/L)*(Y[0])])
    return (Yprime)

def Euler(f,t,Y0,N,h):
    Ye=np.zeros((N,Y0.size))
    Ye[0,:]=Y0.reshape(2)
    for k in range (N-1):
        Ye[k+1,:]=Ye[k,:] + h*f(Ye[k,:],t[k])
    return (Ye)

    
    
def RK(f,t,Y0,N,h):
    Yrk=np.zeros((N,Y0.size))
    Yrk[0,:]=Y0.reshape(2)
    for k in range (N-1):
        k1=f(Yrk[k,:],t[k])
        k2=f(Yrk[k,:]+(h/2)*k1,t[k]+(h/2))
        k3=f(Yrk[k,:]+(h/2)*k2,t[k]+(h/2))
        k4=f(Yrk[k,:]+h*k3,t[k]+h)
        Yrk[k+1,:]=Yrk[k,:] + (h/6)*(k1+2*k2+2*k3+k4)
    return (Yrk)


from scipy.integrate import odeint

    
def solv_EDO():
    t = np.arange(0,4,0.04)   
    Y0 = np.array([A,0])    
    Yode =odeint(pendule, Y0, t)
    return Yode


def Courbes():
    teta_t()
    Ye= Euler(pendule,t,Y0,N,h)
    Yrk= RK(pendule,t,Y0,N,h)
    Yode =solv_EDO()
    plt.plot(t,Ye[:,0],label='Euler')
    plt.plot(t,Yrk[:,0],label='Runge Kutta')
    plt.plot(t,Yode[:,0],label='Odeint')
    plt.title("Evolution de teta au cours du temps") 
    plt.xlabel('temps (s)')
    plt.ylabel('teta (rad)')
    plt.grid(True)
    plt.legend()
    plt.show()
  

def Portrait_de_phase():
    
    Ye= Euler(pendule,t,Y0,N,h)
    Yrk= RK(pendule,t,Y0,N,h)
    Yode =solv_EDO()
    plt.plot(Ye[:,0],Ye[:,1],label='Euler')
    plt.plot(Yrk[:,0],Yrk[:,1],label='Runge Kutta')
    plt.plot(Yode[:,0],Yode[:,1],label='Odeint')
    plt.title("Portrait de Phase de teta") 
    plt.xlabel('teta (rad)')
    plt.ylabel("teta' (rad)")
    plt.legend()
    plt.grid()
    plt.show()
    
from statistics import mean 
def erreur_rel():
    Ye= Euler(pendule,t,Y0,N,h)
    Yrk= RK(pendule,t,Y0,N,h)
    Yode =solv_EDO()
    prec1=[]
    prec2=[]
    prec3=[]
    teta=[]
    err1=0
    err2=0
    err3=0
    for k in range (len(t)):
        tet=(A)*cos(wo*t[k])
        err1= (abs(Ye[k,0]-tet))/(abs(tet))
        err2= (abs(Yrk[k,0]-tet))/(abs(tet))
        err3= (abs(Yode[k,0]-tet))/(abs(tet))
        
        teta.append(tet)
        prec1.append(err1*100)
        prec2.append(err2*100)
        prec3.append(err3*100)
    return (prec1,prec2,prec3,teta)
# plt.subplot(3,1,1)
# prec1,prec2,prec3,teta=erreur_rel()
# print("Moyenne de l'erreur par Euler:",mean(prec1),"%")
# print("Moyenne de l'erreur par Runge Kutta:",mean(prec2),"%")
# print("Moyenne de l'erreur par Odeint:",mean(prec3),"%")
# plt.scatter(t,prec1,label='erreur Euler')
# plt.title("Evolution de l'erreur au cours du temps") 
# plt.legend()
# plt.ylabel('erreur (%)')
# plt.grid(True)
# plt.subplot(3,1,2)
# plt.scatter(t,prec2,label='erreur Runge Kutta')
# plt.legend()

# plt.ylabel('erreur (%)')
# plt.grid(True)
# plt.subplot(3,1,3)
# plt.scatter(t,prec3,label='erreur Odeint')

# plt.ylabel('erreur (%)')
# plt.xlabel('temps(s)')
# plt.grid(True)
# plt.legend()
# plt.show()

##########################Suspension du Véhicule####################################

M1 = 15 #kg; 
M2 = 200 #kg; 
C2 = 1200 #Ns/m; 
K1 = 50000 #N/m; 
K2 = 5000 #N/m:

x1_0 = 0 
x2_0 = 0 
x1_p0 = 0
x2_p0 = 0

f=-1000
  

def suspension (Y,t):
    Yprime=np.array([Y[2],Y[3], (1/M1)*(C2*Y[3]-C2*Y[2]-(K1+K2)*Y[0]+K2*Y[1]), (1/M2)*(C2*Y[2]-C2*Y[3]+K2*Y[0]-K2*Y[1]+f)])
    return (Yprime)

def solv_EDO2():
    t = np.arange(0,3,0.03)   
    Y0 = np.array([0,0,0,0])   
    Yode =odeint(suspension, Y0, t)
    n=len(Yode)
    lt1=[]
    lt2=[]
    lt3=[]
    lt4=[]
    x1=0.63*Yode[n-1,0]
    x2=0.63*Yode[n-1,1]
    print ("Tx1_95%=",0.95*0.13932,"Tx2_95%=",0.95*0.3140)
    for k in range (n):
        lt1.append(x1)
        lt2.append(x2)
        lt3.append(0.13932)
        lt4.append(0.3140)
    plt.scatter(0.13932,x1,c='r',label='temps caractéristique')
    plt.scatter(0.3140,x2,c='r')   
    plt.plot(lt3,Yode[:,1],c='black') 
    plt.plot(lt4,Yode[:,1],c='black') 
    plt.plot(t,lt1,c='black')    
    plt.plot(t,lt2,c='black')
    plt.plot(t,Yode[:,0],c='g',label='x1')
    plt.plot(t,Yode[:,1],c='r',label='x2')
    plt.title("Evolution de l'amortissement au cours du temps") 
    plt.xlabel('temps (s)')
    plt.ylabel('déplacement x ')
    plt.grid()
    plt.legend()
    plt.show()
   
    



        
  




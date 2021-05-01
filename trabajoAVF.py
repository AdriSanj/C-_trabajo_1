# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:35:55 2021

@author: adria
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

# FUNCIONES A IMPLEMENTAR

# Funciones de la solucion exacta del problema de Riemann

def h(x,t):
    H=np.zeros(len(x))
    for i in range(len(x)):
        if x[i]-L1*t<0:
            H[i]=hl
        elif (x[i]-L1*t)*(x[i]-L2*t)<=0:
            H[i]=alf1r+alf2l
        elif x[i]-L2*t>0:
            H[i]=hr
    return H

def u(x,t):
    U=np.zeros(len(x))
    for i in range(len(x)):
        if x[i]-L1*t<0:
            U[i]=ul
        elif (x[i]-L1*t)*(x[i]-L2*t)<=0:
            U[i]=(-alf1r+alf2l)*sqrt(g/ht)          ## Si le pongo un - sale guay
        elif x[i]-L2*t>0:
            U[i]=ur
    return U

# ------------------------------------------------#
# Funciones para la solucion del problema de Cauchy

def v1_pulso(x,t):          ## IMPORTANTE: buscar otra condicion inicial de Cauchy
    return np.exp(-2*(x-L1*t)**2)

def v2_pulso(x,t):
    return np.exp(-2*(x-L2*t)**2)
    
def h_C(x,t):
    return v1_pulso(x,t)+v2_pulso(x,t)          # Aqui h es un pulso
    
def u_C(x,t):
    return (-v1_pulso(x,t)+v2_pulso(x,t))*sqrt(g/ht)        # En el momento inicial esto es cero


# ------------------------------------------------#
# ------------------------------------------------#

# Elegir el apartado que se quiere representar
apartado="c"

# Constantes del problema
ut = 0.5          # u tilde
ht = 7         # h tilde. Ha de ser siempre positiva
g = 9.8         # gravedad g

# Condiciones iniciales
ur = 0
ul = 0
hr = 5
hl = 10

# Autovalores 
L1= ut-sqrt(ht/g)       # Tiene sentido la forma de la grafica ya que no son autovalores simetricos, a diferencia del caso de la ecuacion de ondas
L2= ut+sqrt(ht/g)


# Valores de las z iniciales
## Se corresponden con las funciones v1 y v2 soluciones de V = P**-1 W
alf1r=0.5*(hr-ur*sqrt(ht/g))            
alf2l=0.5*(hl+ul*sqrt(ht/g))

# Representacion grafica
x=np.linspace(-10,10,100)
t=0

if apartado=='c':
    while t<10:
        plt.subplot(211)
        plt.axis('equal')
        plt.cla()           # plt.cla() borra los ejes correspondientes
        plt.grid()
        plt.plot(x,h(x,t))
        plt.title("h(x,t) vs x")
        plt.ylabel("h (m)")
        plt.xlabel("x (m)")
        plt.draw()
        plt.subplot(212)
        plt.axis('equal')
        plt.cla()
        plt.grid()
        plt.plot(x,u(x,t))
        plt.title("h(x,t) vs x")
        plt.ylabel("u (m*s-1)")
        plt.xlabel("x (m)")
        plt.draw()
        t+=1
        plt.pause(0.1)
elif apartado=='d':
    while t<10:
        plt.subplot(211)
        plt.axis('equal')
        plt.cla()           
        plt.grid()
        plt.plot(x,h_C(x,t))
        plt.title("h(x,t) vs x")
        plt.ylabel("h (m)")
        plt.xlabel("x (m)")
        plt.draw()
        plt.subplot(212)
        plt.axis('equal')
        plt.cla()
        plt.grid()
        plt.plot(x,u_C(x,t))
        plt.title("h(x,t) vs x")
        plt.ylabel("h (m)")
        plt.xlabel("u (m*s-1)")
        plt.draw()
        t+=1
        plt.pause(0.1)
elif apartado=='er':
    u_i=np.zeros(len(x))
    p_i=np.zeros(len(x))
    u_i_1=np.zeros(len(x))
    p_i_1=np.zeros(len(x))
    for i in range(len(x)):
        if x[i]<0:
            u_i_1[i]=ul
            p_i_1[i]=pl
        elif x[i]>=0:
            u_i_1[i]=ur
            p_i_1[i]=pr
    u_i=np.copy(u_i_1)
    p_i=np.copy(p_i_1)
    plt.subplot(211)
    plt.cla()
    plt.grid()
    plt.plot(x,u_i)
    plt.title("u(x,t) vs x ")
    plt.ylabel("u (m * s**-1)")
    plt.draw()
    plt.subplot(212)
    plt.cla()
    plt.grid()
    plt.plot(x,p_i)
    plt.title("p(x,t) vs x ")
    plt.ylabel("p (N*m**-2)")
    plt.xlabel("x (m)")
    plt.draw()
    while t<10:
        for i in range(1,len(x)-1):     ## Estas formulas son donde esta el error
            u_i_1[i]=u_i[i]*(1-Dt*c0/Dx)-(0.5*Dt/Dx)*((p_i[i+1]-p_i[i-1])/rho -c0*(u_i[i+1]+u_i[i-1]))
            p_i_1[i]=p_i[i]*(1-Dt*c0/Dx)-(0.5*Dt/Dx)*((u_i[i+1]-u_i[i-1])*rho*c0**2 -c0*(p_i[i+1]+p_i[i-1]))
        # Condiciones de contacto sobre los extremos
        u_i_1[0]=np.copy(u_i_1[1])
        u_i_1[len(x)-1]=np.copy(u_i_1[len(x)-2])
        p_i_1[0]=np.copy(p_i_1[1])
        p_i_1[len(x)-1]=np.copy(p_i_1[len(x)-2])
        plt.subplot(211)
        plt.cla()
        plt.grid()
        plt.plot(x,u_i_1)
        if comparacion=="si":
            plt.plot(x,u(x,t))
        plt.title("h(x,t) vs x ")
        plt.ylabel("h (m)")
        plt.draw()
        plt.subplot(212)
        plt.cla()
        plt.grid()
        plt.plot(x,p_i_1)
        if comparacion=="si":
            plt.plot(x,p(x,t))
        plt.title("u(x,t) vs x ")
        plt.ylabel("u (m*s-1")
        plt.xlabel("x (m)")
        plt.draw()
        plt.pause(0.1)
        u_i=np.copy(u_i_1)
        p_i=np.copy(p_i_1)
        t+=1
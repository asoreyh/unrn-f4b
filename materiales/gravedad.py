#!/usr/bin/python
# -*- coding: utf8 -*-

"""
H. Asorey
Cálculo del campo gravitatorio tridimensional para un objeto central utilizando una clase propia de vectores. Debería usarse numpy pero se trata de un ejemplo para aprender a trabajar con clases en python  
"""
import math

####################### mi clase vector

class vector:
    def __init__(self,lista):
        self.x=[]
        for xi in lista:
          self.x.append(float(xi))
        self.dim=len(self.x)
        self.mod=math.sqrt(self.prod_escalar(self))

    def prod_escalar(self,a):
        if (self.dim == a.dim):
            pesc=0.
            for i in range(0,self.dim):
                pesc += (self.x[i] * a.x[i])
            return pesc
        else:
            print ("El productor escalar se define para vectores de la misma dimension")
            return None

###################### Operaciones entre vectores

def suma(a, b):
    suma=[]
    if (a.dim == b.dim):
       for i in range(0,a.dim):
          suma.append(a.x[i]+b.x[i])
       return vector(suma)
    else:
       return None

def resta(a, b):
    resta=[]
    if (a.dim == b.dim):
       for i in range(0,a.dim):
          resta.append(a.x[i]-b.x[i])
       return vector(resta)
    else:
       return None

def vector_escalar(a, num):
    escalar=[]
    for i in range(0,a.dim):
      escalar.append(num*a.x[i])
    return vector(escalar)


################ Unidades
## Hay que pasar todo a unidades del SI
## Usamos múltiplos comunes (1km = 1000m; 1 dia=86400 s; etc)

kg = 1.
s = 1.
m = 1.
G = 6.67e-11 * m**3 / (kg * s**2)

## algunas definiciones útiles
km = 1000. * m
ua = 1.5e8 * km
# masa solar
ms = 1.989e30 * kg 
# masa tierra
mt = 5.972e24 * kg
# masa júpiter
mj = 1.898e27 * kg

# la masa de prueba
mp = 1. * kg

## las masas de campo
## son las que producen el campo gravitatorio
## por ejemplo, dos estrellas con masas solares separadas por 1 UA centradas en el eje x 
m1 = 1. * ms
r1=vector([0.5*ua, 0, 0])

m2 = 1. * ms 
r2=vector([-0.5*ua, 0, 0])

# auxiliares
# r es el vector posición donde calculo el campo  
r=vector([0, 0, 0])

# dimensiones de la grilla cartesiana. 
# base es la unidad base para la grilla cartesiana. 
base = 1. * ua
# res es la cantidad de puntos por unidad de base. Por ejemplo 10 significa que tomará 10 puntos por unidad base, en dr
res=1/2.
# límite es el máximo (y mínimo) de la grilla en bases. La grilla será un cubo de 2*lim*base de lado
lim = 3.

# paso
dr = base / res

# recorro el espacio
for x in range(int(-lim*res),int((lim*res)+1)):
  for y in range(int(-lim*res),int((lim*res)+1)):
    for z in range(int(-lim*res),int((lim*res)+1)):
      # estoy mirando el punto r
      r=vector([x * dr, y * dr, z * dr])
      # calculo los vectores de r respecto a cada estrella
      d1=resta(r, r1)
      d2=resta(r, r2)

      # si ninguno es 0 (no estoy en la estrella donde el campo diverge), calculo:
      if (d1.mod and d2.mod):
        # Fuerza debida a la estrella 1, va distancia cubo xq uso el vector y no el versor
        F1_mod = -G * m1 * mp / d1.mod**3
        # La fuerza tiene dirección de d1 y módulo F1        
        F1 = vector_escalar(d1, F1_mod)
        # Fuerza debida a la estrella 2, va distancia cubo xq uso el vector y no el versor
        F2_mod = -G * m2 * mp / d2.mod**3
        F2 = vector_escalar(d2, F2_mod)

        # Superposición -> suma vectorial de las fuerzas
        F=suma(F1, F2)

        # ahora las energías potenciales gravitatorias
        U1 = -G * m1 *mp / d1.mod
        U2 = -G * m2 *mp / d2.mod
        U=U1+U2

        # pone el vector en el vector r para dibujar en cada punto
        Faux = suma (r, F)

        # imprimo los resultados: posición, fuerza, módulo y energía para este punto y vuelvo. 
        # Al terminar sale
        for i in range (0, r.dim):
          print (r.x[i], end=',')
        for i in range (0, r.dim): 
          print (Faux.x[i], end=',')
        print (F.mod, end=',')
        print (U)


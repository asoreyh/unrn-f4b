#!/usr/bin/python
# -*- coding: utf8 -*-

"""
H. Asorey
1) Cálculo de la órbita del exoplaneta elegido utilizando una clase propia de vectores. Debería usarse numpy pero se trata de un ejemplo para aprender a trabajar con clases en python  
"""

import math

####################### mi clase vector

class vector:
    def __init__(self,lista):
        self.componentes=[]
        for componente_i in lista:
          self.componentes.append(float(componente_i))
        self.dim=len(self.componentes)
        self.mod=math.sqrt(prod_escalar(self,self))

###################### Operaciones entre vectores

def prod_escalar(a,b):
    if (a.dim == b.dim):
        pesc=0.
        for i in range(0,a.dim):
            pesc += (a.componentes[i] * b.componentes[i])
        return pesc
    else:
        print("El productor escalar se define para vectores de la misma dimension")
        return None

def coseno(a,b):
    if (a.dim == b.dim):
      mods=a.mod*b.mod
      if (mods != 0 ):
        return (prod_escalar(a,b) / mods)
      else:
        print("Para calcular el coseno los vectores no pueden ser nulos")
        return None
    else:
        print("Sólo se calcula el coseno para vectores de la misma dimension")
        return None

def suma(a, b):
    suma=[]
    if (a.dim == b.dim):
       for i in range(0,a.dim):
          suma.append(a.componentes[i]+b.componentes[i])
       return vector(suma)
    else:
       print ("La suma se define para vectores de la misma dimension")
       return None

def resta(a, b):
    resta=[]
    if (a.dim == b.dim):
       for i in range(0,a.dim):
          resta.append(a.componentes[i]-b.componentes[i])
       return vector(resta)
    else:
       print ("La resta se define para vectores de la misma dimension")
       return None

################################## Producto por un escalar

def vector_escalar(a, num):
    escalar=[]
    for i in range(0,a.dim):
      escalar.append(num*a.componentes[i])
    return vector(escalar)

################ Unidades
## Hay que pasar todo a unidades del SI
## Usamos múltiplos comunes (1km = 1000m; 1 dia=86400 s; etc)

#unidades y constantes fundamentales
s=1.
kg=1.
m=1.
G = 6.67e-11 * m**3 / (kg * s**2) 

#múltiplos
km=1000. * m
dia=86400. * s

#referencias
masa_jupiter= 1.898e27 * kg
masa_sol= 1.989e30 * kg
ua=1.5e8 * km
radio_sol = 695500. * km
radio_jup = 69911 * km

################ clase exoplaneta
## ejemplo de clase para contener los datos de los planetas que podríamos necesitar para el cálculo de la órbita

class exoplaneta:
  def __init__(self,lista):
    if (lista[0] != ''):
        self.nombre=lista[0]
    else:
        self.nombre=None

    if (lista[1] != ''):
        self.semieje=float(lista[1])
    else:
        self.semieje=-1.

    if (lista[2] != ''):
        self.periodo=float(lista[2])
    else:
        self.periodo=-1.
    
    if (lista[3] != ''):
        self.exc=float(lista[3])
    else:
        self.exc=-1.
    
    if (lista[4] != ''):
        self.masa=float(lista[4])
    else:
        self.masa=-1.
    
    if (lista[5] != ''):
        self.masae=float(lista[5])
    else:
        self.masae=-1.
    
    if (lista[6] != ''):
        self.radioe=float(lista[6])
    else:
        self.radioe=-1.
    
    if (lista[7] != ''):
        self.tempe=float(lista[7])
    else:
        self.tempe=-1.
    
    if (lista[8] != ''):
        self.radio=float(lista[8])
    else:
        self.radio=-1.

# Características del exoplaneta a simular

nombre         = "Test"
excentricidad  = 0.
semieje_mayor  = 1.  * ua
masa_estrella  = 1.  * masa_sol
masa_planeta   = 1.  * masa_jupiter
radio_extrella = 1.  * radio_sol
radio_planeta  = 1.  * radio_jup

# número de pasos para completar una órbita
intervalos = 1000 
# cantidad de órbitas a completar
n = 1

#lista de características
mi_planeta=[nombre, semieje_mayor, -1., excentricidad, masa_planeta, masa_estrella, radio_extrella, -1., radio_planeta]

# elemento exoplaneta
myexo=exoplaneta(mi_planeta)

##### CONDICIONES INICIALES
# Calculo algunos parámetros orbitales. Trabajo en el SI. Con visviva calculo la velocidad en el apoastro
# Trabajo en 2D porque por Kepler las órbitas están en un plano
foco = myexo.semieje * myexo.exc
apoastro = myexo.semieje + foco
visviva = math.sqrt(G * myexo.masae * ((2./apoastro)-(1./myexo.semieje)))
periodo = 2*math.pi*math.sqrt(myexo.semieje**3 / (G*myexo.masae) )
# intervalo temporal
dt = periodo / intervalos 
# posicion inicial, empiezo en el apoastro
r=vector([-apoastro,0])
r_unit=vector([0,0])
# velocidad inicial, perpendicular a r, lo pongo en la dirección 'y'
v=vector([0, visviva])  
# aceleración inicial, entramos con 0 y se calcula luego del avance v*dt a la nueva posición orbital 
a=vector([0,0])
amod = 0. 

#inicia el loop
for i in range(0, int(n*intervalos)):
    # imprimo el paso i, el tiempo transcurrido en días y la posición del planeta (en unidades astronomicas) 
    print (i, i*dt/dia, r.componentes[0]/ua, r.componentes[1]/ua)
    # Algoritmo Newton-Hooke
    # actualizo la posición r_{i+1} = r_i + dt * v_i 
    r=suma(r, vector_escalar(v,dt))
    # determino un vector con dirección r, sentido hacia dentro y módulo 1
    if (r.mod):
      r_unit=vector_escalar(r,(-1.)/r.mod)
    else:
      print ("Hay un problema: |r| = 0")
      exit()

    # ahora, calculo la aceleración
    # en kepler, la aceleración centrípeta está dada por la gravedad 
    # que tiene dirección radial, sentido hacia dentro y módulo 
    # |a| = G * M / |r|^2
    amod = G * myexo.masae / (r.mod**2)
    # entonces la aceleración es:
    # a = amod * r_unit 
    a=vector_escalar(r_unit,amod)
    # y actualizo la velocidad v_{i+1} = v_i + a_i dt
    v=suma(v,vector_escalar(a,dt))
    # e itero

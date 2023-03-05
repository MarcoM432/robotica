#Programa1_sesion1
import numpy as np

#creacion de un vector
a = np.array([[1,1,0]]) #Vectores fila 
b = np.array([[1],
              [1],
              [0]]) #vector columna

c = np.array([ [1,0,0],
              [0,1,0],
              [0,1,1]])
print(a@b)
print(c)
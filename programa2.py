import numpy as np
import math



R = 0.03 #radio en metros de la rueda 
L = 0.18 #distancia entre las ruedas en metros
wR = 500  # velocidad de la rueda derecha (radianes por segundo)
wL = -500  #velocidad de la rueda izquierda
teta = math.pi/2

a = np.array([  [R*math.cos(teta)/2, R*math.cos(teta)/2],
                [R*math.sin(teta)/2, R*math.sin(teta)/2],
                [R/L          , -R/L         ]  ])

b = np.array([  [wR],
                [wL]    ])

c = a@b
print(c)
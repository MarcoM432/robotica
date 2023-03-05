import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import time
import math
import numpy as np
import matplotlib.pyplot as plt

pi=math.pi
#FUNCIONES MATEMATICAS 
def cos(A):
    return(math.cos(A))

def sin(A):
    return(math.sin(A))

def integrador_euler(y_i , h, d_y_i):
    y_i_1= y_i+h*(d_y_i)
    return y_i_1
def derivador(y0 , h, yf):
    if (h>0):
        dx= (yf-y0)/h
        return dx 
    else:
        dx= (yf-y0)/0.0001
        return dx

#FUNCIONES DE DIBUJO
def Draw_line(A,B,grosor,color,escala):
    r,g,b=color
    glLineWidth(grosor)#grosor de linea
    glBegin(GL_LINES)#inicia dibujado
    glColor(r,g,b)#color de la linea
    glVertex3d(A[0]*escala, A[1]*escala, A[2]*escala)#punto inicial 
    glVertex3d(B[0]*escala, B[1]*escala, B[2]*escala)#punto inicial 
    glEnd()
    glColor(1, 1, 1)
    glLineWidth(1)

def Draw_circuferencia_plano_z(A,radio,divisiones,grosor,color,escala):
    r,g,b=color
    angulo=2*pi/divisiones
    lista=[]
    glLineWidth(grosor)#grosor de linea
    glBegin(GL_LINES)#inicia dibujado
    glColor(r,g,b)#color de la linea
    
    for i in range(divisiones+1):
        x=A[0]+radio*cos(i*angulo)
        y=A[1]+radio*sin(i*angulo)
        z=A[2]
        lista.append([x,y,z])
    for j in range(divisiones):
        glVertex3d(lista[j][0]*escala, lista[j][1]*escala, lista[j][2]*escala)#punto inicial 
        glVertex3d(lista[j+1][0]*escala, lista[j+1][1]*escala, lista[j+1][2]*escala)#punto inicial 
        
    glEnd()
    glColor(1, 1, 1)
    glLineWidth(1)

def Draw_esfera(A,radio,color,escala):
    r,g,b=color
    sphere = gluNewQuadric()
    glColor(r,g,b)
    glPushMatrix()
    glTranslated(A[0]*escala, A[1]*escala, A[2]*escala)#En donde queremos colocal la esfera
    #                 radio,divisiones,divisiones
    gluSphere(sphere, radio*escala, 15, 15) 
    glPopMatrix()
    glColor(1, 1, 1)
    glLineWidth(1)

def Draw_cartesiano(escala):
    Draw_esfera([0,0,0],0.02,(1,1,1),escala)
    Draw_line([0,0,0],[5,0,0],1,(0,1,0),escala)
    Draw_line([0,0,0],[0,5,0],1,(0,0,1),escala)
    Draw_line([0,0,0],[0,0,5],1,(1,0,0),escala)

def Draw_Robug(Robug,escala):
    A=[Robug.x,Robug.y,Robug.z]
    teta=Robug.teta 
    radio=Robug.rb
    rw=Robug.rw 
    color=Robug.color  
    #coordenadas del centroide
    x=A[0];y=A[1];z=A[2]
    #centroide de la primera esfera
    x_1=x+radio*cos(teta)
    y_1=y+radio*sin(teta)
    z_1=z+rw
    #centroide de la segunda esfera
    x_2=x+radio*cos(teta+2*math.pi/3)
    y_2=y+radio*sin(teta+2*math.pi/3)
    z_2=z+rw
    #centroide de la segunda esfera
    x_3=x+radio*cos(teta+2*2*math.pi/3)
    y_3=y+radio*sin(teta+2*2*math.pi/3)
    z_3=z+rw
    #Dibujamos las rectas
    Draw_line([x_1,y_1,z_1],[x,y,z+rw],1.5*escala,(color[0]*.5,color[1]*.5,color[2]*.5),escala)
    Draw_line([x_2,y_2,z_2],[x,y,z+rw],1.5*escala,(color[0]*.5,color[1]*.5,color[2]*.5),escala)
    Draw_line([x_3,y_3,z_3],[x,y,z+rw],1.5*escala,(color[0]*.5,color[1]*.5,color[2]*.5),escala)
    #Ahora dibujaremos el contorno
    Draw_circuferencia_plano_z([x,y,z+rw],radio,20,1.5*escala,color,escala)
    #ahora dibujamos las ruedas como esferas
    Draw_esfera([x_1,y_1,z_1],rw,color,escala)
    Draw_esfera([x_2,y_2,z_2],rw,color,escala)
    Draw_esfera([x_3,y_3,z_3],rw,color,escala)
    Draw_esfera([x_1,y_1,z_1],rw/1.5,(1,1,1),escala)
    
def Draw_robot_diferencial(Robot,escala):
    #coordenadas del centroide
    x=Robot.x;y=Robot.y;z=Robot.z
    teta=Robot.teta
    L=Robot.L
    rw=Robot.radio
    teta=Robot.teta
    color=Robot.color 
    #centroide de la primera esfera (LLANTA DERECHA)
    x_1=x+(L/2)*sin(teta)
    y_1=y-(L/2)*cos(teta)
    z_1=z+rw
    #centroide de la segunda esfera (LLANTA IZQUIERDA)
    x_2=x-(L/2)*sin(teta)
    y_2=y+(L/2)*cos(teta)
    z_2=z+rw
    #Dibujamos las rectas
    Draw_line([x_1,y_1,z_1],[x_2,y_2,z_2],1.5*escala,(color[0]*.5,color[1]*.5,color[2]*.5),escala)
    #Dibujamos la circunferencia
    Draw_circuferencia_plano_z([x,y,z+rw],L/2,20,1.5*escala,color,escala)
    #ahora dibujamos las ruedas como esferas
    Draw_esfera([x_1,y_1,z_1],rw,color,escala)
    Draw_esfera([x_2,y_2,z_2],rw,color,escala)
    Draw_esfera([x+(L/2)*cos(teta),y+(L/2)*sin(teta),z+rw],rw/1.5,(1,1,1),escala)

def draw_drone(Drone1,escala):
    #Primero definimos la geometría del robot
    #4 puntos 
    A=np.array([[Drone1.radio],
              [ 0 ],
              [ 0 ],
              [ 1 ]])#[ .3,0,0] # en metros
    B=np.array([[-Drone1.radio],
              [ 0 ],
              [ 0 ],
              [ 1 ]])#[ .3,0,0] # en metros
    C=np.array([[0  ],
              [ Drone1.radio],
              [ 0 ],
              [ 1 ]])#[ .3,0,0] # en metros
    D=np.array([[ 0 ],
              [-Drone1.radio],
              [ 0 ],
              [ 1 ]])#[ .3,0,0] # en metros
    alfa=Drone1.roll 
    beta=Drone1.pitch 
    gama=Drone1.yaw 
    #Creamos matríz de transformación local
    M_rx=np.array([[1,    0    ,    0      ],
                [0,cos(alfa), -sin(alfa)],
                [0,sin(alfa), cos(alfa) ]])
    M_ry=np.array([[cos(beta),  0 , sin(beta)],
                 [0        ,  1 ,    0     ],
                 [-sin(beta), 0 , cos(beta) ]])
    M_rz=np.array([[cos(gama), -sin(gama), 0  ],
                [sin(gama),  cos(gama), 0  ],
                [  0      ,  0        , 1  ]])
    M_rotacion=M_rx@M_ry@M_rz
    M_Homogenea=np.array([[M_rotacion[0][0],M_rotacion[0][1],M_rotacion[0][2],0],
                      [M_rotacion[1][0],M_rotacion[1][1],M_rotacion[1][2],0],
                      [M_rotacion[2][0],M_rotacion[2][1],M_rotacion[2][2],0],
                      [    0           ,       0        ,        0       ,1]])
    lx=Drone1.x;
    ly=Drone1.y;
    lz=Drone1.z 
    M_Posicion=np.array([[1,0,0,lx],
                     [0,1,0,ly],
                     [0,0,1,lz],
                     [0,0,0,1 ]])
    A=M_Homogenea@A
    A2=M_Posicion@A  
    B=M_Homogenea@B 
    B2=M_Posicion@B
    C=M_Homogenea@C 
    C2=M_Posicion@C
    D=M_Homogenea@D 
    D2=M_Posicion@D


    #puntos de los motores del cuadricoptero
    vectorX=np.array([[A2[0]],
                   [A2[1]],
                   [A2[2]],
                   [1]])
    vectorX2=np.array([[B2[0]],
                   [B2[1]],
                   [B2[2]],
                   [1]])
    vectorX3=np.array([[C2[0]],
                   [C2[1]],
                   [C2[2]],
                   [1]])
    vectorX4=np.array([[D2[0]],
                   [D2[1]],
                   [D2[2]],
                   [1]])
    #Comienza dibujo
    Draw_line(vectorX,vectorX2,1.5*escala,Drone1.color,escala)
    Draw_line(vectorX3,vectorX4,1.5*escala,Drone1.color,escala)
    #Draw_circuferencia_plano_z([Drone1.x,Drone1.y,Drone1.z],Drone1.radio*escala,20,1*escala,Drone1.color,escala)
    alfa=Drone1.roll
    beta=Drone1.pitch 
    gama=Drone1.yaw 
    Draw_circuferencia3D([Drone1.x,Drone1.y,Drone1.z],alfa,beta,gama,Drone1.radio,20,1*escala,(1,1,1),escala)
    Draw_esfera([A2[0],A2[1],A2[2]],0.02,(1,1,1),escala)

def Draw_circuferencia3D(B,alfa,beta,gama,radio,divisiones,grosor,color,escala):
    r,g,b=color
    angulo=2*pi/divisiones
    lista=[]
    glLineWidth(grosor)#grosor de linea
    glBegin(GL_LINES)#inicia dibujado
    glColor(r,g,b)#color de la linea
    #Creamos matríz de transformación local
    M_rx=np.array([[1,    0    ,    0      ],
                       [0,cos(alfa), -sin(alfa)],
                      [0,sin(alfa), cos(alfa) ]])
    M_ry=np.array([[cos(beta),  0 , sin(beta)],
                     [0        ,  1 ,    0     ],
                     [-sin(beta), 0 , cos(beta) ]])
    M_rz=np.array([[cos(gama), -sin(gama), 0  ],
                    [sin(gama),  cos(gama), 0  ],
                    [  0      ,  0        , 1  ]])
    M_rotacion=M_rx@M_ry@M_rz
    M_Homogenea=np.array([[M_rotacion[0][0],M_rotacion[0][1],M_rotacion[0][2],0],
                            [M_rotacion[1][0],M_rotacion[1][1],M_rotacion[1][2],0],
                            [M_rotacion[2][0],M_rotacion[2][1],M_rotacion[2][2],0],
                            [    0           ,       0        ,        0       ,1]])
    M_Posicion=np.array([[1,0,0,B[0]],
                         [0,1,0,B[1]],
                         [0,0,1,B[2]],
                         [0,0,0,1 ]])
            
    for i in range(divisiones+1):
        x=radio*cos(i*angulo)
        y=radio*sin(i*angulo)
        z=0
        A=np.array([[x],
                    [y],
                    [z],
                    [1]])
        A=M_Homogenea@A #Primero Rotamos el punto
        A2=M_Posicion@A #Luego Trasladamos a la posicion del Drone
        #Nueva posición del punto
        vectorX=np.array([[A2[0]*escala],
                          [A2[1]*escala],
                          [A2[2]*escala],
                          [1]])
        lista.append([vectorX[0][0],vectorX[1][0],vectorX[2][0]])
    for j in range(divisiones):
        glVertex3d(lista[j][0], lista[j][1], lista[j][2])#punto inicial 
        glVertex3d(lista[j+1][0], lista[j+1][1], lista[j+1][2])#punto inicial 
        
    glEnd()
    glColor(1, 1, 1)
    glLineWidth(1)

#FUNCIONES DE CONTROL
def cinematica_directa_diferencial(wL,wR,robot):
    R=robot.radio
    L=robot.L
    teta=robot.teta
    A=np.array([[R*cos(teta)/2, R*cos(teta)/2],
                 [R*sin(teta)/2, R*sin(teta)/2],
                 [R/L,-R/L]])
    B=np.array([[wR],
                [wL]])
    C=A@B
    return (C)
def Calcula_angulo_vector(A):
    norma=math.sqrt((A[0][0]*A[0][0]) +(A[1][0]*A[1][0]))
    x=A[0][0]
    y=A[1][0]
    #1er cuadrante
    if ((x>=0 ) & (y>=0)):
        teta=math.acos(x/norma)
        return teta
    #2do cuadrante
    if ((x<=0 ) & (y>=0)):
        fi=math.acos(abs(x)/norma)
        teta=math.pi - fi
        return teta
    #3er cuadrante
    if ((x<=0 ) & (y<=0)):
        fi=math.acos(abs(x)/norma)
        teta=math.pi + fi
        return teta
    #4to cuadrante
    if ((x>=0 ) & (y<=0)):
        fi=math.acos(x/norma)
        teta=(2*math.pi) - fi
        return teta
def filtra_angulo(teta):
    while True:
        if abs(teta)> (2*math.pi):
            if teta>0:
                teta=teta- (2*math.pi)
            else:
                teta=teta + (2*math.pi)
        else:
            break
    #En este momento ya es un valor entre -2pi y 2 pi
    if teta>=0:
        return teta
    else:
        #hay que convertirlo a positivo sumandole una revolucion
        teta=teta+ (2*math.pi)
        return teta 
def norma_vector(A):
    x=A[0][0]
    y=A[1][0]
    return math.sqrt((x*x) +(y*y))
def producto_cruz_2_vectores(A,B):
    cruz=A[0][0]*B[1][0] - (A[1][0]*B[0][0])
    return cruz    

def Control_de_robot_diferencial(Robot,dt,Punto_deseado):
    #Primero Orientacion
    B=np.array([[Punto_deseado[0]],
                [Punto_deseado[1]]])
    A=np.array([[Robot.x],
                [Robot.y]])
    vec_direccion=B-A
    norma=norma_vector(vec_direccion)
    teta_deseado=Calcula_angulo_vector(vec_direccion)
    #Ahora empezaremos nuestro esquema de control

    #Avanzar hacia delante: 
            #calcular d
    vec_ang_deseado=np.array([[cos(teta_deseado)],
                              [sin(teta_deseado)]])
    vec_ang_robot=np.array([[cos(Robot.teta)],
                            [sin(Robot.teta)]])
    d=vec_ang_deseado-vec_ang_robot
    norma_d=norma_vector(d)
    wL=0 #valores iniciales
    wR=0
    if norma_d<= math.sqrt(2): #El movimiento sera hacia adelante
        #Iremos hacia adelante 
            #En que sentido???
        sentido=producto_cruz_2_vectores(vec_ang_robot,d)
        k_teta=20
        Kp=20
        error_posicion=norma
        if sentido>=0:
            #sentido antihorario
            wL=-norma_d*k_teta 
            wR=norma_d*k_teta 
        else:
            wL=norma_d*k_teta 
            wR=-norma_d*k_teta
        if abs(norma)<0.005:
            norma=0
        if (norma_d<0.3):
            wL=Kp*norma
            wR=Kp*norma
    
    else: #El movimiento sera hacia atras
        #Iremos hacia atras
            #En que sentido???
        vec_ang_deseado=np.array([[cos(teta_deseado+pi)],
                              [sin(teta_deseado+pi)]])
        d=vec_ang_deseado-vec_ang_robot
        norma_d=norma_vector(d)
        sentido=producto_cruz_2_vectores(vec_ang_robot,d)
        k_teta=20
        Kp=20
        error_posicion=norma
        if sentido>=0:
            #sentido antihorario
            wL=-norma_d*k_teta 
            wR=norma_d*k_teta 
        else:
            wL=norma_d*k_teta 
            wR=-norma_d*k_teta
        if abs(norma)<0.005:
            norma=0
        if (norma_d<0.3):
            wL=-Kp*norma
            wR=-Kp*norma
        
        #nos desplazamos 
                
    dx=cinematica_directa_diferencial(wL,wR,Robot)
    Robot.x=integrador_euler(Robot.x,dt,dx[0][0])
    Robot.y=integrador_euler(Robot.y,dt,dx[1][0])
    Robot.teta=integrador_euler(Robot.teta,dt,dx[2][0])
    #print (Robot.x,Robot.y,Robot.teta)
def cinematica_inversa_omni(dx_local,dy_local,dteta_local,rb,rw):
    Beta1=0;
    Beta2=Beta1 + 2*pi/3
    Beta3=Beta2 + 2*pi/3
    
    A=np.array([[sin(Beta1), -cos(Beta1),-rb],
                 [sin(Beta2), -cos(Beta2),-rb],
                 [sin(Beta3), -cos(Beta3),-rb]])
    B=np.array([[dx_local],
                [dy_local],
                [dteta_local]])
    C=(1/rw)*(A@B)

    return (C)
def transforma_global_a_local_omnidireccional(Robot_color,dx_global):
    RGL=np.array([[cos(Robot_color.teta), -sin(Robot_color.teta),0],
                  [sin(Robot_color.teta), cos(Robot_color.teta) , 0],
                  [           0      , 0                  , 1]])
    #Transformacion
    velocidades_locales=np.linalg.inv (RGL) @ dx_global
    return velocidades_locales

def cinematica_directa_omni(rb,rw,dfi_1,dfi_2,dfi_3):
    Beta1=0;
    Beta2=Beta1 + 2*pi/3
    Beta3=Beta2 + 2*pi/3
    A=np.array([[sin(Beta1), -cos(Beta1),-rb],
                 [sin(Beta2), -cos(Beta2),-rb],
                 [sin(Beta3), -cos(Beta3),-rb]])
    B=np.array([[dfi_1],
                [dfi_2],
                [dfi_3]])
    C=(rw)*np.linalg.inv(A)@B

    return (C)
def Control_PID_OMNIDIRECCIONAL(trayectoria_deseada,Robug,kp,ki,kd,dt):
        Robug.vec_x.append(float(Robug.x))
        Robug.vec_y.append(float(Robug.y))
        Robug.vec_teta.append(float(Robug.teta))

        P_deseada=trayectoria_deseada
        #Calculo del error x
        error_robug_x  = P_deseada[0][0] - Robug.x
        d_error_robug_x= derivador(Robug.error_x_pasado , dt, error_robug_x)
        Robug.error_x_pasado=error_robug_x #actualizamos la memoria
        Robug.i_error_x=integrador_euler(Robug.i_error_x,dt,error_robug_x)
        PID_x=kp*error_robug_x +ki*Robug.i_error_x + kd*d_error_robug_x
        #Calculo del error y
        error_robug_y  = P_deseada[1][0] - Robug.y
        d_error_robug_y= derivador(Robug.error_y_pasado , dt, error_robug_y)
        Robug.error_y_pasado=error_robug_y #actualizamos la memoria
        Robug.i_error_y=integrador_euler(Robug.i_error_y,dt,error_robug_y)
        PID_y=kp*error_robug_y +ki*Robug.i_error_y + kd*d_error_robug_y
        #Calculo del error teta
        error_robug_teta  = P_deseada[2][0] - Robug.teta
        d_error_robug_teta= derivador(Robug.error_teta_pasado , dt, error_robug_teta)
        Robug.error_teta_pasado=error_robug_teta #actualizamos la memoria
        Robug.i_error_teta=integrador_euler(Robug.i_error_teta,dt,error_robug_teta)
        PID_teta=kp*error_robug_teta +ki*Robug.i_error_teta + kd*d_error_robug_teta

        PID=np.array([[PID_x],
                      [PID_y],
                      [PID_teta]])
        #El resultado del PID son velocidades cartesianas globales
        #Las transformamos a locales para el robot 
        dx_locales = transforma_global_a_local_omnidireccional(Robug,PID)
        #ahora las transformamos a angulares y despues aplicamos CD para mover al robot (INECESARIO PARA SIMULACION)
        #d_angulares=cinematica_inversa_omni(dx_locales[0][0],dx_locales[1][0],dx_locales[2][0],Robug.rb,Robug.rw)
        #dx=cinematica_directa_omni(Robug.rb,Robug.rw,d_angulares[0][0],d_angulares[1][0],d_angulares[2][0])
        
        #Movemos al robot #integrando las velocidades
        Robug.x=    integrador_euler(Robug.x,dt,dx_locales[0][0])
        Robug.y=    integrador_euler(Robug.y,dt,dx_locales[1][0])
        Robug.teta= integrador_euler(Robug.teta,dt,dx_locales[2][0])
def Control_P_OMNIDIRECCIONAL(trayectoria_deseada,d_trayec,Robug,k_x,k_y,k_teta,dt):
        Robug.vec_x.append(float(Robug.x))
        Robug.vec_y.append(float(Robug.y))
        Robug.vec_teta.append(float(Robug.teta))

        P_deseada=trayectoria_deseada
        #Calculo del error vectorial
        error_robug  = P_deseada - np.array([[Robug.x],
                                             [Robug.y],
                                             [Robug.teta]])
        matriz_k=np.array([[k_x, 0 ,0],
                           [0  ,k_y,0],
                           [ 0 , 0 ,k_teta]])

        velocidades_globales=d_trayec - matriz_k@(np.array([[Robug.x],
                                             [Robug.y],
                                             [Robug.teta]])-trayectoria_deseada)

        #ahora las transformamos a angulares y despues aplicamos CD para mover al robot (INECESARIO PARA SIMULACION)
        #Las transformamos a locales para el robot 
        dx_locales = transforma_global_a_local_omnidireccional(Robug,velocidades_globales)
        d_angulares=cinematica_inversa_omni(dx_locales[0][0],dx_locales[1][0],dx_locales[2][0],Robug.rb,Robug.rw)
        dx=cinematica_directa_omni(Robug.rb,Robug.rw,d_angulares[0][0],d_angulares[1][0],d_angulares[2][0])
        
        #Movemos al robot #integrando las velocidades
        Robug.x=    integrador_euler(Robug.x,dt,dx[0][0])
        Robug.y=    integrador_euler(Robug.y,dt,dx[1][0])
        Robug.teta= integrador_euler(Robug.teta,dt,dx[2][0])
        

#FUNCIONES DE CONTROL DE DRONES
def model_simplificado_drone_part1(g,mass,ixx,iyy,izz,u1,u2,u3,u4,d_z,d_roll,d_pitch,d_yaw):
    ds_1=d_z #velocidad en z
    ds_2=-g +(u1/mass) #aceleración en z
    ds_3=d_roll #velocidad angular en x
    ds_4=u2/ixx #aceleración angular en x
    ds_5=d_pitch #velocidad angular en x
    ds_6=u3/iyy #aceleración angular en x
    ds_7=d_yaw #velocidad angular en x
    ds_8=u4/izz #aceleración angular en x
    return(ds_1,ds_2,ds_3,ds_4,ds_5,ds_6,ds_7,ds_8)
def model_simplificado_drone_part2(mass,roll,pitch,yaw,u1):
    dd_x=(sin(yaw)*sin(roll) +cos(yaw)*sin(pitch)*cos(roll))*(u1/mass)
    dd_y=(-cos(yaw)*sin(roll) +sin(yaw)*sin(pitch)*cos(roll))*(u1/mass)
    return (dd_x,dd_y)
def Dinamica_Dron(dt,g,Drone1,U,d_x_0,d_y_0,d_z_0,d_roll_0,d_pitch_0,d_yaw_0):
    u1=U[0]
    u2=U[1]
    u3=U[2]
    u4=U[3]
    d_z,dd_z,d_r,dd_r,d_p,dd_p,d_yw,dd_yw=model_simplificado_drone_part1(g,Drone1.mass,Drone1.ixx,Drone1.iyy,Drone1.izz,u1,u2,u3,u4,d_z_0,d_roll_0,d_pitch_0,d_yaw_0)
    #integramos
    Drone1.z=integrador_euler(Drone1.z,dt,d_z)
    Drone1.roll=integrador_euler(Drone1.roll,dt,d_r)
    Drone1.pitch=integrador_euler(Drone1.pitch,dt,d_p)
    Drone1.yaw=integrador_euler(Drone1.yaw,dt,d_yw)
    d_z_0    =integrador_euler(d_z_0,dt,dd_z)
    d_roll_0 =integrador_euler(d_roll_0,dt,dd_r)
    d_pitch_0=integrador_euler(d_pitch_0,dt,dd_p)
    d_yaw_0  =integrador_euler(d_yaw_0,dt,dd_yw)
    ddx,ddy  =model_simplificado_drone_part2(Drone1.mass,Drone1.roll,Drone1.pitch,Drone1.yaw,u1)
    d_x_0=integrador_euler(d_x_0,dt,ddx)
    d_y_0=integrador_euler(d_y_0,dt,ddy)
    Drone1.x=integrador_euler(Drone1.x,dt,d_x_0)
    Drone1.y=integrador_euler(Drone1.y,dt,d_y_0)
    return(d_x_0,d_y_0,d_z_0,d_roll_0,d_pitch_0,d_yaw_0)
def control_z(Drone1,t_acumulado,dt,z_deseado,error_z_pasado,i_error_z):
    if t_acumulado>0:
        pass
    else:
        z_deseado=0
    error_z=z_deseado -Drone1.z
    if abs(error_z)<0.0001:
        error_z=0
    d_error_z=derivador(error_z_pasado , dt, error_z)
    i_error_z=integrador_euler(i_error_z,dt,error_z)
    kp=10;
    kd=4;
    ki=.2
    u1=11.78+kp*error_z +kd*d_error_z+ki*i_error_z
    #actualizacion de datos de control
    error_z_pasado=error_z
    return(u1,error_z_pasado,i_error_z)
def control_x(Drone1,t_acumulado,dt,x_deseado,error_x_pasado,i_error_x):
    #Comienza el control de posición X , Y ------------------- 
    if t_acumulado>0:
        pass#.5*cos(t_acumulado)#*cos(Drone1.yaw)
    else:
        x_deseado=Drone1.x
    error_x=x_deseado - (Drone1.x)
    if abs(error_x)<0.0001:
        error_x=0   
    d_error_x=derivador(error_x_pasado,dt,error_x)
    i_error_x=integrador_euler(i_error_x,dt,error_x)
    kp_x=4
    kd_x=1.5
    ki_x=1
    pitch_deseado=(kp_x*error_x +kd_x*d_error_x + ki_x*i_error_x)
    error_x_pasado=error_x 
    return(pitch_deseado,error_x_pasado,i_error_x)
def control_y(Drone1,t_acumulado,dt,y_deseado,error_y_pasado,i_error_y):
    if t_acumulado>0:
        pass#.5*sin(t_acumulado)
    else:
        y_deseado=Drone1.y
    error_y=y_deseado - Drone1.y
    if abs(error_y)<0.0001:
        error_y=0
    d_error_y=derivador(error_y_pasado,dt,error_y)
    i_error_y=integrador_euler(i_error_y,dt,error_y)
    kp_y=4
    kd_y=1.5
    ki_y=1
    roll_deseado=-(kp_y*error_y +kd_y*d_error_y + ki_y*i_error_y)
    error_y_pasado=error_y
    return(roll_deseado,error_y_pasado,i_error_y) 
def control_Yaw(Drone1,t_acumulado,dt,yaw_deseado,error_yaw_pasado,i_error_yaw):
    #Comienza control PID para YAW ------------------------------------
    error_yaw=yaw_deseado-Drone1.yaw
    if abs(error_yaw)<0.0001:
        error_yaw=0
    d_error_yaw=derivador(error_yaw_pasado,dt,error_yaw)
    i_error_yaw=integrador_euler(i_error_yaw,dt,error_yaw)
    error_yaw_pasado=error_yaw
    kp_yaw=1
    ki_yaw=0
    kd_yaw=0.49
    u4=kp_yaw*error_yaw +ki_yaw*i_error_yaw +kd_yaw*d_error_yaw
    return(u4,error_yaw_pasado,i_error_yaw)
def control_Roll(Drone1,t_acumulado,dt,roll_deseado,error_roll_pasado,i_error_roll):
    #PID para ROLL --------------------------------------
    error_roll=roll_deseado -Drone1.roll
    if abs(error_roll)<0.0001:
        error_roll=0    
    i_error_roll=integrador_euler(i_error_roll,dt,error_roll)
    d_error_roll=derivador(error_roll_pasado,dt,error_roll)
    kp_roll=1
    ki_roll=0
    kd_roll=.49
    u2=kp_roll*error_roll +ki_roll*i_error_roll  + kd_roll*d_error_roll
    error_roll_pasado=error_roll
    return(u2,error_roll_pasado,i_error_roll)
def control_Pitch(Drone1,t_acumulado,dt,pitch_deseado,error_pitch_pasado,i_error_pitch):
    #PID para PITCH --------------------------------------
    error_pitch=pitch_deseado -Drone1.pitch
    if abs(error_pitch)<0.0001:
        error_pitch=0

    i_error_pitch=integrador_euler(i_error_pitch,dt,error_pitch)
    d_error_pitch=derivador(error_pitch_pasado,dt,error_pitch)
    kp_pitch=1
    ki_pitch=0
    kd_pitch= .49
    u3=kp_pitch*error_pitch +ki_pitch*i_error_pitch  + kd_pitch*d_error_pitch
    error_pitch_pasado=error_pitch
    return(u3,error_pitch_pasado,i_error_pitch)
def Control_dinamico_drones(Drone1,dt,t_acumulado,tray,g,yaw_deseado):
    error_x_pasado=Drone1.error_x_pasado
    error_y_pasado=Drone1.error_y_pasado
    error_z_pasado=Drone1.error_z_pasado
    
    error_roll_pasado =Drone1.error_roll_pasado
    error_pitch_pasado=Drone1.error_pitch_pasado
    error_yaw_pasado  =Drone1.error_yaw_pasado
    
    i_error_x =Drone1.i_error_x
    i_error_y =Drone1.i_error_y
    i_error_z =Drone1.i_error_z
    
    i_error_roll  =Drone1.i_error_roll
    i_error_pitch =Drone1.i_error_pitch
    i_error_yaw   =Drone1.i_error_yaw
    
    d_x_0=Drone1.d_x_0
    d_y_0=Drone1.d_y_0
    d_z_0=Drone1.d_z_0
    
    d_roll_0 =Drone1.d_roll_0
    d_pitch_0=Drone1.d_pitch_0
    d_yaw_0  =Drone1.d_yaw_0
    
    x_deseado=tray[0];y_deseado=tray[1];z_deseado=tray[2];
    u1,error_z_pasado,i_error_z=control_z(Drone1,t_acumulado,dt,z_deseado,error_z_pasado,i_error_z)
    pitch_deseado,error_x_pasado,i_error_x=control_x(Drone1,t_acumulado,dt,x_deseado,error_x_pasado,i_error_x)
    roll_deseado,error_y_pasado,i_error_y =control_y(Drone1,t_acumulado,dt,y_deseado,error_y_pasado,i_error_y)
    #atacaremos el pitch deseado
    #print(roll_deseado,' t ',t_acumulado)

    u4,error_yaw_pasado,i_error_yaw=control_Yaw(Drone1,t_acumulado,dt,yaw_deseado,error_yaw_pasado,i_error_yaw)
    u2,error_roll_pasado,i_error_roll=control_Roll(Drone1,t_acumulado,dt,roll_deseado,error_roll_pasado,i_error_roll)
    u3,error_pitch_pasado,i_error_pitch=control_Pitch(Drone1,t_acumulado,dt,pitch_deseado,error_pitch_pasado,i_error_pitch)
    
    #CALCULO DE LA DINAMICA DEL DRON----------------------------------------------------------
    U=[u1,u2,u3,u4]
    #print(U,t_acumulado)
    d_x_0,d_y_0,d_z_0,d_roll_0,d_pitch_0,d_yaw_0=Dinamica_Dron(dt,g,Drone1,U,d_x_0,d_y_0,d_z_0,d_roll_0,d_pitch_0,d_yaw_0)
    
    Drone1.error_x_pasado=error_x_pasado
    Drone1.error_y_pasado=error_y_pasado
    Drone1.error_z_pasado=error_z_pasado

    Drone1.error_roll_pasado =error_roll_pasado
    Drone1.error_pitch_pasado=error_pitch_pasado
    Drone1.error_yaw_pasado  =error_yaw_pasado

    Drone1.i_error_x =i_error_x
    Drone1.i_error_y =i_error_y
    Drone1.i_error_z =i_error_z
    
    Drone1.i_error_roll  =i_error_roll
    Drone1.i_error_pitch =i_error_pitch
    Drone1.i_error_yaw   =i_error_yaw

    Drone1.d_x_0 =d_x_0
    Drone1.d_y_0 =d_y_0
    Drone1.d_z_0 =d_z_0
    
    Drone1.d_roll_0  =d_roll_0
    Drone1.d_pitch_0 =d_pitch_0
    Drone1.d_yaw_0   =d_yaw_0
    #Guardo valores en vectores de memoria
    Drone1.vec_x.append(float(Drone1.x))
    Drone1.vec_y.append(float(Drone1.y))
    Drone1.vec_z.append(float(Drone1.z))
    
    Drone1.vec_roll.append(float(Drone1.roll))
    Drone1.vec_pitch.append(float(Drone1.pitch))
    Drone1.vec_yaw.append(float(Drone1.yaw))
    return 0
   
    




    
    
#OBJETOS
class Pelota():
    def __init__(self, x,y,z,radio,color,dx,ddx,dy,ddy,mass) :
        self.x   = x
        self.y   = y
        self.z   = z        
        self.radio   = radio
        self.color =color
        self.dx =dx
        self.ddx=ddx 
        self.dy=dy
        self.ddy=ddy
        self.mass=mass 
class Robot_diff():
    def __init__(self, x,y,z,teta,radio,color,L) :
        self.x   = x
        self.y   = y
        self.z   = z        
        self.teta= teta
        self.radio   = radio
        self.color =color
        self.L=L #Distancia entre las ruedas
class Robot_Omni():
    def __init__(self, x,y,z,teta,color,rb,rw,error_x_pasado,error_y_pasado,error_teta_pasado,i_error_x,i_error_y,i_error_teta,vec_x,vec_y,vec_teta) :
        self.x   = x
        self.y   = y
        self.z   = z
        self.teta= teta
        self.color =color
        self.rb=rb #distancia del centro a la rueda
        self.rw= rw #radio de las ruedas
        #variables de memoria para el control
        self.error_x_pasado =error_x_pasado 
        self.error_y_pasado=error_y_pasado
        self.error_teta_pasado=error_teta_pasado
        self.i_error_x =i_error_x
        self.i_error_y =i_error_y
        self.i_error_teta =i_error_teta
        self.vec_x =vec_x
        self.vec_y =vec_y
        self.vec_teta =vec_teta

class Drone_Cuadri():
    def __init__(self, x,y,z,roll,pitch,yaw,color,radio,mass,ixx,iyy,izz,error_x_pasado,error_y_pasado,error_z_pasado,error_roll_pasado,error_pitch_pasado,error_yaw_pasado,i_error_x,i_error_y,i_error_z,i_error_roll,i_error_pitch,i_error_yaw,d_x_0,d_y_0,d_z_0,d_roll_0,d_pitch_0,d_yaw_0,vec_x,vec_y,vec_z,vec_roll,vec_pitch,vec_yaw) :
        #posicion
        self.x   = x
        self.y   = y
        self.z   = z
        #orientacion
        self.roll = roll
        self.pitch= pitch #radio del robot
        self.yaw  = yaw #radio del robot
        #Parametros
        self.color =color
        self.mass  =mass 
        self.radio = radio
        self.ixx   =ixx 
        self.iyy   =iyy 
        self.izz   =izz
        #Memoria errores
        self.error_x_pasado=error_x_pasado
        self.error_y_pasado=error_y_pasado
        self.error_z_pasado=error_z_pasado
        
        self.error_roll_pasado=error_roll_pasado
        self.error_pitch_pasado=error_pitch_pasado
        self.error_yaw_pasado=error_yaw_pasado

        self.i_error_x = i_error_x
        self.i_error_y = i_error_y
        self.i_error_z = i_error_z
        
        self.i_error_roll = i_error_roll
        self.i_error_pitch = i_error_pitch
        self.i_error_yaw = i_error_yaw

        #Velocidades Iniciales
        self.d_x_0= d_x_0
        self.d_y_0= d_y_0
        self.d_z_0= d_z_0
        
        self.d_roll_0= d_roll_0
        self.d_pitch_0= d_pitch_0
        self.d_yaw_0= d_yaw_0

        #vectores de memorias
        self.vec_x=vec_x 
        self.vec_y=vec_y 
        self.vec_z=vec_z 

        self.vec_roll =vec_roll 
        self.vec_pitch=vec_pitch 
        self.vec_yaw  =vec_yaw 


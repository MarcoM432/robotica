#Programa 3 
#Simulacion de un robot diferencial
from Libreria_Global import*
#Para navegar en el simulador
#usamos las flechas del teclado
#la tecla K para vista aerea

#la tecla L para vista 3D


#Variables del simulador
vec_data=[]
vec_t=[]
angulo_z=0
escala=1
inc=0
v_control_1=0
v_control_2=0
t_acumulado=0

#Inicia la simulación
pygame.init()
display = (700,700) #pixeles 
pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
gluPerspective(45, 1.2*(display[0]/display[1]), 0.1, 50.0) 
glTranslatef(0,-1.5, -12)  #Matriz de translacion
glRotatef(-90, 1, 0, 0)
glRotatef(-118-10, 0, 0, 1)
token=0

#Aqui vamos a ingresar nuestra logica o el codigo
t0=time.time()
#Creamos un Robot 
#Robot_1=Robot_diff(x,y,z,teta,radio,color,L)
color=(1,0,0)#RGB

Robot_1=Robot_diff(1,1,0,0,0.03,color,0.3)

    


while token==0:
    tf=time.time()
    dT=tf-t0 #tiempo real de la computadora
    dt=0.004 #Tiempo en paso fijo
    #dt=dT #paso variable
    t0=tf
    t_acumulado=t_acumulado+dt
    #eventos de tecladoo-------------------------------------------------------------  
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            token=1
        if event.type==pygame.KEYDOWN:
            if event.key==pygame.K_LEFT:
                angulo_z=-.1
            if event.key==pygame.K_RIGHT:
                angulo_z=.1
            if event.key == pygame.K_UP:
                inc=0.001
            if event.key==pygame.K_DOWN:
                inc=-0.001
            if event.key==pygame.K_a:
                v_control_1=-0.005
            if event.key==pygame.K_d:
                v_control_1=+0.005
            if event.key==pygame.K_w:
                v_control_2=+0.005
                Objeto.dx=1
                Objeto.dy=-2
                
            if event.key==pygame.K_s:
                v_control_2=-0.005
            if event.key==pygame.K_k:
                glRotatef(128, 0, 0, 1);glRotatef(90, 1, 0, 0 );glTranslatef(0,1.5, 0)  #Matriz de translacion
            if event.key==pygame.K_l:
                glTranslatef(0,-1.5, 0);glRotatef(-90, 1, 0, 0);glRotatef(-128, 0, 0, 1)
        #AL DEJAR DE PRECIONAR
        if event.type==pygame.KEYUP:            
            if event.key==pygame.K_LEFT:
                angulo_z=0
            if event.key==pygame.K_RIGHT:
                angulo_z=0
            if event.key==pygame.K_UP:
                inc=0
            if event.key==pygame.K_DOWN:
                inc=0
            if event.key==pygame.K_a:
                v_control_1=0
            if event.key==pygame.K_d:
                v_control_1=0
            if event.key==pygame.K_w:
                v_control_2=0
            if event.key==pygame.K_s:
                v_control_2=0

    #Aqui me encuentro dentro de la simulacion---------------------
    vec_t.append(t_acumulado)
                
                
    #COMIENZA DIBUJO---------------------------------------------------------        
    escala=escala+inc
    glRotatef(angulo_z, 0, 0, 1)
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    Draw_cartesiano(escala)

    Draw_robot_diferencial(Robot_1,escala)

   

    pygame.display.flip()
    #pygame.time.wait(5)
    if not dt==0:
        #print(1/dt)
        pass
#Aqui salimos del bucle WHILE
print("sali")
pygame.quit()
quit()



import numpy as np

#Colocar as condições iniciais:
def inicial():
    X = int(input("Coloque o tamanho de um lado da caixa"))
    Y = int(input("Coloque o tamanha do outro lada do caixa"))
    dt = float(input("Coloque o passo de tempo"))
    N = int(input("Coloque o numero de particulas"))
    tmax = float(input("Coloque o tempo final"))
    temp = float(input("Coloque a temperatura inicial"))
    l2 = float(input("Coloque a minima distancia para o calculo da forca"))

    #Distribuir uniformemente as particulas

    nx = np.sqrt(N)          #numeros de particulas em X
    ny = np.sqrt(N)          #numeros de particulas em Y

    #O numpy.tile repete a lista quantas vezes o programado
    #o numpy.repeat repete cada elemento quantas vezes o programado

    x = np.tile(np.linspace(0.5,X-0.5,nx),nx)
    y = np.repeat(np.linspace(0.5,Y-0.5,ny),ny)

    #O numpy.linspace cria uma lista aonde o programador da o primeiro numero, o ultimo e quantos numeros ele quer na lista

    vx = np.random(0,0.5,N)
    vy = np.random(0,0.5,N)

    return x,y,vx,vy,X,Y,l2,tmax,dt,N

x,y,vx,vy,X,Y,l2,tmax,dt,N = inicial() #nome da função

#fazer as forcas entre as particulas

def forcas(x,y,X,Y,l2):

    #inicializar as matrizes para as forcas e distancia
    N = len(x)
    f = np.zeros(shape=(N,N))
    v = np.zeros(shape=(N,N))

    #Calcular as distancias
    distx = x[:, np.newaxis] - x[np.newaxis, :]   #np.newaxis ?????
    disty = y[:, np.newaxis] - y[np.newaxis, :]

    #condicao periodica de contorno
    distx = distx - X * np.rint(distx / X) #O np.rint aproxima a divisao distx/X para o inteiro mais proximo.
    disty = disty - Y * np.rint(disty / Y)

    #Junta a distancia x com a y
    R2 = distx*distx + disty*disty

    k = (R2>0.0)&(R2<l2) #calcula o intervalo da distancia para que seja calculado a interacao entre as particulas

    v[k] = 4*((R2[k]**-6.)-(R2[k]**-3.))  #?????

    #calculo da forca

    f[k] = 48*(np.sqrt(R2[k])**(-13.) - 0.5*np.sqrt(R2[k])**(-7.))

    #separando as forcas nas componentes X e Y
    Fx = f*distx
    Fy = f*disty

    #fazer a soma para saber a soma em cada particula.
    v = np.sum(v,axis=1)
    fx = np.sum(Fx,axis=1)
    fy = np.sum(Fy,axis=1)

    return fx,fy,v,R2










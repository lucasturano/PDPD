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

#Calcular os indices das particulas

def verletlist(R2,rv,rc):
    s = (R2<rv*rv)&(R2<0.)
    vlist = [np.nonzero(s[i])[0] for i in range(len(R2))]

    return vlist

#Calcular as forcas com o metodo de verlet

def forcas_verlet(x,y,xold,yold,X,Y,l2,R2,rv,rc,vlist):

    desloc = np.sqrt((xold - x)**2 + (yold-y)**2)

    if (np.any(desloc)>(rv-rc)):
        p = np.c_[x,y]
        size = len(x)
        fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)
        #Calcular os indices das listas de verlet
        vlist = verletlist(R2,rv,rc)

        return fx,fy,V,R2,vlist

    else:
        f = np.zeros(shape=(N,N))
        V = np.zeros(shape=(N,N))
        Fx = np.zeros(shape=(N, N))
        Fy = np.zeros(shape=(N, N))
        for i in range(R2.shape[0]):

            distx = x[i] - x[vlist[i]]
            disty = y[i] - y[vlist[i]]

            #condicao periodica de contorno
            distx = distx - X*np.rint(distx/X)
            disty = disty - Y*np.rint(disty/Y)

            R2[i,:][vlist[i]] = distx**2 + disty**2
            r2 = distx**2 + disty**2

            V[i,:][vlist[i]] = 4*((r2**-6.)-(r2**-3.))
            f[i, :][vlist[i]] = 48*(np.sqrt(r2)**(-13.) - 0.5*np.sqrt(r2)**(-7.))
            #separando as forcas nos componestes x e y
            Fx[i, :][vlist[i]] = f[i,:][vlist[i]]*distx
            Fy[i,:][vlist[i]] = f[i,:][vlist[i]]*disty
            #Somando os resultados para ter a froca resultante em cada particula

        fx = np.sum(Fx,axis=1)
        fy = np.sum(Fy,axis=1)
        V = np.sum(V,axis=1)
        return fx,fy,V,R2,vlist

def integrate(x,y,vx,vy,fx,fy,dt):
    X = x.copy()
    Y = y.copy()
    ax = fx
    ay = fy
    vx += ax*dt
    vy += ay*dt
    X += vx*dt
    Y += vy*dt
    return X,Y,vx,vy

def period(x,y,X,Y):
    x[(x<=0.0)]= X+x[(x<=0.0)]
    y[(y<=0.0)]= Y+y[(y<=0.0)]
    x[(x>=X)]= x[(x>=X)] - X
    y[(y>=Y)]= y[(y>=Y)] - Y
    return x,y

def range_of_steps(x,y,vx,vy,X,Y,l2,dt,s):
    for i in range(s):
        fx,fy,V,R2 = forcas(x,y,X,Y,l2)
        x,y,vx,vy = integrate(x,y,vx,vy,fx,fy,dt)
    return x,y,vx,vy,fx,fy,V,R2


def termo_andersen(vx, vy, dt, Temp, nu):
    if (np.random.uniform(0, 1) < nu * dt):
        sig = np.sqrt(Temp)
        vx = np.random.normal(0, sig, len(vx))
        vy = np.random.normal(0, sig, len(vy))

    return vx, vy











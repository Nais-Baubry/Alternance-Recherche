import numpy as np
from scipy import *
from math import *


def IX(x,y):
    f=x+y*N
    return f



###STRUCTURE CUBE###

#size
N=10
#lenght of the timestep
dt=1
#amount of diffusion
diff=1
#viscosity
visc=1
#density
density=np.zeros(N*N)
#Velocity(x,y)
Vx=np.zeros(N*N)
Vy=np.zeros(N*N)

#Initial velocity (x,y)
Vx0=np.zeros(N*N)
Vy0=np.zeros(N*N)


#Add dye
def FluidCubeAddDensity(x,y,amount):
    density[IX(x, y)] += amount
    
    
    
def FluidCubeAddVelocity(x,y,amountX,amountY):
    index = IX(x, y)
    Vx[index] += amountX
    Vy[index] += amountY

#limite les bords extérieurs des cubes

def set_bnd(b, x, N):  
    for k in range(1,N-1):
        for i in range(N):
            x[IX(i, 0 )] = b #== 2 ?
            x[IX(i, N-1)] = b #== 2 ? 
    
    for k in range(1,N-1):
        for j in range(N):
            x[IX(0  , j)] = b #== 1 ? 
            x[IX(N-1, j)] = b #== 1 ? 
            
        
    x[IX(0, 0)]       = 0.33 * (x[IX(1, 0)] + x[IX(0, 1)]+ x[IX(0, 0)])
    x[IX(0, N-1)]     = 0.33 * (x[IX(1, N-1)] + x[IX(0, N-2)] + x[IX(0, N-1)])
    x[IX(N-1,0)]     = 0.33 * (x[IX(N-2, 0)] + x[IX(N-1, 1)] + x[IX(N-1, 0)])
    x[IX(N-1, N-1)]   = 0.33 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)] + x[IX(N-1, N-1)])


#resolve
    
def lin_solve(b, x, x0, a, c, iter, N):
    cRecip = 1.0 / c
    for k in range(iter):
        for j in range(N-1):
            for i in range(1,N-1):
                x[IX(i, j)] =(x0[IX(i, j)] + a*(x[IX(i+1, j)] + x[IX(i-1, j)] + x[IX(i, j+1)] + x[IX(i, j-1 )] )) * cRecip
    set_bnd(b, x, N)
    
    
def diffuse (b, x, x0, diff, dt, iter, N):
    a = dt * diff * (N - 2) * (N - 2)
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N)
    

def project(velocX, velocY, p, div, iter, N):
    div=np.zeros(N*N)
    for j in range(1,N-1):
        for i in range(1,N-1):
            div[IX(i, j)] = -0.5*(velocX[IX(i+1, j)] - velocX[IX(i-1, j)] + velocY[IX(i, j+1)] -velocY[IX(i, j-1)])/N
            p[IX(i, j)] = 0
    
    set_bnd(0, div, N) 
    set_bnd(0, p, N)
    lin_solve(0, p, div, 1, 6, iter, N)
    
    for j in range(1,N-1):
        for i in range(1,N-1):
            velocX[IX(i, j)] -= 0.5* (  p[IX(i+1, j)]-p[IX(i-1, j)]) * N
            velocY[IX(i, j)] -= 0.5* (  p[IX(i, j+1)]-p[IX(i, j-1)]) * N
              

    set_bnd(1, velocX, N)
    set_bnd(2, velocY, N)

#Moving
    
def advect(b, d, d0, velocX, velocY, dt, N):
    i0, i1, j0, j1=0,0,0,0
    dtx = dt * (N - 2)
    dty = dt * (N - 2)

    s0, s1, t0, t1=0,0,0,0
    tmp1, tmp2, x, y=0,0,0,0
    
    Nfloat = N

    for i in range(1,N-1):
        for j in range(1,N-1):
            jfloat=j
            ifloat=i
            tmp1 = dtx * velocX[IX(i, j)]
            tmp2 = dty * velocY[IX(i, j)]
            x    = ifloat - tmp1; 
            y    = jfloat - tmp2;
                
            if x < 0.5 :
                x = 0.5 
            if x > Nfloat + 0.5
                x = Nfloat + 0.5
            i0 = floor(x) 
            i1 = i0 + 1
            if y < 0.5 :
                y = 0.5 
            if y > Nfloat + 0.5 :
                y = Nfloat + 0.5 
            j0 = floor(y)
            j1 = j0 + 1

                
            s1 = x - i0 
            s0 = 1.0 - s1 
            t1 = y - j0 
            t0 = 1.0 - t1

                
            i0i = i0
            i1i = i1
            j0i = j0
            j1i = j1

                
            d[IX(i, j)] = s0 * ( t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)] ) + s1 * ( t0 * d0[IX(i1i, j0i)] +t1 * d0[IX(i1i, j1i)])

    set_bnd(b, d, N)
    
    
diffuse(1, Vx0, Vx, visc, dt, 4, N)
diffuse(2, Vy0, Vy, visc, dt, 4, N)
    
project(Vx0, Vy0, Vx, Vy, 4, N)
    
advect(1, Vx, Vx0, Vx0, Vy0, dt, N)
advect(2, Vy, Vy0, Vx0, Vy0, dt, N)

    
project(Vx, Vy, Vx0, Vy0, 4, N)
    
diffuse(0, s, density, diff, dt, 4, N)
advect(0, density, s, Vx, Vy, dt, N)


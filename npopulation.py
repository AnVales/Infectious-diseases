# SIR

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1
beta=0.25
gamma=0.1
'''time step'''
TS=1 # time steep
ND=400. # tiempo final
S0=N-N*0.0010 # susceptibles iniciales
I0=N*0.0010 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] / N# V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] / N - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] # V[1] es I0, los infectados, dR/dt = gamma * I
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0], '-b', label='Susceptible')
pl.plot(RES[:,2], '-g', label='Recovered with immunity')
pl.plot(RES[:,1], '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIR epidemic with N = 1 and I0=0,1% of N')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIR

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=100
beta=0.20
gamma=0.005
'''time step'''
TS=1 # time steep
ND=400. # tiempo final
S0=N-N*1e-6 # susceptibles iniciales
I0=N*1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] / N # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] / N - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] # V[1] es I0, los infectados, dR/dt = gamma * I
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0], '-b', label='Susceptible')
pl.plot(RES[:,2], '-g', label='Recovered with immunity')
pl.plot(RES[:,1], '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIR epidemic with N = 100')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIR

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=10000
beta=0.20
gamma=0.005
'''time step'''
TS=1 # time steep
ND=400. # tiempo final
S0=N-N*1e-6 # susceptibles iniciales
I0=N*1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] / N # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] / N - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] # V[1] es I0, los infectados, dR/dt = gamma * I
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0], '-b', label='Susceptible')
pl.plot(RES[:,2], '-g', label='Recovered with immunity')
pl.plot(RES[:,1], '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIR epidemic with N = 10000')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();


# SIR

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1000000
beta=0.20
gamma=0.005
'''time step'''
TS=1 # time steep
ND=400. # tiempo final
S0=N-N*1e-6 # susceptibles iniciales
I0=N*1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] / N# V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] / N  - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] # V[1] es I0, los infectados, dR/dt = gamma * I
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0], '-b', label='Susceptible')
pl.plot(RES[:,2], '-g', label='Recovered with immunity')
pl.plot(RES[:,1], '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIR epidemic with N = 1000000')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();
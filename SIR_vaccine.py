# SIR VACUNA

import scipy.integrate as spi
import numpy as np
import pylab as pl

# SIR #

'''tamaño poblacional'''
N=1
beta = 0.75
gamma=0.0075
'''time step'''
TS=1 # time steep
ND=400. # tiempo final
S0=1-1e-6 # susceptibles iniciales
I0=1e-6 # infectados iniciales
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
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIR epidemic with beta=0.75 & mu=0.0075\nwithout vaccine')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();


# SIR VACUNA#


'''tamaño poblacional'''
N=1
beta=0.75
gamma=0.0075
p=0.05
'''time step'''
TS=1 # time steep
ND=400.0 # tiempo final
S0=1-1e-6 # susceptibles iniciales
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
V0=0.0 # vacunados
INPUT = (S0, I0, R0, V0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((4)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] - p * V[0] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] # V[1] es I0, los infectados, dR/dt = gamma * I
    Y[3] = p * V[0]
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.plot(RES[:,3]*N, '-k', label='Vaccinated population')
pl.legend(loc='best')
pl.title('SIR epidemic with vaccine with p = 0.05')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIR VACUNA#


'''tamaño poblacional'''
N=1
beta=0.75
gamma=0.0075
p=0.01
'''time step'''
TS=1 # time steep
ND=400.0 # tiempo final
S0=1-1e-6 # susceptibles iniciales
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
V0=0.0 # vacunados
INPUT = (S0, I0, R0, V0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((4)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] - p * V[0] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] # V[1] es I0, los infectados, dR/dt = gamma * I
    Y[3] = p * V[0]
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.plot(RES[:,3]*N, '-k', label='Vaccinated population')
pl.legend(loc='best')
pl.title('SIR epidemic with vaccine with p = 0.01')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# VACINNE #

#params
beta=0.75
gamma=0.0075

# p = 0.01

dt = 1            # 6 min
D = 400              # simulate for D days
N = int(D/dt)     # corresponding no of DAYS

from numpy import zeros, linspace
t = linspace(0, N*dt, N+1)
S = zeros(N+1)
V = zeros(N+1)
I = zeros(N+1)
R = zeros(N+1)

# Initial condition
S[0] = 1-1e-6
V[0] = 0
I[0] = 1e-6
R[0] = 0

for n in range(N):

    if n<10:
        pn=0
    else:
        pn=0.05

    S[n+1] = S[n] - dt*beta*S[n]*I[n] - dt*S[n]*pn
    V[n+1] = V[n] + dt*S[n]*pn
    I[n+1] = I[n] + dt*beta*S[n]*I[n] - dt*gamma*I[n]
    R[n+1] = R[n] + dt*gamma*I[n] 

from matplotlib.pyplot import plot, savefig, legend, xlabel, show, ylabel
plot(t, S, '-b', t, I, '-r', t, R, '-g', t, V, '-k')
legend(['Susceptible', 'Infectious', 'Recovered with immunity', 'Vaccinated population'], loc='best')
xlabel('Time')
ylabel('Population')
pl.title('SIR epidemic with vaccine with p = 0.05 since the day 10')

show();

# VACINNE #

#params
beta=0.75
gamma=0.0075

# p = 0.01

dt = 1            # 6 min
D = 400              # simulate for D days
N = int(D/dt)     # corresponding no of DAYS

from numpy import zeros, linspace
t = linspace(0, N*dt, N+1)
S = zeros(N+1)
V = zeros(N+1)
I = zeros(N+1)
R = zeros(N+1)

# Initial condition
S[0] = 1-1e-6
V[0] = 0
I[0] = 1e-6
R[0] = 0

for n in range(N):

    if S[n]<=0:
        pn=0
    else:
        pn=0.01*(S[0]+I[0])

    S[n+1] = S[n] - dt*beta*S[n]*I[n] - pn
    V[n+1] = V[n] + pn
    I[n+1] = I[n] + dt*beta*S[n]*I[n] - dt*gamma*I[n]
    R[n+1] = R[n] + dt*gamma*I[n] 

from matplotlib.pyplot import plot, savefig, legend, xlabel, show, ylabel
plot(t, S, '-b', t, I, '-r', t, R, '-g', t, V, '-k')
legend(['Susceptible', 'Infectious', 'Recovered with immunity', 'Vaccinated population'], loc='best')
xlabel('Time')
ylabel('Population')
pl.title('SIR epidemic with a percentage of 0.01%\nof the initial population vaccinated daily')

show();
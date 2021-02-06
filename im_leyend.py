print('I AM LEYEND')

# A genetically re-engineered measles virus, originally created as a cure for cancer, turns into a lethal strain
# which kills 90% of those it infects, while 9% mutate into predatory, nocturnal and vampiric mutants able
# to transmit the disease. Only 1% of the population is immune. 

# Mutants are immortal but unable to reproduce. When they meet a healthy human, there are three possible
# outcomes: either they are killed by the human, the human is killed (probability above) or the disease is transmitted. 

# Healthy humans can reproduce; newborns are virus-free. Humans die also of other diseases and accidents. 

# Given a sensible set of parameters, estimate the rates of human reproduction and mutant killing that lead to the
# disappearance of the latter. Evaluate the human and mutant death toll.

import scipy.integrate as spi
import numpy as np
import pylab as pl

# more mortality than birth rate 

beta=0.035
birth=0.001
dead=0.002
kill=0.0001

'''time step'''
TS=1 # time steep
ND=1800. # tiempo final
S0=0.97
I0=0.02
M0=0.01
D0=0
K0=0

N= S0 + I0 + M0 + D0

INPUT = (S0, I0, M0, D0, K0) 

def diff_eqs(INP,t):  
    Y=np.zeros((5)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = birth * (V[0] + V[2]) * 0.99 - beta * V[0] * V[1] - dead * V[0]
    Y[1] = beta * V[0] * V[1] * (0.09 / 0.99) - kill * V[1] * (V[0] + V[2])
    Y[2] = birth * (V[0] + V[2]) * 0.01 - dead * V[2]
    Y[3] = beta * V[0] * V[1] *(0.90 / 0.90) + dead * V[0] + dead * V[2]
    Y[4] = kill * V[1] * (V[0] + V[2])
    return Y   # For odeint

# TIME

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.plot(RES[:,3]*N, '-k', label='Dead human')
pl.plot(RES[:,4]*N, '-y', label='Dead mutant')
pl.legend(loc='best')
pl.title('I am legend epidemic with more mortality\n than birth rate')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();


# less mortality than birth rate #

beta=0.035
birth=0.002
dead=0.001
kill=0.0001

'''time step'''
TS=1 # time steep
ND=1800. # tiempo final
S0=0.97
I0=0.02
M0=0.01
D0=0
K0=0

N= S0 + I0 + M0 + D0

INPUT = (S0, I0, M0, D0, K0) 

def diff_eqs(INP,t):  
    Y=np.zeros((5)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = birth * (V[0] + V[2]) * 0.99 - beta * V[0] * V[1] - dead * V[0]
    Y[1] = beta * V[0] * V[1] * (0.09 / 0.99) - kill * V[1] * (V[0] + V[2])
    Y[2] = birth * (V[0] + V[2]) * 0.01 - dead * V[2]
    Y[3] = beta * V[0] * V[1] *(0.90 / 0.90) + dead * V[0] + dead * V[2]
    Y[4] = kill * V[1] * (V[0] + V[2])
    return Y   # For odeint

# TIME

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.plot(RES[:,3]*N, '-k', label='Dead human')
pl.plot(RES[:,4]*N, '-y', label='Dead mutant')
pl.legend(loc='best')
pl.title('I am legend epidemic with less mortality\n than birth rate')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# mismo ratio #

beta=0.035
birth=0.001
dead=0.001
kill=0.0001

'''time step'''
TS=1 # time steep
ND=1800. # tiempo final
S0=0.97
I0=0.02
M0=0.01
D0=0
K0=0

N= S0 + I0 + M0 + D0

INPUT = (S0, I0, M0, D0, K0) 

def diff_eqs(INP,t):  
    Y=np.zeros((5)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = birth * (V[0] + V[2]) * 0.99 - beta * V[0] * V[1] - dead * V[0]
    Y[1] = beta * V[0] * V[1] * (0.09 / 0.99) - kill * V[1] * (V[0] + V[2])
    Y[2] = birth * (V[0] + V[2]) * 0.01 - dead * V[2]
    Y[3] = beta * V[0] * V[1] *(0.90 / 0.90) + dead * V[0] + dead * V[2]
    Y[4] = kill * V[1] * (V[0] + V[2])
    return Y   # For odeint

# TIME

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.plot(RES[:,3]*N, '-k', label='Dead human')
pl.plot(RES[:,4]*N, '-y', label='Dead mutant')
pl.legend(loc='best')
pl.title('I am legend epidemic with the same mortality\n as birth rate')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();


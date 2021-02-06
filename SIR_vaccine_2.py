# SIR #

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1
beta=1.2
gamma=0.20
'''time step'''
TS=1 # time steep
ND=100. # tiempo final
S0=1-1e-6 # susceptibles iniciales
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
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
pl.title('SIR epidemic')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIR VACUNA

'''tamaño poblacional'''
N=1
beta=1.2
gamma=0.20
p=0.02
'''time step'''
TS=0.5 # time steep
ND=50.0 # tiempo final
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
pl.title('SIR epidemic with vaccine')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIRS #
import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1
beta=1.2
gamma=0.20
alfa = 0.025
'''time step'''
TS=0.5 # time steep
ND=50.0 # tiempo final
S0=1-1e-6 # susceptibles iniciales
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[1] + alfa * V[2]# V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I +  alfa * R
    Y[1] = beta * V[0] * V[1] - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] - alfa * V[2]# V[1] es I0, los infectados, dR/dt = gamma * I - alfa * R
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIRS epidemic')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIRS VACUNA

'''tamaño poblacional'''
N=1
beta=1.2
gamma=0.20
p=0.02
alfa = 0.025
'''time step'''
TS=0.5 # time steep
ND=50.0 # tiempo final
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
    Y[0] = - beta * V[0] * V[1] - p * V[0] + alfa * V[2] + alfa * V[3] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] - gamma * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] - alfa * V[2] # V[1] es I0, los infectados, dR/dt = gamma * I
    Y[3] = p * V[0] - alfa * V[3]
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
pl.title('SIRS epidemic with vaccine')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SEIR #

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1
beta=1.2
gamma=0.2
alfa=0.2
'''time step'''
TS=0.5 # time steep
ND=100.0 # tiempo final
S0=1-1e-6 # susceptibles iniciales
E0=0.0 # expuestos
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, E0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((4)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[2] 
    Y[1] = beta * V[0] * V[2] - alfa * V[1] 
    Y[2] = alfa * V[1] - gamma * V[2] 
    Y[3] = gamma * V[2]
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,3]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,2]*N, '-r', label='Infectious')
pl.plot(RES[:,1]*N, '-y', label='Exposed')
pl.legend(loc='best')
pl.title('SEIR epidemic')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SEIR WITH VACCINE #

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1
beta=1.2
gamma=0.2
alfa=0.2
p=0.02
'''time step'''
TS=0.5 # time steep
ND=200.0 # tiempo final
S0=1-1e-6 # susceptibles iniciales
E0=0.0 # expuestos
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
V0=0.0 #vacunados
INPUT = (S0, E0, I0, R0, V0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((5)) # el vector de 5 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = - beta * V[0] * V[2] - p * V[0]
    Y[1] = beta * V[0] * V[2] - alfa * V[1] 
    Y[2] = alfa * V[1] - gamma * V[2] 
    Y[3] = gamma * V[2]
    Y[4] = p * V[0]
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,3]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,2]*N, '-r', label='Infectious')
pl.plot(RES[:,1]*N, '-y', label='Exposed')
pl.plot(RES[:,4]*N, '-k', label='Vaccinated population')
pl.legend(loc='best')
pl.title('SEIR epidemic with vaccine')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIR DEMOGRAFÍA CON OSCILACIONES #

import scipy.integrate as spi
import numpy as np
import pylab as pl
#parámetros
mu=1/(70*365.0)
beta=1.2
gamma=0.14
TS=1.0
ND=70*365
#condiciones iniciales
S0=0.1
I0=1e-4
R0=1-S0-I0
INPUT = (S0, I0, R0)
#ecuaciones diferenciales
def diff_eqs(INP,t):  
    '''The main set of equations'''
    Y=np.zeros((3))
    V = INP    
    Y[0] = mu - beta * V[0] * V[1] - mu * V[0]
    Y[1] = beta * V[0] * V[1] - gamma * V[1] - mu * V[1]
    Y[2] = gamma * V[1] - mu * V[2]
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

#Ploting
pl.subplot(311)
pl.plot(RES[:,0], '-b', label='Susceptible')
pl.title('SIR bith and death')
pl.ylabel('Susceptible')
pl.subplot(312)
pl.plot(RES[:,1], '-r', label='Infectious')
pl.ylabel('Infectious')
pl.subplot(313)
pl.plot(RES[:,2], '-g', label='Recovered with immunity')
pl.xlabel('Tiempo')
pl.ylabel('Recovered \nwith immunity')
pl.show();

# SIR CON DEMOGRAFÍA #

import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1
beta=1.2
gamma=0.20
mu=1/(70*365.0)
'''time step'''
TS=0.5 # time steep
ND=50.0 # tiempo final
S0=1-1e-6 # susceptibles iniciales
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = mu - beta * V[0] * V[1] - mu * V[0] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I
    Y[1] = beta * V[0] * V[1] - gamma * V[1] - mu * V[1] # V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] - mu * V[2] # V[1] es I0, los infectados, dR/dt = gamma * I
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIR epidemic with demografy')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIRS DEMOGRAFÍA CON OSCILACIONES #

import scipy.integrate as spi
import numpy as np
import pylab as pl
#parámetros
alfa=0.0001
mu=1/(70*365.0)
beta=1.2
gamma=0.14
TS=1.0
ND=70*365
#condiciones iniciales
S0=0.1
I0=1e-4
R0=1-S0-I0
INPUT = (S0, I0, R0)
#ecuaciones diferenciales
def diff_eqs(INP,t):  
    '''The main set of equations'''
    Y=np.zeros((3))
    V = INP    
    Y[0] = mu - beta * V[0] * V[1] - mu * V[0] + alfa * V[2]
    Y[1] = beta * V[0] * V[1] - gamma * V[1] - mu * V[1]
    Y[2] = gamma * V[1] - mu * V[2] - alfa * V[2]
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

#Ploting
pl.subplot(311)
pl.plot(RES[:,0], '-b', label='Susceptible')
pl.title('SIRS bith and death')
pl.ylabel('Susceptible')
pl.subplot(312)
pl.plot(RES[:,1], '-r', label='Infectious')
pl.ylabel('Infectious')
pl.subplot(313)
pl.plot(RES[:,2], '-g', label='Recovered with immunity')
pl.xlabel('Tiempo')
pl.ylabel('Recovered \nwith immunity')
pl.show();

# SIRS DEMOGRAFÍA #
import scipy.integrate as spi
import numpy as np
import pylab as pl
'''tamaño poblacional'''
N=1
alfa=0.0001
mu=1/(70*365.0)
beta=1.2
gamma=0.14
'''time step'''
TS=0.5 # time steep
ND=50.0 # tiempo final
S0=1-1e-6 # susceptibles iniciales
I0=1e-6 # infectados iniciales
R0=0.0 # recuperados
INPUT = (S0, I0, R0) # input de susceptible, infectado y recuperado


def diff_eqs(INP,t):  
    Y=np.zeros((3)) # el vector de 3 zeros para las respuestas
    V = INP # vector V
    '''Las ecuaciones diferenciales'''
    V = INP # condiciones iniciales
    Y[0] = mu - beta * V[0] * V[1] + alfa * V[2] - mu * V[0]# V[0] es S0, los susceptibles, V[1] es I0, los infectados; dS/dt = beta * S * I +  alfa * R
    Y[1] = beta * V[0] * V[1] - gamma * V[1] - mu * V[1]# V[0] es S0, los susceptibles, V[1] es I0, los infectados; dI/dS = beta * S * I - gamma * I
    Y[2] = gamma * V[1] - alfa * V[2] - mu * V[2]# V[1] es I0, los infectados, dR/dt = gamma * I - alfa * R
    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc) # Arange: Return evenly spaced values within a given interval, [start, stop, steep) 
RES = spi.odeint(diff_eqs,INPUT,t_range) # Integrate a system of ordinary differential equations.

#Gráfica
pl.plot(RES[:,0]*N, '-b', label='Susceptible')
pl.plot(RES[:,2]*N, '-g', label='Recovered with immunity')
pl.plot(RES[:,1]*N, '-r', label='Infectious')
pl.legend(loc='best')
pl.title('SIRS epidemic')
pl.xlabel('Time')
pl.ylabel('Population')
# pl.savefig('sirpy')
pl.show();

# SIER DEMOGRAFÍA CON OSCILACIONES #

import scipy.integrate as spi
import numpy as np
import pylab as pl
#parámetros
alfa=0.0001
mu=1/(70*365.0)
beta=1.2
gamma=0.14
alfa=0.2
TS=1.0
ND=70*365
#condiciones iniciales
S0=0.1
I0=1e-4
E0=1e-4
R0=1-S0-I0-E0
INPUT = (S0,E0, I0, R0)
#ecuaciones diferenciales
def diff_eqs(INP,t):  
    '''The main set of equations'''
    Y=np.zeros((4))
    V = INP    
    Y[0] = mu - beta * V[0] * V[2] - mu * V[0]
    Y[1] = beta * V[0] * V[2] - alfa * V[1] - mu * V[1]
    Y[2] = alfa * V[1] - gamma * V[2] - mu * V[2]
    Y[3] = gamma * V[2] - mu * V[3]

    return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

#Ploting
pl.subplot(411)
pl.plot(RES[:,0], '-b', label='Susceptible')
pl.title('SIER bith and death')
pl.ylabel('Susceptible')
pl.subplot(412)
pl.plot(RES[:,1], '-y', label='Exposed')
pl.ylabel('Exposed')
pl.subplot(413)
pl.plot(RES[:,2], '-r', label='Infectious')
pl.ylabel('Infectious')
pl.subplot(414)
pl.plot(RES[:,3], '-g', label='Recovered with immunity')
pl.xlabel('Tiempo')
pl.ylabel('Recovered \nwith immunity')
pl.show();
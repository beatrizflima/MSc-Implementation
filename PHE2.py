import pandas as pd
from pandas import Series as s
import matplotlib.pyplot as plt
import seaborn as sns
from cmath import pi, log, sqrt
from FunctionsUFH import Reynolds, pressuredrop
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import pygfunction as gt
from CoolProp.CoolProp import PropsSI
from interpolateAMRcold import T_ci

def reynolds_plate(m_dot,rho,L,N,mu,D_h):
    u_m = m_dot/(rho*L**2/2)
    G = m_dot / ((N-1)*(gap*L/3))
    Re_D = G*D_h/mu
    return Re_D

def equations(p,T_hi,T_co,T_ln_in,T_ci):
    x = p
    return np.divide(T_hi-T_co-(x-T_ci),np.log((T_hi-T_co)/(x-T_ci)))-T_ln_in

D_i = 0.026*2
m = 0.3
rough = 1e-6
mu = 0.002
N=2
L = 120*2
fluid = gt.media.Fluid('MPG',8)
cp_f = fluid.cp
rho = fluid.rho
v = m/rho
vel = v/(D_i**2*pi/4)
Re = Reynolds(m,mu,D_i)
f = 0.316/Re**0.25
#f=64/Re
p = pressuredrop(f,L,D_i,rho,v)
h_w = 7.54*0.582
pump = p*N*m*N/(rho*0.5)

print(rho,Re,f,p,vel,pump)

outputs_BHE = pd.read_csv(r'C:\Users\beali\Documents\TUDelft\pygfunction-master\T_Q_outputs.csv')  #reads values obtained from the g-function model by Massimo Cimmino
Date = outputs_BHE.Time
T_initial = outputs_BHE.T_inlet
T_hi = outputs_BHE.T_outlet
q = outputs_BHE.Q_total
T_co = 8
N = 70
gap = 0.0015
max_mass_MCHP = 0.96655
density_HTF = PropsSI('D','T',8+273.15,'P',101325,'INCOMP::MEG-20%')
mu_HTF = PropsSI('V','T',8+273.15,'P',101325,'INCOMP::MEG-20%')
lambda_HTF = PropsSI('L','T',8+273.15,'P',101325,'INCOMP::MEG-20%')
h_htf = 7.54 * lambda_HTF
#----------set empty list----------#
T_ho = [0]*len(q)
T_ln = [0]*len(q)
ratio = [0]*len(q)
T_ln_needed = [0]*len(q)
m_HE = [0]*len(q)
m_valve = [0]*len(q)
m_total = [0]*len(q)
h_ho = [0]*len(q)
h_hi = [0]*len(q)
h_wi = [0]*len(q)
T_wi = [0]*len(q)
h_should = [0]*len(q)
Re_HP = [0]*len(q)

#-----------PHE geometry----------#

for i in range(len(q)):

    T_ln[i] = (T_hi[i] - T_co - (T_initial[i] - T_ci[i]))/log((T_hi[i] - T_co)/(T_initial[i] - T_ci[i])).real

    ratio[i] = q[i]/T_ln[i]
for i in range(len(q)):
    if  ratio[i]==max(ratio):
        print(i,q[i],T_ln[i])
        print('temperature inlet BHE',T_hi[i])  
        print('temperature outlet BHE',T_initial[i])  

print('the max ratio is: ',max(ratio))
UA = max(ratio)
L = (-366*N*gap*h_w*h_htf+366*gap*h_w*h_htf+sqrt((366*N*gap*h_w*h_htf-366*gap*h_w*h_htf)**2-4*(122*N*h_w*h_htf-122*h_w*h_htf)*(-600*gap*UA*h_w-600*gap*UA*h_htf)))/(2*(122*N*h_w*h_htf-122*h_w*h_htf))
D_h = 4 * (gap*L/3)/(2*(gap+L/3))
print('max Re MCHP',reynolds_plate(max_mass_MCHP,density_HTF,L,N,mu_HTF,D_h))

for i in range(len(q)):
    
    T_ln_needed[i] = q[i] / max(ratio)
    T_ho[i] = fsolve(equations,T_ci[i]+0.0001,args=(T_hi[i],T_co,T_ln_needed[i],T_ci[i]))
    m_HE[i] = q[i]/(cp_f*(T_hi[i]-T_ho[i]))
    m_valve[i] = m*2 - m_HE[i]
    h_ho[i] = PropsSI('H','T',T_ho[i]+273.15,'P',101325,'INCOMP::MEG-60%')
    h_hi[i] = PropsSI('H','T',T_hi[i]+273.15,'P',101325,'INCOMP::MEG-60%')
    h_wi[i] = (m_valve[i]*h_hi[i] + m_HE[i]*h_ho[i])/(m*2)
    h_should[i] = PropsSI('H','T',T_initial[i]+273.15,'P',101325,'INCOMP::MEG-60%')
    T_wi[i] = PropsSI('T','H',h_wi[i],'P',101325,'INCOMP::MEG-60%')-273.15
    Re_HP[i] = reynolds_plate(m_HE[i],rho,L,N,mu,D_h)

print('max Re',max(Re_HP))
print('L',L)
print('D',D_h)


#----polyfit----#
'''
q = np.array(q)
q=q.astype('float64')

m_HE = np.array(m_HE)
m_HE = m_HE.astype('float64')

m_valve = np.array(m_valve)
m_valve = m_valve.astype('float64')

a1,b1,c1,d1,e1 = np.polyfit(q,m_HE,4)
a2,b2,c2,d2,e2 = np.polyfit(q,m_valve,4)

h_ho = np.array(h_ho)
h_ho = h_ho.astype('float64')
h_hi = np.array(h_hi)
h_hi = h_hi.astype('float64')
h_wi = np.array(h_wi)
h_wi = h_wi.astype('float64')

a3,b3,c3 = np.polyfit(q,h_ho,2)
a4,b4,c4 = np.polyfit(q,h_hi,2)
a5,b5,c5 = np.polyfit(q,h_wi,2)
'''

#----graphs----#

plt.figure()
plt.plot(T_ln_needed)
plt.plot(T_ln)
plt.figure()

#plt.plot(T_wi,c='pink',label='Temperatura que entra no borehole')
plt.plot(T_hi, c='green',label='Temperatura que entra do PHE')

plt.plot(T_ho, c='green',label='Temperatura que sai do PHE')

plt.figure()
plt.scatter(q,m_HE,s=1,c='red',label='Mass flow rate through PHE2')
plt.scatter(q,m_valve,s=1,c='cyan',label='Mass flow rate directly to T2')
#plt.plot(q, a1*q**4 + b1*q**3 + c1*q**2 + d1*q + e1)
#plt.plot(q, a2*q**4 + b2*q**3 + c2*q**2 + d2*q + e2)
plt.xlabel('Heating Power [W]', fontsize=24)
plt.ylabel('Mass Flow Rate [kg/s]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(markerscale=6,fontsize=15)

plt.figure()
plt.scatter(T_wi,m_HE,s=1,c='red',label='HE')
plt.scatter(T_wi,m_valve,s=1,c='pink',label='valve')
plt.figure()
plt.plot(q, T_ho)
plt.figure()

plt.scatter(Date, T_hi, s=1,c='orange',label='HTF leaving BHE')
plt.scatter(Date,T_wi,s=1,c='pink',label='HTF entering BHE')
plt.scatter(Date,T_ho,s=1, c='green',label='HTF leaving PHE')
#plt.plot(T_initial,c='blue',label='Temperatura que devia entrar no borehole')
plt.legend(markerscale=6,fontsize=15)
plt.xlabel('Date', fontsize=24)
plt.ylabel('Temperature [Celsius]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.figure()
plt.scatter(q,h_ho,s=1,c='orange',label='heat exchanger')
plt.scatter(q,h_hi,s=1,c='pink',label='valve')
plt.scatter(q,h_wi,s=1, c='green',label='entering borehole')
#plt.plot(q, a3*q**2 + b3*q**1 + c3,c='red')
#plt.plot(q, a4*q**2 + b4*q**1 + c4,c='magenta')
#plt.plot(q, a5*q**2 + b5*q**1 + c5,c='blue')

#plt.scatter(q,h_should,s=1, c='blue',label='h_blue')
plt.legend()

plt.show()

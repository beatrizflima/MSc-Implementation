import pandas as pd
from FunctionsUFH import convHT, Prandtl
from matplotlib import pyplot as plt
from cmath import log, sqrt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from sympy import Symbol, solve, log
from sympy import symbols, Eq, solve
from CoolProp.CoolProp import PropsSI
from interpolateAMRhot import T_hi

#------------------------------Functions------------------------------#
def reynolds_plate(m_dot,rho,L,N,mu,D_h):
    u_m = m_dot/(rho*L**2/2)
    G = m_dot / ((N-1)*(gap*L/3))
    Re_D = G*D_h/mu
    return Re_D

def T_log(T_hi,T_co,T_ho,T_ci):
    T_log = np.divide(T_hi-T_co-(T_ho-T_ci),np.log((T_hi-T_co)/(T_ho-T_ci)))
    return T_log

def equations(p,T_hi,T_ho,T_ln_in,T_ci):
    x = p
    return np.divide(T_hi-x-(T_ho-T_ci),np.log((T_hi-x)/(T_ho-T_ci)))-T_ln_in


#------------------------Inputs from DataFrame------------------------#
df_HE_temp = pd.read_csv('ext_floor.csv', delimiter=',')
df_HE_temp['Time'] = pd.to_datetime(df_HE_temp['Time'])

Q_needed = df_HE_temp['Q_total']
T_ci = df_HE_temp['T_wo']
T_co = df_HE_temp['T_wi']
T_wi_initial = df_HE_temp['T_wi']
Time = df_HE_temp['Time']
#plt.scatter(Time,T_wi_initial,s=1,c='blue',label='should underfloor')

#---------------------------Constant Inputs---------------------------#
T_hi = T_hi                        #C
T_ho = 35                          #C
lambda_w = 0.61                    #W/mK
m_UFH = 0.08                       #kg/s
m_HP = [0]*len(Q_needed)
rho = 994                          #kg/m^3
mu_UFH = 0.0008                    #Pa.s
mu_HP = 0.000705                   #Pa.s
D_i = 0.012                        #m
c_w = 4190                         #J/(kgK)
N = 70
#D_h = 2*(L/N)
gap=0.0015 #m

Re_HP = [0]*len(Q_needed)
T_ln = [0]*len(Q_needed)
UA = [0]*len(Q_needed)
m_HE = [0.08]*len(Q_needed)
m_valve = [0]*len(Q_needed)
h_ci = [0]*len(Q_needed)
h_co = [0]*len(Q_needed)
h_wi = [0]*len(Q_needed)
T_wi = [0]*len(Q_needed)
mu = 0.001 
h_w = 7.54*lambda_w

max_mass_MCHP = 0.9665530980904599
density_HTF = PropsSI('D','T',35+273.15,'P',101325,'INCOMP::MEG-20%')/1000
mu_HTF = PropsSI('V','T',35+273.15,'P',101325,'INCOMP::MEG-20%')
lambda_HTF = PropsSI('L','T',35+273.15,'P',101325,'INCOMP::MEG-20%')
h_htf = 7.54 * lambda_HTF

for i in range(len(Q_needed)):
    if Q_needed[i] == 3000 :
        T_ln_in = T_log(T_hi[i],T_co[i],T_ho,T_ci[i])
        UA = Q_needed[i]/T_ln_in

L = (-366*N*gap*h_w*h_htf+366*gap*h_w*h_htf+sqrt((366*N*gap*h_w*h_htf-366*gap*h_w*h_htf)**2-4*(122*N*h_w*h_htf-122*h_w*h_htf)*(-600*gap*UA*h_w-600*gap*UA*h_htf)))/(2*(122*N*h_w*h_htf-122*h_w*h_htf))
D_h = 4 * (gap*L/3)/(2*(gap+L/3))

print('max Re MCHP',reynolds_plate(max_mass_MCHP,density_HTF,L,N,mu_HTF,D_h))
print('L',L)
print('D_h',D_h)

for i in range(len(Q_needed)):
    if Q_needed[i] > 0 :
        T_ln[i] = np.divide(Q_needed[i],UA)
        T_co[i] = fsolve(equations,T_hi[i]-0.0001,args=(T_hi[i],T_ho,T_ln[i],T_ci[i]))
        m_HE[i] = Q_needed[i]/(c_w*(T_co[i]-T_ci[i]))
        m_valve[i] = m_UFH - m_HE[i]
        h_co[i] = PropsSI('H','T',T_co[i]+273.15,'P',101325,'Water')
        h_ci[i] = PropsSI('H','T',T_ci[i]+273.15,'P',101325,'Water')
        h_wi[i] = (m_valve[i]*h_ci[i] + m_HE[i]*h_co[i])/m_UFH
        T_wi[i] = PropsSI('T','H',h_wi[i],'P',101325,'Water')-273.15
        Re_HP[i] = reynolds_plate(m_HE[i],rho,L,N,mu,D_h)
        
print('max Re UFH',max(Re_HP))

print('L',L)
print('D_h',D_h)


plt.scatter(Time,T_co,s=1,c='magenta',label='HTF leaving the PHE')
plt.scatter(Time,T_ci,s=1,c='green',label='HTF leaving the UFH')
plt.scatter(Time,T_wi,s=1,c='red',label='HTF entering UFH')
plt.xlabel('Date', fontsize=24)
plt.ylabel('Temperature [Celsius]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(markerscale=6,fontsize=15)
plt.figure()
plt.scatter(T_ci,T_co,s=1,c='pink',label='Floor Temperature')
plt.figure()
plt.scatter(Q_needed,m_HE,s=1,c='red',label='Mass flow rate through PHE1')
plt.scatter(Q_needed,m_valve,s=1,c='cyan',label='Mass flow rate directly to T1')
plt.xlabel('Heating Power [W]', fontsize=24)
plt.ylabel('Mass Flow Rate [kg/s]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(markerscale=6,fontsize=15)
plt.figure()
plt.scatter(Q_needed,h_co,s=1,c='red',label='heat exchanger')
plt.scatter(Q_needed,h_ci,s=1,c='pink',label='valve')
plt.scatter(Q_needed,h_wi,s=1,c='green',label='entering underfloor')
plt.legend(markerscale=4,fontsize=15)
plt.figure()
plt.scatter(Q_needed,T_co,s=1,c='red',label='heat exchanger')
plt.scatter(Q_needed,T_ci,s=1,c='pink',label='valve')
plt.scatter(Q_needed,T_hi,s=1,c='orange',label='leaving MCHP')
plt.scatter(Q_needed,[T_ho]*len(T_co),s=1,c='yellow',label='entering MCHP')
plt.scatter(Q_needed,T_wi,s=1,c='green',label='entering underfloor')
plt.legend(markerscale=4,fontsize=15)
plt.figure()
plt.scatter(Q_needed,T_ln,s=1,c='pink',label='Floor Temperature')
plt.scatter(Q_needed,T_wi,s=1,c='pink',label='Floor Temperature')


plt.legend(markerscale=4,fontsize=15)
#plt.show()

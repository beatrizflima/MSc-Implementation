import seaborn as sns
from cmath import sqrt
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from FunctionsUFH import iterationTemperatures, water_PL, Nusselt, N_T_U, eff, Reynolds, pressuredrop,darcy,heat_room_ext, Prandtl,Nusselt,convHT
from TimevsTemp import Temperature09_10, Date09_10, Time09_10
from ThermalResistances import PL_resis,pipe_resis,water_resis

## --------------------------------------------Inputs-----------------------------------------------##
N = 8                               #number of circuits
Q_dot_max = 3000/N                  #W
T_ext = Temperature09_10            #ºC
T_room = 20                         #ºC
T_floor_max = 29                    #ºC
rho = 994                          #kg/m^3
m_dot = 0.01                       #kg/s
c_w = 4190                          #J/(kgK)
lambda_w = 0.61                     #W/mK
mu=0.001                            #Pa.s
D_i = 0.012                         #m
delta_pipe = 0.002                  #m
D_o = D_i + 2*delta_pipe
lambda_pipe = 0.41                  #W/mK
w = 0.15                            #m
A_rad = 100/N                       #m^2
L_pipe = A_rad/(w)                  #m
rough = 7*10**(-6)                  #relative roughness of the pipes
lambda_screed = 1.4                 #W/mK
delta_up = 0.045                     #m
t_step = 300                        #s
k_floor = 0.14                      #W/mK
delta_floor=0.05                    #m
rho_floor = 650                     #kg/m3
vol_floor = delta_floor*A_rad       #m3
cp_floor = 1200                     #J/kgK
rho_screed = 1200                   #kg/m3
vol_screed = delta_up*A_rad         #m3
cp_screed = 840                     #J/kgK
k_screed = 0.41                     #W/mK
D = 4*A_rad/(sqrt(A_rad)*4)         #hydraulic diameter
rho_air = 1.204                     #kg/m^3
cp_air = 1005                       #J/kgK
vol_air = A_rad*2.6                 #m3

##------------------------------------Calculated Variables-----------------------------------------##
Q_dot, T_room, UA_house = heat_room_ext(21,T_ext,-8,Q_dot_max, Time09_10)    #heating power through each circuit
T_wi = [35]*len(Q_dot)               
T_wo = [25]*len(Q_dot)
T_pl = [28]*len(Q_dot)
Q_total = N*Q_dot

for i in range(len(Q_dot)):
    if Q_dot[i]>0:   
        if T_ext[i] <= -8:
            T_wi[i] = 35.7
        else:
            T_wi[i] = (-13.7/26)*T_ext[i] + 31.484615384615

        T_wo[i] = -Q_dot[i]/(m_dot*c_w)+T_wi[i]
        
        R_w = water_resis(w,D_i,delta_pipe,L_pipe,m_dot)
        R_pipe = pipe_resis(w,D_i,delta_pipe,lambda_pipe).real
        R_x = PL_resis(w,D_i,lambda_screed).real

        v_dot = m_dot/rho       #m^3/s
        Re = Reynolds(m_dot,mu,D_i) 
        f = darcy(Re)
        Pr = Prandtl(lambda_w,mu,c_w)
        Nu = Nusselt(Pr,Re,f)
        h_conv = convHT(Nu,D_i,lambda_w)
        p_drop = pressuredrop(f,L_pipe,D_i,rho,v_dot)               #Pascal
        UA_pipe = 1/(R_w+R_pipe+R_x)*A_rad
        NTU = N_T_U(UA_pipe,m_dot,c_w)
        epsilon = eff(NTU)
        T_pl[i] = water_PL(epsilon,T_wi[i],T_wo[i])

    elif Q_dot[i]==0:
        T_wi[i] = 19.5
        T_wo[i] = -Q_dot[i]/(0.000001*c_w)+T_wi[i]

        R_w = water_resis(w,D_i,delta_pipe,L_pipe,0.000001)
        R_pipe = pipe_resis(w,D_i,delta_pipe,lambda_pipe).real
        R_x = PL_resis(w,D_i,lambda_screed).real

        UA_pipe = 1/(R_w+R_pipe+R_x)*A_rad
        NTU = N_T_U(UA_pipe,0.000001,c_w)
        epsilon = eff(NTU)
        T_pl[i] = water_PL(NTU,T_wi[i],T_wo[i])

pump = p_drop*N*m_dot*N/(rho*0.5)
T_inter,T_floor,T_pl,T_room = iterationTemperatures(UA_house,T_ext,rho_air,cp_air,vol_air,T_pl,Q_dot,A_rad,k_floor,delta_floor,t_step,delta_up,T_room,rho_floor,vol_floor,cp_floor,rho_screed,vol_screed,cp_screed,k_screed,D)

#--------------------------------------Get COP and Q_BHE---------------------------------------------#

COP = pd.read_csv(r'C:\Users\beali\Downloads\pie_script (1)\list_res.csv', delimiter=',')
COP = COP['COP']
W_HP = [0]*COP
Q_BHE = [0]*COP

for i in range(len(COP)):
    if COP[i] != 0:
        W_HP[i] = Q_total[i]/COP[i]
        Q_BHE[i] = Q_total[i] - W_HP[i]
print('max COP',max(COP))
plt.scatter(Q_total,COP,s=1)
plt.xlabel('Heating Power (W)', fontsize=24)
plt.ylabel('COP', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.figure()


#----------------------------------------Values to CSV-----------------------------------------------#

list_ext_floor = list(zip(Time09_10, Date09_10,T_ext,T_floor,T_room,T_inter,T_wo,T_wi,T_pl,Q_total,Q_BHE, COP))
list_ext_floor = pd.DataFrame(list_ext_floor, columns=['timestamp','Time','T_ext','T_floor','T_room','T_inter','T_wo','T_wi','T_pl','Q_total','Q_BHE','COP'])
list_ext_floor['T_floor'] = np.real(list_ext_floor['T_floor'])


Q_total.to_csv('qh_needed.csv',index=False, header=False)
list_ext_floor.to_csv('ext_floor.csv')
#---------------------------------Hourly Values for the BHE-----------------------------------------#
hourly_entries = pd.DataFrame([el for el in list_ext_floor['Time'] if el.minute == 0],columns=['Time'])

hourly_entries = pd.merge(hourly_entries, list_ext_floor, on = ['Time'], how='left')
hourly_entries.to_csv('hour_power.csv')

#print(hourly_entries)

##-----------------------------------------Plots-------------------------------------------------##
#plt.scatter(Time,T_wo,s=1,c='blue',label='Temperature water out')

plot_winter = list_ext_floor.loc[list_ext_floor['timestamp'] >= 76319]
Date09_10 = plot_winter['Time']
T_ext = plot_winter['T_ext']
T_floor = plot_winter['T_floor']
T_room = plot_winter['T_room']
T_inter = plot_winter['T_inter']
T_wo = plot_winter['T_wo']
T_wi = plot_winter['T_wi']
T_pl = plot_winter['T_pl']
Q_total = plot_winter['Q_total']
Q_BHE = plot_winter['Q_BHE']

print('The average power is ',np.average(Q_total))
print('The average room temperature is ',np.average(T_room))
print('The min room temperature is ',np.min(T_room))
print('The max room temperature is ',np.max(T_room))
        
min_temp = plot_winter.loc[plot_winter['T_ext'] >= -8]
T_room_min = min_temp['T_room']
print('The min room temperature above -8 is ',np.min(T_room_min))

#plt.scatter(Date09_10,T_wi,s=1,c='red',label='Inlet Water Temperature')
#plt.scatter(Date09_10,T_wo,s=1,c='green',label='Outlet Water Temperature')
#plt.scatter(Date09_10,T_floor,s=1,c='pink',label='Floor Temperature')
plt.scatter(Date09_10,T_room,s=1,c='blue',label='Room Temperature')
plt.legend(markerscale=6,fontsize=15)
#plt.scatter(Date09_10,[21]*len(Date09_10),s=1,c='red',label='21')
#plt.scatter(Time,T_inter,s=1,c='red',label='Temperature between screed and floor')
#plt.scatter(Time,T_pl,s=1,c='green',label='Temperature between screed and floor')
#plt.legend()
plt.xlabel('Date', fontsize=24)
plt.ylabel('Room Temperature [Celsius]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.figure()
plt.scatter(Date09_10,Q_total,s=1,c='orange')
#plt.scatter(Date09_10,Q_BHE,s=1,c='pink',label='Heat load on the borehole side')
#plt.scatter(Time,Q_total,s=1,c='blue',label='Total heat load')
plt.xlabel('Date', fontsize=24)
plt.ylabel('Heating Power [W]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend()
plt.figure()

sns.histplot(data=T_room, stat="percent",binwidth = 0.25,binrange=(18.0,21),kde=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Temperature bin [Celsius]',fontsize=24)
plt.ylabel('Percentage [%]',fontsize=24)
plt.figure()

plt.scatter(T_ext,Q_total/N,s=1)
plt.xlabel('Outside Temperature [Celsius]')
plt.ylabel('Total Heating Load [W]')
plt.figure()

plt.scatter(T_ext,T_floor,s=1,c='green',label='Temperature Floor Surface')
plt.scatter(T_ext,T_inter,s=1,c='pink',label='Temperature Interface between Floor cover and screed')
plt.scatter(T_ext,T_room,s=1,c='red',label='Temperature room')

#plt.scatter(T_ext,T_pipe_ext,s=1,c='orange',label='Temperature of the exterior of the pipe')
#plt.scatter(T_ext,T_pipe_int,s=1,c='red',label='Temperature of the interior of the pipe')
#plt.scatter(T_ext,T_wm,s=1,c='black',label='Mean Temperature of the water')
plt.scatter(T_ext,T_wi,s=1,c='grey',label='Temperature of the water in')
plt.scatter(T_ext,T_wo,s=1,c='blue',label='Temperature of the water out')
plt.scatter(T_ext,T_pl,s=1,c='cyan',label='Temperature of the PL')
plt.legend(markerscale=6,fontsize=15)
plt.figure()

plt.show()

##-----------------------------------------Prints------------------------------------------------##

print("UA of the house: ", UA_house)
print("The length of each circuit is: ", L_pipe)
print("The Reynolds number is: ", Re)
print("The Prandtl number is: ", Pr)
print("The Nusselt number is: ", Nu)
print("The convective heat transfer is: ", h_conv)
print("The pressure drop is", p_drop)
print("The pump power is", pump)
print("The NTU is: ", NTU)
print("The effectiveness is: ", epsilon)






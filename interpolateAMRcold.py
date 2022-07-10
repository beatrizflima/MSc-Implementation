import pandas as pd
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI
import numpy as np

res = pd.read_csv(r'C:\Users\beali\Downloads\pie_script (1)\list_res.csv', delimiter=',')
q_needed = pd.read_csv("ext_floor.csv",delimiter=',')
qc_needed = q_needed['Q_BHE']
qh_needed = q_needed['Q_total']
date = pd.to_datetime(q_needed['Time'])
density = PropsSI('D','T',8+273.15,'P',101325,'INCOMP::MEG-20%')/1000
cp = PropsSI('C','T',8+273.15,'P',101325,'INCOMP::MEG-20%')
lambda_HTF = PropsSI('L','T',8+273.15,'P',101325,'INCOMP::MEG-20%')

flow = res['mass']
v_flow = [0]*len(qc_needed)
T_ci = [0]*len(qc_needed)


for i in range(len(qc_needed)):
    if qh_needed[i]>820:
        v_flow[i] = flow[i]*(2*68/7)/60
    elif 500<qh_needed[i]<=820:
        v_flow[i] = flow[i]*(2*34/7)/60
    elif 0<qh_needed[i]<=500:
        v_flow[i] = flow[i]*(2*17/7)/60

for i in range(len(qc_needed)):
    if qc_needed[i]==0:
        T_ci[i] = 8
    else:
        T_ci[i] = 8 - qc_needed[i] / (v_flow[i]*density*cp)
        

list_Tci_BHE = list(zip(date,qh_needed,qc_needed,T_ci))
list_Tci_BHE = pd.DataFrame(list_Tci_BHE, columns=['Time','Q_total','Q_BHE','T_ci'])

hourly_entries = pd.DataFrame([el for el in list_Tci_BHE['Time'] if el.minute == 0],columns=['Time'])
hourly_entries = pd.merge(hourly_entries, list_Tci_BHE, on = ['Time'], how='left')

T_ci = hourly_entries['T_ci']

plt.scatter(hourly_entries['Time'],T_ci, s=1)
plt.xlabel('Date', fontsize=24)
plt.ylabel('Temperature [Celsius]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()
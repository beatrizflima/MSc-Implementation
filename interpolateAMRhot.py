from cProfile import label
import pandas as pd
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI
import numpy as np

res = pd.read_csv(r'C:\Users\beali\Downloads\pie_script (1)\list_res.csv', delimiter=',')
qh_needed = np.loadtxt(open("qh_needed.csv","rb"))
q_needed = pd.read_csv("ext_floor.csv",delimiter=',')
date = pd.to_datetime(q_needed['Time'])


density = PropsSI('D','T',35+273.15,'P',101325,'INCOMP::MEG-20%')/1000
lambda_HTF = PropsSI('L','T',35+273.15,'P',101325,'INCOMP::MEG-20%')
cp = PropsSI('C','T',35+273.15,'P',101325,'INCOMP::MEG-20%')


flow = res['mass']
cop = res['COP']
frequency = res['f']
v_flow = [0]*len(qh_needed)
m_flow = [0]*len(qh_needed)

T_hi = [0]*len(qh_needed)
n_reg = [0]*len(qh_needed)


for i in range(len(qh_needed)):
    if qh_needed[i]>820:
        v_flow[i] = flow[i]*(2*68/7)/60
        n_reg[i] = 68
    elif 500<qh_needed[i]<=820:
        v_flow[i] = flow[i]*(2*34/7)/60
        n_reg[i] = 34

    elif 0<qh_needed[i]<=500:
        v_flow[i] = flow[i]*(2*17/7)/60
        n_reg[i] = 17


for i in range(len(qh_needed)):
    m_flow[i] = v_flow[i]*density
    T_hi[i] = qh_needed[i] / (m_flow[i]*cp) + 35

list_q_cop_m_f = list(zip(qh_needed,m_flow,frequency,n_reg))
list_q_cop_m_f = pd.DataFrame(list_q_cop_m_f, columns=['Q_total','m_flow','Freq','Reg'])
list_q_cop_m_f.to_csv('list_q_cop_m_f.csv')


plt.scatter(date,T_hi,s=1,c='red')
plt.xlabel('Date', fontsize=24)
plt.ylabel('Temperature [Celsius]', fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.figure()
plt.scatter(qh_needed,frequency,s=1,label='Frequency [Hz]')
plt.scatter(qh_needed,m_flow,s=1,label='Mass flow [kg/s]')
plt.scatter(qh_needed,cop,s=1,label='COP')
plt.legend(markerscale=6,fontsize=20)
plt.xlabel('Heating power [W]',fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


plt.show()
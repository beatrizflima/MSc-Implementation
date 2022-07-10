from email import header
import pandas as pd
from pandas import read_csv
import math
import numpy as np
from scipy.interpolate import interpolate,Rbf
import csv
import matplotlib.pyplot as plt

# configuration
qh_needed = np.loadtxt(open("pie_script\qh_needed.csv","rb"))
qh_headers = np.loadtxt(open("pie_script\qh.csv","rb"), delimiter=",")
cop_headers = np.loadtxt(open("res.csv","rb"), delimiter=",")
qc = np.loadtxt(open("pie_script\qc.csv","rb"), delimiter=",")

#configure q_needed:
qh_reg = [0]*len(qh_needed)

for i in range(len(qh_needed)):
    if qh_needed[i]>820:
        qh_reg[i] = qh_needed[i]/68
    elif 500<qh_needed[i]<=820:
        qh_reg[i] = qh_needed[i]/34
    elif 0<qh_needed[i]<=500:
        qh_reg[i] = qh_needed[i]/17

qh_reg = pd.DataFrame(qh_reg)
qh_reg.to_csv('pie_script\qh_needed_reg.csv',header=True,index=False)
#flow rate and frequency
m_dot = qh_headers[0,1:]
f = qh_headers[1:,0]


#remove flow rate row and frequency column
qh = qh_headers[:,1:len(qh_headers)]
qh = qh[1:len(qh_headers),:]

qc = qc[:,1:len(qh_headers)]
qc = qc[1:len(qh_headers),:]

cop = cop_headers[:,1:len(cop_headers)]
cop = cop[1:len(cop_headers),:]

#calculate number of regenerators needed
#reg = np.max(qh_needed)/np.max(qh)
#reg = math.ceil(reg)                    #rounds up the number to next unit number

#divide total qh_needed by number of reg
#qh_needed = qh_needed/reg



#X, Y = np.meshgrid(m_dot,f)

interp = interpolate.interp2d(m_dot,f,qh,kind='linear')
xnew = np.linspace(m_dot.min(), m_dot.max(), 100)
ynew = np.linspace(f.min(), f.max(), 100)
znew = interp(xnew, ynew)
print(interp(1.1,1)) #I'm obtaining the q for a value of m_dot and f, but I want the other way


plt.contourf(m_dot, f, qc,levels = np.linspace(0,np.max(qc),100))
plt.colorbar()
plt.xlabel('Flow Rate (L/min)',fontsize=20)
plt.ylabel('Frequency (Hz)',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.figure()
plt.contourf(m_dot, f, qh, levels = np.linspace(0,np.max(qh),100))
plt.colorbar()
plt.xlabel('Flow Rate (L/min)',fontsize=20)
plt.ylabel('Frequency (Hz)',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.figure()
plt.contourf(m_dot, f, cop,levels = np.linspace(0,np.max(cop),100))
plt.colorbar()
plt.xlabel('Flow Rate (L/min)',fontsize=20)
plt.ylabel('Frequency (Hz)',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.figure()
plt.scatter(qh*68,cop)
plt.show()
print('max qc',np.max(qc),'max qh',np.max(qh),'max cop',np.max(cop))
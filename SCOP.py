from matplotlib.pyplot import plot
import pandas as pd
import math
from scipy.interpolate import interp1d
import numpy as np


values_power_temp = pd.read_csv(r'hour_power.csv', delimiter=',')
values_m_f = pd.read_csv(r'list_q_cop_m_f.csv', delimiter=',')

plot_winter = values_power_temp.loc[values_power_temp['timestamp'] >= 76295]
plot_m_f_winter = values_m_f.loc[values_m_f.index >= 35424]

plot_winter = pd.merge(plot_winter, plot_m_f_winter, on = ['Q_total'], how='left')
del plot_winter['Unnamed: 0_y']
del plot_winter['Unnamed: 0_x']

plot_winter=plot_winter.drop_duplicates(subset='Time')
print(plot_winter)
#print('hours 68 reg ',plot_winter.groupby(['Reg']).count())

heat_demand =  plot_winter['Q_total']
t_ext = plot_winter['T_ext']

plot_winter = plot_winter.assign(bins = plot_winter['T_ext'].apply(np.floor))

f = interp1d(t_ext,heat_demand)

plot_winter_mean = plot_winter.groupby('bins')[['Q_total','COP','m_flow','Reg']].mean()
hours = plot_winter.groupby('bins')[['Time']].count()
print(hours)
#print(plot_winter_mean)
plot_winter_mean = pd.merge(hours,plot_winter_mean, on = ['bins'], how='left')
plot_winter_mean['Q_total'] = plot_winter_mean['Q_total']/1000
plot_winter_mean = plot_winter_mean.assign(El_cons_h = plot_winter_mean['Q_total']/plot_winter_mean['COP']).fillna(0)
plot_winter_mean = plot_winter_mean.assign(Total_heat = plot_winter_mean['Q_total']*plot_winter_mean['Time']).fillna(0)
plot_winter_mean = plot_winter_mean.assign(Total_El_cons = plot_winter_mean['El_cons_h']*plot_winter_mean['Time']).fillna(0)
plot_winter_mean['m_flow'] = plot_winter_mean['m_flow']
plot_winter_mean['Reg'] = plot_winter_mean['Reg']


total_heating_demand = plot_winter_mean['Total_heat'].sum()
total_electricity_year = plot_winter_mean['Total_El_cons'].sum()
total_hours = plot_winter_mean['Time'].sum()
SCOP = total_heating_demand/total_electricity_year

rounded = plot_winter_mean.round(3)

print(rounded)
print('the total heating demand during the year is ', total_heating_demand, 'Wh')
print('the total electricity consumed during the year is ', total_electricity_year, 'Wh')
print('the SCOP of the winter season is ', SCOP)
print('the total number of hours is ', total_hours)

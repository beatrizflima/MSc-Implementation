from operator import index
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
import numpy as np

data_temp = pd.read_csv('TimevsDatavsSoil.csv', delimiter=';')
data_temp['Time'] = pd.to_datetime(data_temp['Time'], dayfirst=True)

plt.scatter(data_temp['Time'],data_temp['Temperature'],s=1)
plt.xlabel('Date',fontsize=24)
plt.ylabel('Outdoor Temperature [Celsius]',fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.figure()
time_stamp = data_temp.set_index('Time')['Temperature'].resample('5T').interpolate('linear')
columns_name = ['Temperature']
time_stamp = pd.DataFrame(time_stamp, columns=columns_name)
time_stamp = time_stamp.index

Time = np.linspace(0,87646,num=87647, endpoint=True)
Temperature = data_temp['Temperature']
f = interp1d(Time,Temperature)

Time_5min = np.linspace(0,87646,num=1051753, endpoint=True)
Temperature_5min = f(Time_5min)
plt.plot(Time, Temperature, 'o', Time_5min, f(Time_5min), '-')
plt.legend(['data', 'linear'],loc='best')


new_data = list(zip(Time_5min,time_stamp,f(Time_5min)))
new_data = pd.DataFrame(new_data, columns=['Time','Date','Temperature'])

winter09_10 = new_data.loc[new_data['Time'] >= 73343]
winter09_10 = winter09_10.loc[winter09_10['Time'] < 82103]
winter09_10=winter09_10.reset_index()
Temperature09_10 = winter09_10['Temperature']
Date09_10 = winter09_10['Date']
Time09_10 = winter09_10['Time']
#print(winter09_10)
winter09_10.to_csv('vizualize_temp_list.csv')
plt.figure()
plt.scatter(Date09_10,Temperature09_10,s=1)
plt.xlabel('Date',fontsize=24)
plt.ylabel('Outdoor Temperature [Celsius]',fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.show()

##-----------------------------------------Air Temperature------------------------------------------##
#Time = data_temp.index
#Temperature = data_temp['Temperature']
T_min_ext=min(Temperature)
T_max_ext=max(Temperature)

##-----------------------------------------Soil Temperature------------------------------------------##

data_temp_soil = pd.read_csv('Soil2001_2010.csv', delimiter=';')
data_temp_soil['Time'] = pd.to_datetime(data_temp_soil['Time'])

time_stamp_soil = data_temp_soil.set_index('Time')['Soil'].resample('5T').interpolate('linear')
columns_name_soil = ['Soil']
time_stamp_soil = pd.DataFrame(time_stamp_soil, columns=columns_name_soil)
time_stamp_soil = time_stamp_soil.index

Time_soil = np.linspace(0,14608,num=14609, endpoint=True)
Soil_Temperature = data_temp_soil['Soil']
g = interp1d(Time_soil,Soil_Temperature)

Time_soil_5min = np.linspace(0,14608,num=1051753, endpoint=True)
Soil_Temperature_5min = g(Time_soil_5min)
#plt.figure()
#plt.plot(Time_soil, Soil_Temperature, 'o', Time_soil_5min, g(Time_soil_5min), '-')
#plt.legend(['data', 'linear'],loc='best')
plt.show()


new_data_soil = list(zip(Time_soil_5min,time_stamp_soil,g(Time_soil_5min)))
new_data_soil = pd.DataFrame(new_data_soil, columns=['Time','Date','Soil'])
new_data_soil.to_csv('vizualize_temp_soil.csv')

winter09_10_soil = new_data_soil.loc[new_data_soil['Time'] >= 12224.112269812656]
winter09_10_soil = winter09_10_soil.loc[winter09_10_soil['Time'] < 13684.145585651371]


winter09_10_soil=winter09_10_soil.reset_index()
Temperature09_10_soil = winter09_10_soil['Soil']
Date09_10_soil = winter09_10_soil['Date']
winter09_10.to_csv('vizualize_winter.csv')

'''
##------------------------------------------------Plots----------------------------------------------##

#plt.plot(Time, dif_T)
#plt.figure()

n, bins, patches = plt.hist(Temperature, color = 'blue', edgecolor = 'black',bins=range(-14,36))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title('Histogram of Temperatures per hour',fontsize=20)
plt.xlabel('Temperature bin [Celsius]',fontsize=24)
plt.ylabel('Number of Hours',fontsize=24)
plt.figure()

sns.distplot(Temperature, hist=True, kde=True, bins=range(-15,35),color = 'darkblue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 2})
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title('Histogram of Temperatures per hour',fontsize=20)
plt.xlabel('Temperature bin [Celsius]',fontsize=24)
plt.ylabel('Frequency',fontsize=24)
plt.figure()

sns.histplot(Temperature,cumulative = True, bins=range(-15,35))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title('Cumulative Histogram of Temperatures per hour',fontsize=20)
plt.xlabel('Temperature bin [Celsius]',fontsize=24)
plt.ylabel('Number of Hours',fontsize=24)


plt.figure()
print(T_min_ext,T_max_ext)
plt.scatter(Time_5min,f(Time_5min),s=1)

#plt.scatter(Time,Soil_Temperature,s=1)

#plt.show()
#print(n, bins, sum(n[0:6]))
print(len(Temperature_5min))
print(len(Time_5min))
print(len(time_stamp))

#print(len(Soil_Temperature))
'''
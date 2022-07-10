import matplotlib.pyplot as plt
from numpy import average
from TimevsTemp import Temperature09_10_soil, Date09_10_soil
from CoolProp.CoolProp import PropsSI


def Temp_depth(depth,Temperature09_10_soil):
    T_depth = Temperature09_10_soil + 20*depth
    return T_depth


depth = 0.120 #km
T_depth = Temp_depth(depth,Temperature09_10_soil)   #celsius
T_mean_soil = (T_depth + Temperature09_10_soil)/2   #celsius - average undisturbed soil temperature till specified meters depth
print(min(T_depth),min(T_mean_soil), average(Temperature09_10_soil),average(T_mean_soil),average(T_depth))
T_mean_soil.to_csv('T_depth_ground.csv')
#plt.plot(Date09_10_soil, T_depth)
'''av = [average(Temperature09_10_soil)]*len(Date09_10_soil)
plt.figure()
plt.plot(Date09_10_soil, Temperature09_10_soil, c='blue',label='Instantaneous temperature')
plt.plot(Date09_10_soil,av, c='orange',label='Yearly average temperature')
plt.xlabel('Date',fontsize=18)
plt.ylabel('Ground Temperature (Celsius)',fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

plt.show()'''

from cmath import log, pi, e, sqrt
from operator import abs
from scipy.optimize import fsolve
import numpy as np


def Reynolds(m_dot,mu,D_i):
    Re = (4*m_dot)/(mu * pi * D_i )

    return Re

def Prandtl(k,mu,c_w):
    Pr = mu*c_w/k

    return Pr

def Nusselt(Pr,Re,f):
    #Nu = 0.023*Re**(4/5)*Pr**(1/3)
    if Re>= 3000:
        Nu = ((f/8)*(Re-1000)*Pr)/(1+12.7*abs(f/8)**0.5*(Pr**(2/3)-1))
    elif Re<3000:
        Nu = 3.66
    return Nu

def convHT(Nu,D_i,k):
    h_conv = Nu*k/D_i

    return h_conv

def N_T_U(UA_pipe,m_dot,c_w):
    NTU = abs(UA_pipe/(m_dot*c_w))

    return NTU

def eff(NTU):
    epsilon = 1 - e**(-NTU)

    return epsilon

def darcy(Re):
    #f = 0.316/Re**0.25
    if Re>= 3000:
        f = abs((0.790*log(Re)-1.64)**(-2))
    elif Re<3000:
        f = 64/Re
    return f

def pressuredrop(f,L_pipe,D_i,rho,v_dot):
    p_drop = f*L_pipe/D_i*rho*(v_dot**2/(pi*D_i**2/4)**2)/2
    return p_drop

def heat_room_ext(T_room,T_ext,T_ext_min,Q_dot_max,Time09_10):
    UA_room_ext = Q_dot_max/(T_room-T_ext_min)
    Q_room_ext = UA_room_ext*(T_room-T_ext)
    T_room = [T_room]*len(T_ext)

        
    for i in range(len(T_ext)):
        if i>0:
            if T_ext[i] > 18:
                Q_room_ext[i] = 0
            
            if T_ext[i] <= T_ext_min:
                Q_room_ext[i] = Q_dot_max
                T_room[i] = Q_room_ext[i]/UA_room_ext + T_ext[i]
            
            if 73343<= Time09_10[i] < 76295 :
                Q_room_ext[i] = 0
            

    return Q_room_ext, T_room, UA_room_ext

def iterationTemperatures(UA_house,T_ext,rho_air,cp_air,vol_air,T_pl,Q_dot,A_rad,k_floor,delta_floor,t_step,delta_up,T_room,rho_floor,vol_floor,cp_floor,rho_screed,vol_screed,cp_screed,k_screed,D):
    T_inter = [25]*len(Q_dot)
    #T_inter[0] = -Q_dot[0]/(A_rad*k_screed/delta_up) + T_pl[0]
    T_floor = [23]*len(Q_dot)
    #T_floor[0]=-Q_dot[0]/(A_rad*k_floor/delta_floor)+T_inter[0]
    h=10.8

    for i in range(len(Q_dot)):
        if i > 0:
            T_inter[i] = (T_inter[i-1]*(delta_floor*delta_up*rho_screed*vol_screed*cp_screed-A_rad*(delta_floor*k_screed+k_floor*delta_up)*t_step)+2*A_rad*(delta_floor*T_pl[i]*k_screed+T_floor[i]*k_floor*delta_up)*t_step)/(delta_floor*delta_up*rho_screed*vol_screed*cp_screed+A_rad*(delta_floor*k_screed+k_floor*delta_up)*t_step)

    for i in range(len(Q_dot)):
        if i > 0:
        #Crank Nicolson
            T_floor[i]=(k_floor*(T_floor[i-1]*(rho_floor*vol_floor*cp_floor-A_rad*h*t_step)+2*A_rad*h*T_room[i]*t_step)-A_rad*t_step*(T_floor[i-1]-2*T_inter[i]))/(k_floor*(rho_floor*vol_floor*cp_floor+A_rad*h*t_step)+A_rad*t_step)
    
    for i in range(len(Q_dot)):
        if i > 0:
            T_room[i]=(2*A_rad*t_step*h*T_floor[i]-A_rad*t_step*h*T_room[i-1]-UA_house*t_step*T_room[i-1]+2*T_ext[i]*UA_house*t_step+2*rho_air*cp_air*vol_air*T_room[i-1])/(2*rho_air*cp_air*vol_air+A_rad*t_step*h+UA_house*t_step)
    
    for i in range(len(Q_dot)):
        if i > 0:
            T_inter[i] = (T_inter[i-1]*(delta_floor*delta_up*rho_screed*vol_screed*cp_screed-A_rad*(delta_floor*k_screed+k_floor*delta_up)*t_step)+2*A_rad*(delta_floor*T_pl[i]*k_screed+T_floor[i]*k_floor*delta_up)*t_step)/(delta_floor*delta_up*rho_screed*vol_screed*cp_screed+A_rad*(delta_floor*k_screed+k_floor*delta_up)*t_step)

    for i in range(len(Q_dot)):
        if i > 0:
        #Crank Nicolson
            T_floor[i]=(k_floor*(T_floor[i-1]*(rho_floor*vol_floor*cp_floor-A_rad*h*t_step)+2*A_rad*h*T_room[i]*t_step)-A_rad*t_step*(T_floor[i-1]-2*T_inter[i]))/(k_floor*(rho_floor*vol_floor*cp_floor+A_rad*h*t_step)+A_rad*t_step)
    
    for i in range(len(Q_dot)):
        if i > 0:
            T_room[i]=(2*A_rad*t_step*h*T_floor[i]-A_rad*t_step*h*T_room[i-1]-UA_house*t_step*T_room[i-1]+2*T_ext[i]*UA_house*t_step+2*rho_air*cp_air*vol_air*T_room[i-1])/(2*rho_air*cp_air*vol_air+A_rad*t_step*h+UA_house*t_step)
    

    for i in range(len(Q_dot)):
        if Q_dot[i] == 0:
            T_pl[i] = (T_pl[i-1]*(2*delta_up*rho_screed*cp_screed*vol_screed-A_rad*k_screed*t_step)+2*A_rad*k_screed*T_inter[i]*t_step)/(2*delta_up*rho_screed*cp_screed*vol_screed+A_rad*k_screed*t_step)

    for i in range(len(Q_dot)):
        if i > 0:
            T_inter[i] = (T_inter[i-1]*(delta_floor*delta_up*rho_screed*vol_screed*cp_screed-A_rad*(delta_floor*k_screed+k_floor*delta_up)*t_step)+2*A_rad*(delta_floor*T_pl[i]*k_screed+T_floor[i]*k_floor*delta_up)*t_step)/(delta_floor*delta_up*rho_screed*vol_screed*cp_screed+A_rad*(delta_floor*k_screed+k_floor*delta_up)*t_step)

    for i in range(len(Q_dot)):
        if i > 0:
        #Crank Nicolson
            T_floor[i]=(k_floor*(T_floor[i-1]*(rho_floor*vol_floor*cp_floor-A_rad*h*t_step)+2*A_rad*h*T_room[i]*t_step)-A_rad*t_step*(T_floor[i-1]-2*T_inter[i]))/(k_floor*(rho_floor*vol_floor*cp_floor+A_rad*h*t_step)+A_rad*t_step)
    
    for i in range(len(Q_dot)):
        if i > 0:
            T_room[i]=(2*A_rad*t_step*h*T_floor[i]-A_rad*t_step*h*T_room[i-1]-UA_house*t_step*T_room[i-1]+2*T_ext[i]*UA_house*t_step+2*rho_air*cp_air*vol_air*T_room[i-1])/(2*rho_air*cp_air*vol_air+A_rad*t_step*h+UA_house*t_step)

    
    return T_inter,T_floor,T_pl,T_room
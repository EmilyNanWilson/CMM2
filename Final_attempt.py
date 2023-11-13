# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:45:39 2023

@author: emily
"""

# importing modules
import numpy as np
import math
from scipy import optimize
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp



#--------------------------------------
# importing module with aerodynamics model 
# as tables of CL CM vs alpha
#    tables of CD vs CL
# and tables of CL_el and CM_el vs delta_el 
#--------------------------------------
import sys
sys.path.append('C:/Users/emily/OneDrive - University of Edinburgh/Third Year/Computational Methods/Design Project')
import aero_table
import env
import vehicle
#--------------------------------------
# Lift vs alpha

rho = env.air_density
V = 100              #[m/s]
S = vehicle.Sref
W = vehicle.acMass *env.gravity          #[N]
gamma = np.radians(2.864789)       #[rad]

# Initial guesses for the fitting
# i.e., initial values of a and b
CL_0 = 0.0410
CL_alpha = 0.1

# Functional form to fit to data
def CL_a_func(x, a, b):
    return a + b * x

# Fitting (find a and b of the above function)  
# using the python module optimize from scipy.
# params contains the fitted values of a and b
# params_covariance contains a measure of the achieved 
# accuracy 
params, params_covariance = optimize.curve_fit(CL_a_func, aero_table.alpha, aero_table.CL,
        p0=[CL_0, CL_alpha])

CL_0 = params[0]
CL_alpha = params[1]
#--------------------------------------

x_values_CL_a = np.linspace(min(aero_table.alpha), max(aero_table.alpha), 10)
y_values_CL_a = CL_a_func(x_values_CL_a,CL_0,CL_alpha)

#--------------------------------------
#Lift vs delta_elevator
CL_delta = 0.003

def CL_d_func(x, a):
    return a * x

params, params_covariance = optimize.curve_fit(CL_d_func, aero_table.delta_el, aero_table.CL_el,
        p0=[CL_delta])

CL_delta = params[0]
#--------------------------------------

x_values_CL_d = np.linspace(min(aero_table.delta_el), max(aero_table.delta_el), 10)
y_values_CL_d = CL_d_func(x_values_CL_d,CL_delta)


#--------------------------------------
# CD vs CL
CD_0 = 0.026
CD_k = 0.045

def CD_CL_func(x, a, b):
    return a + b * x**2.0

params, params_covariance = optimize.curve_fit(CD_CL_func, aero_table.CL, aero_table.CD,
        p0=[CD_0, CD_k])

CD_0 = params[0]
CD_k = params[1]
#--------------------------------------
x_values_CD_CL = np.linspace(min(aero_table.CL), max(aero_table.CL), 10)
y_values_CD_CL = CD_CL_func(x_values_CD_CL,CD_0,CD_k)


#--------------------------------------
# Moment vs alpha
CM_0 = 0.
CM_alpha = -0.01

# TO BE COMPLETED HERE
#--------------------------------------

def CM_alpha_func(x, a, b):
    return a + b * x

params, params_covariance = optimize.curve_fit(CM_alpha_func, aero_table.alpha, aero_table.CM,
        p0=[CM_0, CM_alpha])

CM_0 = params[0]
CM_alpha = params[1]

x_values_CM_a = np.linspace(max(aero_table.alpha), min(aero_table.alpha), 10)

y_values_CM_a = CM_alpha_func(x_values_CM_a, CM_0, CM_alpha)

#--------------------------------------
#Moment vs delta_elevator
CM_delta = -0.004

# TO BE COMPLETED HERE
#--------------------------------------
def CM_delta_func(x, a):
    return a*x

params, params_covariance = optimize.curve_fit(CM_delta_func, aero_table.delta_el, aero_table.CM_el, p0 = [CM_delta])

CM_delta =params[0]

#--------------------------------------
x_values_CM_d = np.linspace(min(aero_table.delta_el), max(aero_table.delta_el), 10)
y_values_CM_d = CM_delta_func(x_values_CM_d,CM_delta)

#--------------------------------------
# Write results on screen (check)
print("CL0, ", CL_0, " CL_alpha, ",CL_alpha, " CL_delta, ",CL_delta)
print("CD_0 ", CD_0, " CD_k ",CD_k)
print("CM_0 ", CM_0, " CM_alpha ",CM_alpha, " CM_delta, ", CM_delta)
#--------------------------------------

#Plotting Results of A1
def plot(x_data, y_data, x_label, y_label, title, x_fitted, y_fitted):
    plt.figure()
    plt.scatter(x_data, y_data, color='purple')
    plt.plot(x_fitted, y_fitted, color='red')
    plt.xlabel(x_label, fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.title(title, fontsize=16)
    plt.show()

plot(aero_table.alpha, aero_table.CL, 'alpha', 'CL_wing',
          'CL_wing = C_L_0 + C_L_alpha*alpha', x_values_CL_a, y_values_CL_a)

plot(aero_table.delta_el, aero_table.CL_el, 'delta_el', 'CL_el',
          'CL_el = CL_d * delta_el', x_values_CL_d, y_values_CL_d)

plot(aero_table.CL, aero_table.CD, 'CL', 'CD',
          'CD = CD_0 + K*(CL^2)', x_values_CD_CL, y_values_CD_CL)

plot(aero_table.alpha, aero_table.CM, 'alpha', 'CM_wing',
          'CM_wing = CM_0 + CM_alpha*alpha', x_values_CM_a, y_values_CM_a)

plot(aero_table.delta_el, aero_table.CM_el, 'delta_el', 'CM_delta',
          'CM_el = CM_delta*delta_el', x_values_CM_d, y_values_CM_d)

#----------------------
# PART A2
#----------------------

def pitch_function(alpha,V,CM):
    return 0.5*rho*(V**2)*S*CM

def CM_function(alpha):
    return CM_0 + (CM_alpha*alpha)+(CM_delta*delta)

def delta_function(alpha): 
    return (-((CM_0 + (CM_alpha*alpha))/CM_delta))

def Cl_function(alpha): 
    delta = delta_function(alpha)
    return CL_0 + (CL_alpha*alpha) + (CL_delta*delta)

def lift_function(alpha):
    Cl = Cl_function(alpha)
    return 0.5*rho*(V**2)*S*Cl

def CD_function(Cl): 
    return CD_0 + (CD_k * (Cl**2))

def drag_function(alpha):
    CD = CD_function(alpha)
    return 0.5*rho*(V**2)*S*CD



def A2(alpha, V, gamma):
    
    alpha = alpha[0]
    alpha_deg = np.degrees(alpha)
    Cl = Cl_function(alpha_deg)
    CD = CD_function(Cl)
    L = 0.5*rho*(V**2)*S*Cl
    D = 0.5*rho*(V**2)*S*CD
    return - (L * np.cos(alpha)) - (D * np.sin(alpha)) + (W * np.cos(alpha + gamma))

alpha_solved = fsolve(A2, 0,args =  (V,gamma))

alpha_deg = np.degrees(alpha_solved)
print ("Alpha is:...... ", alpha_solved, "   [rad]")


def thrust_function(alpha_solved,gamma):
    alpha_deg = np.degrees(alpha_solved)
    Cl = Cl_function(alpha_deg)
    CD = CD_function(Cl)
    L = 0.5*rho*(V**2)*S*Cl
    D = 0.5*rho*(V**2)*S*CD
    return (vehicle.acMass * env.gravity) * math.sin(alpha_solved + gamma) + D * math.cos(alpha_solved) - L * math.sin(alpha_solved)



delta = np.radians(delta_function(alpha_deg))

Thrust = thrust_function(alpha_solved, gamma)
#- lift_function(np.degrees(alpha_solved)) * np.sin(np.degrees(alpha_solved)) + drag_function(np.degrees(alpha_solved)) * np.cos(np.degrees(alpha_solved)) + W * np.sin(np.degrees(alpha_solved) + np.degrees(gamma))
q = 0
theta = (delta) + (alpha_solved)
u_b = V * np.cos(alpha_solved)
w_b = V * np.sin(alpha_solved)


print('delta_e is', delta)
print('Thrust is', Thrust)
print('q is', q)
print('theta is', theta)
print('u_b is', u_b)
print('w_b is', w_b)



#--------------------------------
# PART A3
#--------------------------------
new_delta = 0

#Part A3: define function for dub/dt

def model(t,y):
    if t>=100:
        new_delta = (delta*1.1)*180/math.pi
        pass
    else:
        new_delta = (delta)*180/math.pi
        pass
    '''
    y1 = y[0] # u_b
    y2 = y[1] # w_b
    y3 = y[2] # q
    y4 = y[3] # theta
    y5 = y[4] # Xe
    y6 = y[5] # Ze
    alpha_deg = np.degrees(np.arctan(w_b/u_b))
    V_new = np.sqrt(u_b**2 + w_b**2)
    CM =  CM_0 + (CM_alpha*alpha_deg)+(CM_delta*new_delta)
    Cl = CL_0 + (CL_alpha*alpha_deg) + (CL_delta*new_delta)
    CD = CD_0 + (CD_k * (Cl**2))
    L_new = 0.5*rho*(V_new**2)*S*Cl
    D_new = 0.5*rho*(V_new**2)*S*CD
    M_new = 0.5*rho*(V_new**2)*S*CM*vehicle.cbar
    f1 = ((L_new/vehicle.acMass)*np.sin(alpha_deg))-((D_new/vehicle.acMass)*np.cos(alpha_deg))-(y3*y2-(W/vehicle.acMass)*np.sin(theta)) +(Thrust/vehicle.acMass)
    f2 = ((-L_new/vehicle.acMass)*np.cos(alpha_deg))-((D_new/vehicle.acMass)*np.sin(alpha_deg))+(y3*y1)+((W/vehicle.acMass)*np.cos(theta))
    f3 = (M_new/vehicle.inertia_yy)
    f4 = y3
    f5 = (y1*np.cos(y4))+(y2*np.sin(y4))
    f6 = (-y1*np.sin(y4)) + (y2*np.cos(y4))
    return f1,f2,f3,f4,f5,f6
'''
    u_b, w_b, theta,q, x_e,z_e = y

    alpha = np.arctan(w_b/u_b) 
    V = np.sqrt(u_b**2 + w_b**2)

    W = (vehicle.acMass * env.gravity)
    S = vehicle.Sref
    rho = env.air_density
    cbar = vehicle.cbar

    alpha_deg = np.degrees(alpha)
    delta_deg = np.degrees(delta)

    CL = CL_0 + (CL_alpha * alpha_deg) + (CL_delta * delta_deg)
    CD = CD_0 + CD_k*(CL**2)
    CM = CM_0 + (CM_alpha * alpha_deg) + (CM_delta * delta_deg)

    L = 0.5 * rho * (V**2) * S * CL
    D = 0.5 * rho * (V**2) * S * CD
    M = 0.5 * rho * (V**2) * S * CM * cbar

    dq_dt = (M/vehicle.inertia_yy)
    dtheta_dt = q

    du_dt = (L*math.sin(alpha) - D*math.cos(alpha) - vehicle.acMass*q*w_b - W*math.sin(theta) + Thrust)/vehicle.acMass
    dw_dt = (-L*math.cos(alpha)-D*math.sin(alpha)+vehicle.acMass*q*u_b + W*math.cos(theta))/vehicle.acMass

    dx_dt = u_b*math.cos(theta) + w_b*math.sin(theta)
    dz_dt = - u_b*math.sin(theta) + w_b*math.cos(theta)

    return du_dt,dw_dt,dtheta_dt,dq_dt,dx_dt,dz_dt

t0 = 0
y1_0 = 1.646 
y2_0 = 99.986 
y3_0 = 0
y4_0 = 0.0164
y5_0 = 0
y6_0 = 0

t_final = 300

t_eval = np.linspace(t0,t_final,5000)
# ------------------------------------------------------
# total solution interval
# step size
# not needed here. The solver solve_ivp 
# will take care of finding the appropriate step 
# ------------------------------------------------------

# ------------------------------------------------------
# Apply solve_ivp method

soln1 = solve_ivp(model, [t0,t_final] , [u_b,w_b,theta,0,0,0])
y = soln1.y



  
fig1, ax1 = plt.subplots() 
ax1.plot(soln1.t,soln1.y[0, :]) 
ax1.set_title("u_b") 

fig2, ax2 = plt.subplots() 
ax2.plot(soln1.t,soln1.y[1, :]) 
ax2.set_title("w_b") 
  
fig3, ax3 = plt.subplots() 
ax3.plot(soln1.t,soln1.y[2, :]) 
ax3.set_title("q") 

fig4, ax4 = plt.subplots() 
ax4.plot(soln1.t,soln1.y[3, :]) 
ax4.set_title("theta") 


fig5, ax5 = plt.subplots() 
ax5.plot(soln1.t,soln1.y[4, :]) 
ax5.set_title("xe") 

fig6, ax6 = plt.subplots() 
ax6.plot(soln1.t,soln1.y[5, :]) 
ax6.set_title("ze") 
# Combine all the operations and display 
plt.show() 

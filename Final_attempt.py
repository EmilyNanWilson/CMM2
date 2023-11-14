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
sys.path.append('C:/Users/emily/OneDrive - University of Edinburgh/Third Year/Computational Methods/Design Project 2')
import aero_table
import env
import vehicle
#--------------------------------------
# Lift vs alpha

rho = env.air_density
V = 100              #[m/s]
S = vehicle.Sref
W = vehicle.acMass *env.gravity          #[N]
gamma = 0.05       #[rad]

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
params, params_covariance = optimize.curve_fit(CL_a_func, np.radians(aero_table.alpha), aero_table.CL,
        p0=[CL_0, CL_alpha])

CL_0 = params[0]
CL_alpha = params[1]
#--------------------------------------

#--------------------------------------
#Lift vs delta_elevator
CL_delta = 0.003

def CL_d_func(x, a):
    return a * x

params, params_covariance = optimize.curve_fit(CL_d_func, np.radians(aero_table.delta_el), aero_table.CL_el,
        p0=[CL_delta])

CL_delta = params[0]
#--------------------------------------

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


#--------------------------------------
# Moment vs alpha
CM_0 = 0.
CM_alpha = -0.01

# TO BE COMPLETED HERE
#--------------------------------------

def CM_alpha_func(x, a, b):
    return a + b * x

params, params_covariance = optimize.curve_fit(CM_alpha_func, np.radians(aero_table.alpha), aero_table.CM,
        p0=[CM_0, CM_alpha])

CM_0 = params[0]
CM_alpha = params[1]


#--------------------------------------
#Moment vs delta_elevator
CM_delta = -0.004

# TO BE COMPLETED HERE
#--------------------------------------
def CM_delta_func(x, a):
    return a*x

params, params_covariance = optimize.curve_fit(CM_delta_func,np.radians( aero_table.delta_el), aero_table.CM_el, p0 = [CM_delta])

CM_delta =params[0]


#--------------------------------------
# Write results on screen (check)
print("CL0, ", CL_0, " CL_alpha, ",CL_alpha, " CL_delta, ",CL_delta)
print("CD_0 ", CD_0, " CD_k ",CD_k)
print("CM_0 ", CM_0, " CM_alpha ",CM_alpha, " CM_delta, ", CM_delta)
#--------------------------------------


#----------------------
# PART A2
#----------------------



def delta_function(alpha):
    return -((CM_0 + (CM_alpha * alpha))  / CM_delta)

def CL_function(alpha):
    delta = delta_function(alpha) 
    return CL_0 + (CL_alpha * alpha) + (CL_delta * delta)

def CD_function(alpha):
    CL = CL_function(alpha)
    return CD_0 + (CD_k * (CL**2))

def CM_function(alpha):
    delta = delta_function(alpha)
    return CM_0 + (CM_alpha * alpha) + (CM_delta * delta)



def A2(alpha, V, gamma):
    CL = CL_function(alpha)
    CD = CD_function(alpha)
    L = 0.5 * env.air_density * V**2 * vehicle.Sref * CL
    D = 0.5 * env.air_density * V**2 * vehicle.Sref * CD
    W = vehicle.acMass * env.gravity
    return -L * np.cos(alpha) - D * np.sin(alpha) + W * np.cos(alpha + gamma)


alpha_radians = fsolve(A2, x0 = 0, args=(V, gamma))

print ("Alpha is:...... ", alpha_radians, "   [rad]")


delta_radians = delta_function(alpha_radians)

CL = CL_0 + CL_alpha*alpha_radians + CL_delta * delta_radians
CD = CD_0 + CD_k*(CL**2)


L = 0.5 * rho * V**2 * S * CL
D = 0.5 * rho * V**2 * S * CD


theta_radians = gamma + alpha_radians

Thrust = -L*np.sin(alpha_radians)+D*np.cos(alpha_radians)+W*np.sin(theta_radians)
q = 0
theta = (gamma) + (alpha_radians)

u_b = V * np.cos(alpha_radians)
w_b = V * np.sin(alpha_radians)


print('delta_e is   ', delta_radians, "   [rad]")
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

def A3(t,y):
    y1 = y[0] # u_b
    y2 = y[1] # w_b
    y3 = y[2] # q
    y4 = y[3] # theta
    y5 = y[4] # Xe
    y6 = y[5] # Ze
    if t>=100:
        delta_A3 = (delta_radians*1.1)
        
    else:
        delta_A3 = (delta_radians)
    
    
    thrust_A3 = Thrust
    M = vehicle.acMass
    W = env.gravity * M
    
    alpha_A3 = np.arctan(y2/y1)
    V_A3 = np.sqrt((y1**2) + (y2**2))
    
    Cl_A3 = CL_0 + (CL_alpha*alpha_A3) + (CL_delta*delta_A3)
    CM_A3 =  CM_0 + (CM_alpha*alpha_A3)+(CM_delta*delta_A3)
    CD_A3 = CD_0 + (CD_k * (Cl_A3**2))
    
    L_A3 = 0.5*rho*(V_A3**2)*S*Cl_A3
    D_A3 = 0.5*rho*(V_A3**2)*S*CD_A3
    M_A3 = 0.5*rho*(V_A3**2)*S*CM_A3*vehicle.cbar
    

    
    y1_dt = ((L_A3 / M_A3) * np.sin(alpha_A3)) - ((D_A3 / M_A3) * np.cos(alpha_A3)) - (y3 * y2) - ((W / M_A3) * np.sin(y4)) + (thrust_A3 / M_A3)
    y2_dt =  ((-L_A3 / M_A3) * np.cos(alpha_A3)) - ((D_A3 / M_A3) * np.sin(alpha_A3)) + (y3 * y1) + (W / M_A3) * np.cos(y4)
    y3_dt = (M_A3/vehicle.inertia_yy)
    y4_dt = y3
    y5_dt = (y1 * np.cos(y4)) + (y2 * np.sin(y4))
    y6_dt = -(y1 * np.sin(y4)) + (y2 * np.cos(y4))
    return y1_dt,y2_dt,y3_dt,y4_dt,y5_dt,y6_dt


t0 = 0
y1_0 = u_b 
y2_0 = w_b
y3_0 = 0
y4_0 = theta_radians
y5_0 = 0
y6_0 = -2000

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

soln1 = solve_ivp(A3, [t0,t_final] , [y1_0,y2_0,y3_0,y4_0,y5_0,y6_0])
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

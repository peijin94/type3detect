'''
    File name: radioTools.py
    Author: Peijin Zhang 张沛锦
    Date : 2022-4-26
    
    handy tools to convert amoung frequency density heliocentric distance
'''


import numpy as np
import math

from scipy.optimize import fsolve as fsolver 
const_c = 2.998e10

R_S = 6.96e10         # the radius of the sun 
c   = 2.998e10        # speed of light
c_r = c/R_S           # [t]

def omega_pe_r(ne_r,r):
    # plasma frequency density relationship
    return 8.93e3* (ne_r(r))**(0.5) * 2 * np.pi

def saito77(r):
    return 1.36e6 * r**(-2.14) + 1.68e8 * r**(-6.13)

def leblanc98(r):
    return 3.3e5* r**(-2.)+ 4.1e6 * r**(-4.)+8.0e7* r**(-6.)

def parkerfit(r):
    h0=144.0/6.96e5
    h1=20.0/960.
    nc=3e11*np.exp(-(r-1.0e0)/h1)
    return  4.8e9/r**14. + 3e8/r**6.+1.39e6/r**2.3+nc

def dndr_leblanc98(r):
    return -2.*3.3e5* r**(-3.) -4.*4.1e6 * r**(-5.) -6.*8.0e7* r**(-7.)

def newkirk(r):
    return 4.2e4*10. **(4.32/r)

def f_Ne(N_e): 
    # in Hz
    return 8.93e3 * (N_e)**(0.5)

def Ne_f(f):
    # in cm-3
    return (f/8.93e3)**2.0


def freq_to_R(f_pe, ne_r = parkerfit):
    """
    Starting height for wave frequency f
    """
    func  = lambda R : f_pe - (omega_pe_r(ne_r,R)) /2/np.pi
    R_solution = fsolver(func, 1.5) # solve the R
    return R_solution # [R_s]


def R_to_freq(R,ne_r = parkerfit):
    """
    Wave frequency from R
    """
    return omega_pe_r(ne_r,R) /2/np.pi 



def freq_drift_f_t(t,v,t0,dm = parkerfit ):
    # t0 is the time at 100MHz
    # t [day]
    # f [MHz]
    # v [c]
    r0 =  freq_to_R(300e6,dm)
    dt = (t-t0)
    dR = dt * v *  c_r#/np.sqrt(1-v**2)
    R_new = r0+dR
    return  R_to_freq(R_new, 
                ne_r = dm)/1e6



def freq_drift_t_f(f,v,t0,dm = parkerfit):
    # t0 is the time at 100MHz
    # t [day]
    # f [MHz]
    # v [c]
    
    R_new =  freq_to_R(f*1e6, 
                ne_r = dm)
    dt = (R_new- freq_to_R(300e6,dm))/v/ c_r
    return t0+dt
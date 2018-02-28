#!/usr/bin/python
######################################################################
# Name:         controller
# Author:       A. Marocchino
# Date:			25-02-2016
# Purpose:      controller for PWFA utilities: input list and plots
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy, pylab
import scipy as sci
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.constants import codata
#from mpl_toolkits.mplot3d import Axes3D
# - #
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes','Plasma_PyCalculator','formulary'))
sys.path.append(os.path.join(home_path,'Codes','Plasma_PyCalculator','PWFA'))
from plasma_basic_parameters import *
from PWFA_calculator import *
# --- #

# --- #
%matplotlib inline


n0=1.e16*1e6
alpha=1.83
print(3e10/electron_plasma_frequency(n0))


#--- ---#
n0=1.2e18*1e6
alpha=1.83
print(Q_from_alpha_bigaussian(alpha,n0,3.9e-6,3.9e-6,6.9e-6)*1e12)
print(Q_from_alpha_bigaussian(4*alpha,1e16*1e6,20.5e-6,20.5e-6,75e-6)*1e12)
print(reduced_charge(581e-12,n0,3.9e-6,3.9e-6,6.9e-6))
print(reduced_charge(5.82e-9,1e16*1e6,20.5e-6,20.5e-6,75e-6))

print(Q_from_alpha_bigaussian(alpha,n0,3.9e-6,3.9e-6,6.9e-6)*3e8/(np.sqrt(2.*np.pi)*6.9e-6)/1e3)
lambda_p = electron_plasma_wavelength(n0)
kp= electron_plasma_wavenumber(n0)
print(electron_plasma_wavelength(n0)*1e6)
print(electron_plasma_wavenumber(n0)*1e-6)
print(betatron_wavelength(1e18*1e6,200.)*1e3)
print(Q_from_alpha_bigaussian(130.0,n0,0.017/kp,0.027/kp,0.13/kp)*1e12)
print(10.*10.*np.sqrt(n0/1e16/1e6))
print(9.1e-31*3e8*electron_plasma_frequency(n0)/1.6e-19/1e9)
print(me*c*wp/e)
os.quit()

#--- Simple Example ---#
np=3e13*1e6 #density in m^-3
kp    = electron_plasma_wavenumber(np)
lambda_p = electron_plasma_wavelength(np)
print(kp*30e-6)

Q=Q_from_alpha_bigaussian(15.0,n0,1.5e-6,1.5e-6,7.5e-6)
I=Q/(np.sqrt(2.*np.pi)*7.5e-6/3e8)
print(I/1e3,Q*1e9)

print(Q_from_alpha_bigaussian(16.05,1e17*1e6,8.4e-6,8.4e-6,16.8e-6))
print(matching_condition_transverse(1e17*1e6,1e-6,2000)*1e6)
print(Q_from_alpha_trapezoidalZ_gaussianR(350.00,800.00,1e17*1e6,1.e-6,1.e-6,30e-6)*1e9)


print('-----')
Vol1 = (2*pi)**1.5 * 2.6**2*40.
Vol2 = (2*pi)**1.0 * 1.5**2*20.
Sup1 = pi * 2.6 * 40.
Sup2 = sqrt(2.*pi) * 1.5 * 20.
print(Vol1)
print(Vol2)
print(Sup1)
print(Sup2)
print(electron_plasma_wavelength(1.3e16*1e6)*1e6/2)
print('-----')

print('-----')
np=1e16*1e6 #density in m^-3
kp    = electron_plasma_wavenumber(np)
omegap= electron_plasma_frequency(np)
lambda_p = electron_plasma_wavelength(np)
print(matching_condition_transverse(np,4e-6,2000)*1e6)
print('-----')


#--- Shock injection parameters ---#
np=1e18*1e6 #density in m^-3
kp    = electron_plasma_wavenumber(np)
omegap= electron_plasma_frequency(np)
lambda_p = electron_plasma_wavelength(np)
print(kp*1e-6)
print(480*0.2e-6*kp)
print(520*0.2e-6*0.75*kp)
print(0.08e-6/3e8*omegap)
print(Q_from_alpha_bigaussian(3.3,np,7.5e-6,4e-6,4e-6)*1e9)
print(Q_from_alpha_bigaussian(3.3,np,7.5e-6,4e-6,4e-6) / (sqrt(2.*3.14)*7.5e-6/c) )
print(electron_plasma_wavelength(1e18*1e6)*1e6)
print(electron_plasma_wavelength(1.75*1e18*1e6)*1e6)
print(electron_plasma_wavelength(1.25*1e18*1e6)*1e6)
print(electron_plasma_wavelength(1.07*1e18*1e6)*1e6)
print(0.025*33)
print(100e9*3e-6/1e6)
print('normalization for total bunch density')
Qnorm=1.6e-19*np/kp**3
Q=0.125e-9/Qnorm
rho=2.*Q/(kp*3e-6)
print(Qnorm)
print(Q)
print(rho)
#--- --- --- --- --- --- --- --- --- ---#



#--- Dimensionless parameter Woosley ---#
Te = 10*1.16e4
ni = 0.59e6
Z=1
A=1.0079
L = 3e18
U = 1e9
kv = kinematic_viscosity(Te,ni,Z,A)
chi=kinematic_coefficient_thermal_diffusivity(Te,ni,Z,A)
Re=U*L/kv
Pe=U*L/chi
print(Re)
print(Pe)
L=100e-6
U=2e4
Te=5e5
ni=3e17*1e6
Z=18
A=39.9480
print(U*L/kinematic_coefficient_thermal_diffusivity(Te,ni,Z,A))
L=100e-6
U=3e4
Te=9e4
ni=3e18*1e6
Z=18
A=39.9480
print(U*L/kinematic_coefficient_thermal_diffusivity(Te,ni,Z,A))

msys.exit('stopped script')

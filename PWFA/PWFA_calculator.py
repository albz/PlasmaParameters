#!/usr/bin/python
######################################################################
# Name:         PWFA utilities
# Author:       A. Marocchino
# Date:			24-11-2015
# Purpose:      matching condition for PWFA
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy, pylab
from scipy import *
import numpy as np
from scipy.constants import codata
# - #
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/PlasmaParameters'))
from plasma_basic_parameters import *
# --- #

### --- constants --- ###
physical_const = codata.physical_constants #SI
eps0   = physical_const['electric constant'][0]
try:
	mu0    = physical_const['mag. constant'][0]
except:
	mu0    = physical_const['magnetic constant'][0]
c      = physical_const['speed of light in vacuum'][0]
me     = physical_const['electron mass'][0]
qe     = physical_const['elementary charge'][0]
proton_mass             = physical_const['proton mass'][0]
boltzmann_constant_JK   = physical_const['Boltzmann constant'][0]
boltzmann_constant_eVK  = physical_const['Boltzmann constant in eV/K'][0]


### --- ###
def density_ramp(z,Lramp):
	if z<=Lramp:
		rho = (Lramp-z)/Lramp
	else:
		rho = 0.
	return rho

def fext(z,Lramp,n0):
	f_ext = n0 * qe**2 / (me*eps0*c**2) * density_ramp(z,Lramp)
	return f_ext

def Kext(z,Lramp,n0,gamma):
	K_ext = fext(z,Lramp,n0) / (2.*gamma)
	return K_ext

def f_rk(z,u1,u2,Lramp,n0,gamma,eps):
	return -Kext(z,Lramp,n0,gamma)*u1+eps**2/gamma**2/u1**3

def g_rk(z,u1,u2,Lramp,n0,gamma,eps):
	return u2


### --- Input Condition to have matching at the end of the plasma ramp --- ###
def matching_with_ramp(sx,Lramp,Lvacuum,n0,gamma,emittance):
	h=.1e-6 #integration step
	z=0.
	u1=sx
	u2=0. #D(sz)
	while(z<Lramp+Lvacuum):
		k1=f_rk(z     ,u1        ,u2        ,Lramp,n0,gamma,emittance)
		g1=g_rk(z     ,u1        ,u2        ,Lramp,n0,gamma,emittance)

		k2=f_rk(z+h/2.,u1        ,u2+h/2.*k1,Lramp,n0,gamma,emittance)
		g2=g_rk(z+h/2.,u1+h/2.*g1,u2        ,Lramp,n0,gamma,emittance)

		k3=f_rk(z+h/2.,u1        ,u2+h/2.*k2,Lramp,n0,gamma,emittance)
		g3=g_rk(z+h/2.,u1+h/2.*g2,u2        ,Lramp,n0,gamma,emittance)

		k4=f_rk(z+h   ,u1        ,u2+h*k3   ,Lramp,n0,gamma,emittance)
		g4=g_rk(z+h   ,u1+h*g3   ,u2        ,Lramp,n0,gamma,emittance)

		z=z+h
		u2=u2+h*(k1+2.*k2+2.*k3+k4)/6.
		u1=u1+h*(g1+2.*g2+2.*g3+g4)/6.

	sigma=u1
	alphaTwiss=u1*u2/emittance*gamma
	betaTwiss=u1**2/emittance/1e-6

	return sigma,alphaTwiss,betaTwiss


### ---  matching conditions --- ###
def matching_condition(n0,emittance,gamma):
	k0 = electron_plasma_wavenumber(n0)
	sigma_matching = np.sqrt(np.sqrt(2./gamma)) * np.sqrt(emittance/k0)
	return sigma_matching
def generalised_matching_condition(alpha,n0,emittance,gamma):
	k0    = electron_plasma_wavenumber(n0)
	if( alpha > 1.):
		sigma_matching = np.sqrt(np.sqrt(2./gamma)) * np.sqrt(emittance/k0)
	else:
		k0 *= alpha
		sigma_matching = np.sqrt(np.sqrt(2./gamma)) * np.sqrt(emittance/k0)
	return sigma_matching
def matching_condition_longitudinal(n0):
	k0       = electron_plasma_wavenumber(n0)
	lambda_p = electron_plasma_wavelength(n0)
	sz = lambda_p/6.
	return sz
def matching_condition_transverse(charge,sz,n0,emittance,gamma):
	k0    = electron_plasma_wavenumber(n0)
	stored_charge=0.
	for alpha in arange(0.01,1000.,0.001):
		sx = generalised_matching_condition(alpha,n0,emittance,gamma)
		calculated_charge = Q_from_alpha(alpha,n0,sx,sx,sz)
		if( calculated_charge > charge ):
			return alpha,sx


### --- betatron oscillation wavelength --- ###
def betatron_wavelength(n0,gamma):
	lp        = electron_plasma_wavelength(n0)
	lbetatron = np.sqrt(2.*gamma)*lp
	return lbetatron

### --- calculate total charge from alpha --- #
def Q_from_alpha(alpha,n0,sx,sy,sz):
	volume	= (2.*np.pi)**(3./2.) * (sx*sy*sz)
	charge	= (alpha*n0) * volume * qe #total charge in [C]
	return charge

### --- determine alpha from total charge --- #
def alpha_from_Q(charge,n0,sx,sy,sz):
	volume	= (2.*np.pi)**(3./2.) * (sx*sy*sz)
	alpha =  charge/(n0*volume*qe)
	return alpha

### --- determining Q~ factor or reduced charge factor ---#
def reduced_charge(charge,n0,sx,sy,sz):
	volume	= (2.*np.pi)**(3./2.) * (sx*sy*sz)
	alpha   = alpha_from_Q(charge,n0,sx,sy,sz)
	kp      = electron_plasma_wavenumber(n0)
	reduced_charge = alpha * volume * kp**3
	return reduced_charge

### --- determining Q~ factor or reduced charge factor ---#
def bubble_radius(alpha,sx):
	dimension = 2.*np.sqrt(alpha)*sx
	return dimension

	### --- calculate final particle velocity : constant electric field --- ###
def particle_final_velocity(q,m,Ez,distance):
	Pz=0.0
	run_distance=0.0
	Dt=1e-16
	while run_distance<distance:
		Pz -= Dt * q*Ez
		gamma=np.sqrt(1.+(Pz/m/c)**2)
		run_distance += Dt * (Pz/m/gamma)
	return Pz/m/gamma,gamma

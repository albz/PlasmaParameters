#!/usr/bin/python
######################################################################
# Name:         PlasmaParameters
# Author:       A. Marocchino
# Date:         2015-11-17
# Purpose:      series of python function to calculate basic plasma paramteres
# Source:       python
#####################################################################

### loading shell commands
import os, sys, shutil, time
import numpy as np
import scipy as scy
import pylab as pyl
from scipy.constants import codata
# - #
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/Plasma_PyCalculator/formulary'))
from csg_basic_constants import *
### --- ###


#--- physical constants ---#
physical_const = codata.physical_constants #SI
# most common
electron_mass	 	= physical_const['electron mass'][0]
electron_charge 	= physical_const['elementary charge'][0]
proton_mass 		= physical_const['proton mass'][0]
boltzmann_constant_JK 	= physical_const['Boltzmann constant'][0]
boltzmann_constant_eVK 	= physical_const['Boltzmann constant in eV/K'][0]
c                   = physical_const['speed of light in vacuum'][0]
eps0                = physical_const['electric constant'][0]
try:
	mu0    = physical_const['mag. constant'][0]
except:
	mu0    = physical_const['magnetic constant'][0]
#--- *** ---#


### --- --- --- --- --- --- ###
# All inputs are in SI Units  #
### --- --- --- --- --- --- ###

def coulomb_logarithm(Te,ne,Z):
	TeeV  = Te*boltzmann_constant_eVK
	ne_cc = ne/1e6 #from m^-3 to cm^-3

	if TeeV >= 10.*Z**2:
		CL = 24.0 - 0.5*np.log(ne_cc) + np.log(TeeV)
	else:
		CL = 23.0 - 0.5*np.log(ne_cc) - np.log(Z) + 3./2.*np.log(TeeV)
	CL = max(CL,1.)
	return CL



def electron_plasma_frequency(ne):
	return np.sqrt(ne * electron_charge**2 / electron_mass / eps0)



def electron_plasma_wavenumber(ne):
	return electron_plasma_frequency(ne)/c



def electron_plasma_wavelength(ne):
	return 2.*np.pi/electron_plasma_wavenumber(ne)


def electron_cycloton_frequency(B):
	B_cgs=B*1e4
	return electron_charge_csg*B_cgs /c_cgs /electron_mass_cgs

def ion_cycloton_frequency(B,Z,mu):
	B_cgs=B*1e4
	return Z*electron_charge_csg*B_cgs /c_cgs /(mu*proton_mass_cgs)



def electron_collision_rate(ne,Te,Z):
	ne_cc	= ne/1e6
	TeeV 	= Te*boltzmann_constant_eVK
	CL		= coulomb_logarithm(Te,ne,Z)
	return 2.91e-6 * ne_cc * CL * TeeV**1.5

def electron_thermal_velocity(Te):
	return np.sqrt(boltzmann_constant_ergK*Te/electron_mass_cgs)

def ion_thermal_velocity(Ti,mu):
	return np.sqrt(boltzmann_constant_ergK*Ti/mu/proton_mass_cgs)




def electron_collision_times(Te,ne,Z):
		TeeV = Te*boltzmann_constant_eVK
		ne_cc = ne/1e6 #from m^-3 to cm^-3
		return 3.5e5/coulomb_logarithm(Te,ne,Z) * TeeV**1.5 /Z/ne_cc

def ion_collision_times(Ti,Te,Z,ne,mu):
		TieV = Ti*boltzmann_constant_eVK
		ne_cc = ne/1e6 #from m^-3 to cm^-3
		return 2.09e7/coulomb_logarithm(Te,ne,Z) * TieV**1.5 * np.sqrt(mu) /Z/ne_cc



def K_THermalCondution_Electron_Parallel(Te,ne,Z):
		Teerg = Te*boltzmann_constant_ergK
		ne_cc = ne/1e6 #from m^-3 to cm^-3
		return 3.16 * ne_cc * Teerg * electron_collision_times(Te,ne,Z) / electron_mass_cgs

def K_THermalConduction_Electron_Perpendicular(Te,ne,Z,Bfield):
		Teerg = Te*boltzmann_constant_ergK
		ne_cc = ne/1e6 #from m^-3 to cm^-3
		return 4.66 * ne_cc * Teerg / electron_mass_cgs / electron_collision_times(Te,ne,Z) / electron_cycloton_frequency(Bfield) / electron_cycloton_frequency(Bfield)

def K_THermalConduction_Electron_Cross(Te,ne,Z,Bfield):
		Teerg = Te*boltzmann_constant_ergK
		ne_cc = ne/1e6 #from m^-3 to cm^-3
		return 2.5 * ne_cc * Teerg / electron_mass_cgs / electron_cycloton_frequency(Bfield)


# print ( '%8.5e' % (electron_cycloton_frequency(3e5/1e4)*electron_collision_times(400*1.16e4,1e21*1e6,1.)))
# print ( '%8.5e' % (ion_cycloton_frequency(3e5/1e4,1.,2.5)*ion_collision_times(300*1.16e4,400*1.16e4,1.,1e21*1e6,2.5)))

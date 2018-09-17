# -*- coding: utf-8 -*-
"""
Updated on:13/09/18
@author: Steven Sheppard

Measuring Wavelength of Laser with Ruler!
ALL units are in m 
ALL angles are in Radians

Calculate the angles alpha and beta
Determine the wavelength
Perform error analysis.
"""
from pylab import *

def angle(S,L):
    return arctan(S/L)
def wl(n,beta,alpha):
    return (d/n) * (cos(alpha)-cos(beta))
def percdif(a,b):
    return (a-b) / b
def deviation(ave, num):
    return abs(num-ave)
def frac_unc_beta(dSi,Si):
    return 

unc_Si = 0.0001#0.0129
unc_S0 = 0.0001#0.013
unc_L  = 0.0002#
unc_d  = 0.00005

'''
Data From First Experiment:
'''
L1      = 3.45
d       = 0.001
S01     = 0.1005
alpha1  = angle(S01,L1)

#State Uncertainties:
phi1 = S01/L1
frac_phi1 = unc_S0/S01 + 2*unc_L/L1
dphi1 = phi1 * frac_phi1 

#Create data arrays:
S_i1   = zeros(4)
B_i1   = zeros(4)
wl_i1  = zeros(4)
omega1 = zeros(4)
frac_omega1 = zeros(4)
domega1     = zeros(4)
total_unc1  = zeros(4)
stdsum1     = zeros(4)

#Fill array with data:
S_i1[0] = 0.153
S_i1[1] = 0.192
S_i1[2] = 0.222
S_i1[3] = 0.250

for n in range(4):
    B_i1[n] = angle(S_i1[n],L1) 
    wl_i1[n] = wl(n+1,B_i1[n], alpha1)
    frac_omega1 = unc_Si / S_i1[n] + 2*unc_L/L1
    omega1[n] = S_i1[n]/L1
    domega1 = omega1 * frac_omega1 
    total_unc1[n] =( (unc_d * (1/(n+1)) * (cos(alpha1) - cos(B_i1[n])))**2 + \
                 (dphi1*(-d/(n+1))*sin(arctan(phi1))*1/(1+phi1**2))**2 + \
                 (domega1[n]*(-d/(n+1))*sin(arctan(omega1[n]))*1/(1+omega1[n]**2))**2)**.5
              
#Calculate the <wl> and sigma_wl:
mean_wl1 = sum(wl_i1) / len(wl_i1)
for n in range(4):
    stdsum1[n]= (wl_i1[n]-mean_wl1)**2
std1        = sqrt((1/3) * sum(stdsum1) )
sigma_mean1 = std1 / sqrt(4)
'''
Data From Second Experiment:
'''
#Measurments:
L2     = 3.44
d      = 0.001
S02    = .153
alpha2 = angle(S02,L2)

#State Uncertainties:
phi2 = S01/L1
frac_phi2 = unc_S0/S02 + 2*unc_L/L2
dphi2 = phi2 * frac_phi2 

#Build Empty arrays:
S_i2   = zeros(10)
B_i2   = zeros(10)
wl_i2  = zeros(10)
omega2 = zeros(10)
frac_omega2 = zeros(10)
domega2     = zeros(10)
total_unc2  = zeros(10)
stdsum2     = zeros(10)

#fill with data
S_i2[0] = .191
S_i2[1] = .221 
S_i2[2] = .251
S_i2[3] = .274
S_i2[4] = .296
S_i2[5] = .315
S_i2[6] = .335 
S_i2[7] = .354
S_i2[8] = .373
S_i2[9] = .39

for n in range(10):
    B_i2[n] = angle(S_i2[n],L2) 
    wl_i2[n] = wl(n+1,B_i2[n], alpha2)
    frac_omega2 = unc_Si / S_i2[n] + 2*unc_L/L2
    omega2[n] = S_i2[n]/L2
    domega2 = omega2 * frac_omega2 
    total_unc2[n] =( (unc_d * (1/(n+1)) * (cos(alpha2) - cos(B_i2[n])))**2 + \
                 (dphi2*(-d/(n+1))*sin(arctan(phi2))*1/(1+phi2**2))**2 + \
                 (domega2[n]*(-d/(n+1))*sin(arctan(omega2[n]))*1/(1+omega2[n]**2))**2)**.5

mean_wl2 = sum(wl_i2) / len(wl_i2)    
for n in range(10):
    stdsum2[n]= (wl_i2[n]-mean_wl2)**2
std2 = sqrt((1/10) * sum(stdsum2) )
sigma_mean2 = std2 / sqrt(4)

'''
Plotting Routines:
'''
plott =1
set1 = linspace(1,4,4)
set2 = linspace(1,10,10)
wl   = full(10,534)

if plott == 1:
    figure(1)
    plot(set2, wl, 'r', linewidth=3)
    errorbar(set1, wl_i1*10**9, yerr=total_unc1*10**9, xerr = None,color='g', marker='o',linestyle='none',label= 'First Experiment')
    errorbar(set2, wl_i2*10**9, yerr=total_unc2*10**9, xerr = None, color='b',linestyle='none', marker='o', label= 'Second Experiment')
    ylabel('Wavelength [nm]')
    xlabel('Data Point [Si]')
    title('Uncertainty using General Formula')
    legend()

    figure(2)
    plot(set2, wl, 'r', linewidth=3)
    errorbar(set1, wl_i1*10**9, yerr=sigma_mean1*10**9, xerr = None,color='g', marker='o',linestyle='none',label= 'First Experiment')
    errorbar(set2, wl_i2*10**9, yerr=sigma_mean2*10**9, xerr = None, color='b',linestyle='none', marker='o', label= 'Second Experiment')
    ylabel('Wavelength [nm]')
    xlabel('Data Point [Si]')
    title('Uncertainty using SDOM')
    legend()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "Source" class
=============================

  Source.py -- set "relative" path to code

"""

import os
import shutil
import math  as mth
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import date
from scipy.optimize import curve_fit

import BeamLineElement   as BLE
import BeamLine          as BL
import Particle          as Prtcl
import PhysicalConstants as PhysCnst

print()
print("SourceTest: initialising:")

print("     ----> MatPlotLib options:")
#.. Use LaTeX for plot text only if it is installed:
mpl.rc('text', usetex=(shutil.which('latex') is not None))
mpl.rcParams['text.latex.preamble']="\\usepackage{bm}"
mpl.rcParams["figure.autolayout"]=True
mpl.rcParams['figure.constrained_layout.use'] = True

cm = 1./2.54  # centimeters in inches
print("     <---- MatPlotLib done.")

#.. Physical Constants
print("     ----> Physical constants:")
constants_instance = PhysCnst.PhysicalConstants()
protonMASS         = constants_instance.mp()
print("     <---- Physical constants done.")

##! Start:
print("========  Source: tests start  ========")

##! Test creation and built-in methods:
SourceTest = 1
print()
print("SourceTest:", SourceTest, \
      " check creation and built-in methods.")

#.. __init__
print("    __init__:")
try:
    Src = BLE.Source()
except:
    print('      ----> Correctly trapped no argument exception.')
rStrt = np.array([0.,0.,0.])
vStrt = np.array([[np.pi/2.,np.pi/2.],[0.,0.]])
drStrt = np.array([0.,0.,0.])
dvStrt = np.array([0.,0.,0.])
try:
    Src = BLE.Source(rStrt, vStrt, drStrt, dvStrt)
except:
    print('      ----> Correctly trapped no Source paramters exception.')
#.. Create valid instance:
Mode  = 0
Param = [ 0.8, 2.5E14, 0.6, 1.5E-6, 2.8E-14, 10., 1., 20., 20., 15., -9999. ]
Src = BLE.Source("Source0", rStrt, vStrt, drStrt, dvStrt, Mode, Param)

#.. __repr__
print("    __repr__:")
print("      ---->", repr(Src))
print("    <---- __repr__ done.")
#.. __str__
print("    __str__:")
print(str(Src))
print("    <---- __str__ done.")

print(" <---- Creation and built in method tests done!  --------  --------")

BLE.Source.cleaninstances()

##! Next: check paramterised laser-driven source distribution:
HOMEPATH = os.getenv('HOMEPATH')
filename = os.path.join(HOMEPATH, \
                        '11-Parameters/LIONBeamLine-Params-LsrDrvn.csv')
print("     ----> Load beamline parameters from:", \
      filename)
BL.BeamLine(filename)

Src = BLE.BeamLineElement.getinstances()[1]

iRefPrtcl = Prtcl.ReferenceParticle.getinstances()
print("     ----> Load beamline parameters done:")
print(BL.BeamLine.getinstances())

SourceTest += 1
print()
print("SourceTest:", SourceTest, \
      " test parameterised laser-driven source distribution", \
      "check.")
print(Src)
print(" Test generation:")
print("     ----> First particle: KE, cosThetaPhi:", \
      Src.getParticle())

#! Next: Get data for plots:
PrtclX   = np.array([])
PrtclY   = np.array([])

PrtclKE  = np.array([])
PrtclcT  = np.array([])
PrtclT   = np.array([])
PrtclPhi = np.array([])

SrcX     = np.array([])
SrcY     = np.array([])
SrcXp    = np.array([])
SrcYp    = np.array([])
SrcRp    = np.array([])
SrcZ     = np.array([])
SrcE     = np.array([])

g_E       = Src.getLaserDrivenProtonEnergyProbDensity()
E_min_MeV = Src.getParameters()[6]
E_max_MeV = Src.getParameters()[7]

iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
p0 = iRefPrtcl.getMomentumIn(0)

print("     ----> Generate many particles:")
for i in range(100000):
    TrcSpcFrmSrc = Src.getParticleFromSource()
    PhsSpcFrmSrc = Prtcl.Particle.RPLCTraceSpace2PhaseSpace(TrcSpcFrmSrc)
    
    SrcX         = np.append(SrcX,  TrcSpcFrmSrc[0])
    SrcY         = np.append(SrcY,  TrcSpcFrmSrc[2])
    SrcZ         = np.append(SrcZ,  TrcSpcFrmSrc[4])
    SrcE         = np.append(SrcE,  TrcSpcFrmSrc[5])

    p        = mth.sqrt(np.dot(PhsSpcFrmSrc[1], PhsSpcFrmSrc[1]))
    Etot     = mth.sqrt(protonMASS**2 + \
                        np.dot(PhsSpcFrmSrc[1], PhsSpcFrmSrc[1]))
    
    KE       = Etot - protonMASS
    PrtclKE  = np.append(PrtclKE , KE)

    cTheta   = PhsSpcFrmSrc[1][2] / p
    PrtclcT  = np.append(PrtclcT , cTheta)

    Theta    = mth.acos(cTheta) * 180. / mth.pi
    PrtclT   = np.append(PrtclT , Theta)

    Phi      = mth.atan2(PhsSpcFrmSrc[1][0], PhsSpcFrmSrc[1][1])
    if Phi < 0.: Phi += 2. * mth.pi
    Phi *= 180. / mth.pi
    PrtclPhi = np.append(PrtclPhi, Phi)

    PrtclX   = np.append(PrtclX , PhsSpcFrmSrc[0][0])
    PrtclY   = np.append(PrtclY , PhsSpcFrmSrc[0][1])

    xp           = TrcSpcFrmSrc[1] * p0 / p
    yp           = TrcSpcFrmSrc[3] * p0 / p
    rp           = mth.cos(Phi)*xp + mth.sin(Phi)*yp
    SrcXp        = np.append(SrcXp, xp)
    SrcYp        = np.append(SrcYp, yp)
    SrcRp        = np.append(SrcRp, rp)
    
    
print("     <---- Done.")
    
##! Next: Make plots:
SourceTest += 1
print()
print("SourceTest:", SourceTest, \
      " plot distributions and compare with expectation.")

print("     ----> Make plots:")
today = date.today().strftime("%d/%m/%Y")

#.. ----> Kinetic energy:
print("         ----> Kinetic energy:")

fig = plt.figure(figsize=(18.*cm, 6.*cm))

n, bins, patches = plt.hist(PrtclKE, \
                            bins=100, color='k', \
                            histtype='step', label='Generated Distribution')

Ee, g_E = Src.getLaserDrivenProtonEnergyProbDensity()

dE2     = (E_max_MeV - E_min_MeV) / 100. / 2.
Ee[-1] -= dE2

# Normalise:
hist_heights, bin_edges = np.histogram(PrtclKE, bins=100)
histCntnts = np.sum(hist_heights)
g_scaled   = g_E * histCntnts

# Plot required distribution:
plt.plot(Ee, g_scaled, color='r', label='Required Distribution', linewidth=1)

y_max_cutoff = g_scaled[-1]
x            = E_max_MeV
print(" x, ymax:", x, y_max_cutoff)
plt.vlines(x, ymin=0, ymax=y_max_cutoff, color='r', \
           linestyle='-', linewidth=1)  # [MeV]

plt.xlabel('Kinetic energy (MeV)', loc='right')
plt.ylabel('Entries', loc='top')
plt.yscale("log")
plt.legend(loc="best")

#plt.title('ExpSourceTst (' + today + '): Kinetic energy distribution',
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_K.pdf')
plt.close()

print("         <---- Done.")


#.. ----> Cumulative probability:
print("         ----> Cumulative probability:")

cumPROB = []
sigT    = []
for eps in Ee:
    eps1 = eps * (1.6e-19 * 1e6)
    cumPROB.append(Src.getLaserCumProb(eps1))
    sigT.append(Src.g_theta(eps))

plt.plot(Ee, cumPROB, color='k', label='Required Distribution', linewidth=1)

plt.xlabel('Kinetic energy (MeV)', loc='right')
plt.ylabel('Cumulative probability', loc='top')
plt.yscale("linear")
plt.legend(loc="best")
#plt.title('ExpSourceTst (' + today + '): cumulative probability', \
#          fontname="Times New Roman",  size=12)

plt.savefig('99-Scratch/SourceTst_cumulativePDF.pdf')
plt.close()

#.. ----> cos(theta) and theta distribution:
print("         ----> cos(theta) and theta distribution:")

fig = plt.figure(figsize=(9.*cm, 9.*cm))
n, bins, patches = plt.hist(PrtclcT, \
                            bins=50, color='k', \
                            histtype='step')

plt.xlabel('$\\cos\\theta_S$', loc='right')
plt.ylabel('Entries', loc='top')
#plt.title('ExpSourceTst (' + today + \
#          '): $\\cos\\theta_S$ distribution', \
#          fontname="Times New Roman",  size=12)

plt.savefig('99-Scratch/SourceTst_costheta.pdf')
plt.close()

fig = plt.figure(figsize=(9.*cm, 9.*cm))
n, bins, patches = plt.hist(PrtclT, \
                            bins=50, color='k', \
                            histtype='step')

plt.xlabel('$\\theta_S$ ($^\\circ$)', loc='right')
plt.ylabel('Entries', loc='top')
#plt.title('ExpSourceTst (' + today + \
#          '): $\\theta_S$ distribution', \
#          fontname="Times New Roman",  size=12)

plt.savefig('99-Scratch/SourceTst_theta.pdf')
plt.close()

plt.plot(Ee, sigT, color='r', label='Required dependence', linewidth=1)

plt.xlabel('Kinetic energy (MeV)', loc='right')
plt.ylabel('$\\sigma_{\\theta_S}$ ($^\\circ$)', loc='top')
#plt.title('ExpSourceTst (' + today + \
#          '): $\\sigma_{\\theta_S}$ versus kinetic energy',\
#          fontname="Times New Roman",  size=12)
plt.legend(loc="best")

plt.savefig('99-Scratch/SourceTst_sigTK.pdf')
plt.close()

fig = plt.figure(figsize=(9.*cm, 8.*cm))
plt.hist2d(PrtclKE, PrtclT, \
           bins=50, norm=mpl.colors.LogNorm())
plt.colorbar()

plt.xlabel('Kinetic energy (MeV)', loc='right')
plt.ylabel('$\\theta_S$ ($^\\circ$)', loc='top')
#plt.title('ExpSourceTst (' + today + \
#          '): $(\\theta_S, K)$ distribution', \
#          fontname="Times New Roman",  size=12)

plt.savefig('99-Scratch/SourceTst_thetaK.pdf')
plt.close()

print("         <---- Done.")

#.. ----> phi distribution:
print("         ----> phi distribution:")

fig = plt.figure(figsize=(9.*cm, 6.*cm))
n, bins, patches = plt.hist(PrtclPhi, \
                            bins=50, color='k', \
                            histtype='step')

plt.xlabel('$\\phi_S$ ($^\\circ$)', loc='right')
plt.ylabel('Entries', loc='top')
#plt.title('ExpSourceTst (' + today + \
#          '): $\\phi_S$ distribution', \
#          fontname="Times New Roman",  size=12)

plt.savefig('99-Scratch/SourceTst_phi.pdf')
plt.close()

print("         <---- Done.")

#.. ----> x,y distribution:
print("         ----> (x, y) distribution:")

fig = plt.figure(figsize=(9.*cm, 8.*cm))

plt.hist2d(PrtclX, PrtclY, \
           bins=50, norm=mpl.colors.LogNorm())
plt.colorbar()

plt.xlabel('x (m)')
plt.ylabel('y (m)')
#plt.title('ExpSourceTst (' + today + \
#          '): (x, y) distribution', \
#          fontname="Times New Roman",  size=12)

plt.savefig('99-Scratch/SourceTst_xy.pdf')
plt.close()

#.. ----> Trace space
print("         ----> Trace space distributions:")

fig = plt.figure(figsize=(9.*cm, 9.*cm))

plt.hist(SrcX, bins=50, color='k', histtype='step')
plt.xlabel('$x$ (m)')
plt.ylabel('Entries')
#plt.title('ExpSourceTst (' + today + \
#          '): $x$ distribution', \
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_x.pdf')
plt.close()

fig = plt.figure(figsize=(9.*cm, 9.*cm))

plt.hist(SrcY, bins=50, color='k', histtype='step')
plt.xlabel('$y$ (m)')
plt.ylabel('Entries')
#plt.title('ExpSourceTst (' + today + \
#          '): $x$ distribution', \
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_y.pdf')
plt.close()

fig = plt.figure(figsize=(9.*cm, 6.*cm))
plt.hist(SrcXp, bins=50, color='k', histtype='step')
plt.xlabel('$x^{\\prime}$ (m)')
plt.ylabel('Entries')
#plt.title('ExpSourceTst (' + today + \
#          '): $x^{\\prime}$ distribution', \
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_xp.pdf')
plt.close()

fig = plt.figure(figsize=(9.*cm, 6.*cm))
plt.hist(SrcYp, bins=50, color='k', histtype='step')
plt.xlabel('$y^{\\prime}$ (m)')
plt.ylabel('Entries')
#plt.title('ExpSourceTst (' + today + \
#          '): $y^{\\prime}$ distribution', \
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_yp.pdf')
plt.close()

fig = plt.figure(figsize=(9.*cm, 6.*cm))
plt.hist(SrcRp, bins=50, color='k', histtype='step')
plt.xlabel('$r^{\\prime}$ (m)')
plt.ylabel('Entries')
#plt.title('ExpSourceTst (' + today + \
#          '): $r^{\\prime}$ distribution', \
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_rp.pdf')
plt.close()

plt.hist(SrcZ, bins=50, color='k', histtype='step')
plt.xlabel('$z$ (m)')
plt.ylabel('Entries')
#plt.title('ExpSourceTst (' + today + \
#          '): $z$ distribution', \
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_z.pdf')
plt.close()

print("         <---- Done.")

fig = plt.figure(figsize=(9.*cm, 9.*cm))
plt.hist(SrcE, bins=50, color='k', histtype='step')
plt.xlabel('$\\delta$ (m)')
plt.ylabel('Entries')
#plt.title('ExpSourceTst (' + today + \
#          '): $\\delta$ distribution', \
#          fontname="Times New Roman",  size=12)
plt.savefig('99-Scratch/SourceTst_d.pdf')
plt.close()

print("         <---- Done.")

##! Next: Tabulate parameters:
SourceTest += 1
print()
print("SourceTest:", SourceTest, \
      " test tabulation of parameters.")

print("     ----> Tabulate paramters:")
Src.tabulateParameters('99-Scratch/SourceTst_ParameterTable.tex')


##! Complete:
print()
print("========  Source: tests complete  ========")

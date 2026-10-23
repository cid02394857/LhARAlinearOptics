#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for handling unstable particles
===========================================

"""

import os
import shutil
import sys
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import BeamLine        as BL
import Particle        as Prtcl
import BeamLineElement as BLE

Debug = False

##! Start:
print("========  Decay-chain handling: tests start  ========")

print("decayCHAINTest:", "initialising:")

print("     ----> MatPlotLib options:")
#.. Use LaTeX for plot text only if it is installed:
mpl.rc('text', usetex=(shutil.which('latex') is not None))
mpl.rcParams['text.latex.preamble']="\\usepackage{bm}"
mpl.rcParams["figure.autolayout"]=True
mpl.rcParams['figure.constrained_layout.use'] = True

cm = 1./2.54  # centimeters in inches
print("     <---- MatPlotLib done.")

print(" <---- Done.")

##! Create test beam line:
decayCHAINTest = 0
print()
print("decayCHAINTest:", decayCHAINTest, " create beam line")
HOMEPATH = os.getenv('HOMEPATH')
filename = os.path.join(HOMEPATH, \
                        '11-Parameters/decayCHAIN.csv')
iBL  = BL.BeamLine(filename)
if Debug == True:
    print(iBL)
print("     ----> Beam line:", BLE.Facility.getinstances().getName())

##! Propagate default beam:
decayCHAINTest += 1
print()
print("decayCHAINTest:", decayCHAINTest, " propagate default beam")

iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
print("     ----> Reference particle species:", iRefPrtcl.getSpecies())

print("     ----> Generate events:")
iBL.trackBeam(10, None, None, None, False)
print("         <----", len(Prtcl.Particle.getinstances())-1, \
      "events generated.")

print("     ----> Analyse events:")
dfltMAXs = np.array([])
Cnt = 0
for iPrtcl in Prtcl.Particle.getinstances():
    if isinstance(iPrtcl, Prtcl.ReferenceParticle):
        if Debug: print("     ----> Skipping reference particle")
        continue
    Cnt += 1
    dfltMAXs = np.append(dfltMAXs, iPrtcl.gets()[-1])
    if Debug:
        print("         ----> n, Species, sMAX:", \
              Cnt, iPrtcl.getSpecies(), dfltMAXs[-1])

fig = plt.figure(figsize=(18.*cm, 6.*cm))
n, bins, patches = plt.hist(dfltMAXs, \
                            bins=101, range=(-0.5, 100.5), color='k', \
                            histtype='step', label='Default particle max s')


plt.xlim(-0.5, 101.5)
plt.yscale("log")
plt.savefig('99-Scratch/decayLENGTHdefault.pdf')
plt.close()
        
BL.BeamLine.cleaninstance()
BLE.BeamLineElement.cleaninstances()
Prtcl.Particle.cleanAllParticles()

##! Propagate pion beam
decayCHAINTest += 1
print()
print("decayCHAINTest:", decayCHAINTest, " propagate pion beam")

filename = os.path.join(HOMEPATH, \
                        '11-Parameters/decayCHAINpion.csv')
iBL  = BL.BeamLine(filename)
if Debug == True:
    print(iBL)
print("     ----> Beam line:", BLE.Facility.getinstances().getName())

iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
print("     ----> Reference particle species:", iRefPrtcl.getSpecies())

print("     ----> Generate events:")
iBL.trackBeam(10000, None, None, None, False)
print("         <----", len(Prtcl.Particle.getinstances())-1, \
      "events generated.")

print("     ----> Analyse events:")
pionMAXs = np.array([])
Cnt = 0
for iPrtcl in Prtcl.Particle.getinstances():
    if isinstance(iPrtcl, Prtcl.ReferenceParticle):
        if Debug: print("     ----> Skipping reference particle")
        continue
    Cnt += 1
    pionMAXs = np.append(pionMAXs, iPrtcl.gets()[-1])
    if Debug:
        print("         ----> n, Species, sMAX:", \
              Cnt, iPrtcl.getSpecies(), pionMAXs[-1])

fig = plt.figure(figsize=(18.*cm, 6.*cm))
n, bins, patches = plt.hist(pionMAXs, \
                            bins=101, range=(-0.5, 100.5), color='k', \
                            histtype='step', label='Default particle max s')


plt.xlim(-0.5, 101.5)
plt.yscale("log")
plt.savefig('99-Scratch/decayLENGTHpion.pdf')
plt.close()
        
BL.BeamLine.cleaninstance()
BLE.BeamLineElement.cleaninstances()
Prtcl.Particle.cleanAllParticles()
exit()

##! Test built-in methods:
decayCHAINTest = 1
print()
print("decayCHAINTest:", decayCHAINTest, \
      " check built-in methods.")

#.. __init__
print("    __init__:")
PrtcltInst = Prtcl.decayCHAIN("pion")

PrtcltInst.setLocation("Place 1")
PrtcltInst.setz(1.1)
PrtcltInst.sets(1.2)
TrcSpc = np.array([0.1, 0.002, 0.2, 0.004, 0., 18.])
PrtcltInst.setTraceSpace(TrcSpc)
PrtcltInst.setLocation("Place 2")
PrtcltInst.setz(2.1)
PrtcltInst.sets(2.2)
TrcSpc = np.array([0.15, 0.0025, 0.25, 0.0045, 0., 18.5])
PrtcltInst.setTraceSpace(TrcSpc)
#.. __repr__
print("    __repr__:")
print("      ---->", repr(PrtcltInst))
print("    <---- __repr__ done.")
#.. __str__
print("    __str__:")
print(PrtcltInst)
print("    <---- __str__ done.")


##! Get second particle - check decay time is different:
#   check when decay time for the particle is set.
decayCHAINTest = 2
print()
print("decayCHAINTest:", decayCHAINTest, " check if decay lengths are random")
PrtcltInst2 = Prtcl.decayCHAIN("pion")
PrtcltInst2.initremainingPath(3.0)
PrtcltInst3 = Prtcl.decayCHAIN("pion")
PrtcltInst3.initremainingPath(3.0)

PrtcltInst2.print()
PrtcltInst3.print()
#print(PrtcltInst2)
#print(PrtcltInst3)


#if (abs(PrtcltInst2.getremainingPath() - PrtcltInst3.getremainingPath()) < 0.0001):
#    print ("Path lengths of the two particles are too close")
#    print ("fail on test ", decayCHAINTest)
#    sys.exit(0)
#else:
#    print ("Path lengths of the two particles are different - test passed")

decayCHAINTest = 3
print()
print("decayCHAINTest:", decayCHAINTest, " check proton is stable and path infinite." )
ProtonInst = Prtcl.decayCHAIN("proton")
ProtonInst.initremainingPath(3.0)
#Prtcl.decayCHAIN.setDebug(True)
#ProtonInst.print()
print(ProtonInst)
#Prtcl.decayCHAIN.setDebug(False)

print()
decayCHAINTest = 4
print("decayCHAINTest:", decayCHAINTest, " check response to unknown particle type.")
try:
    uKInst = Prtcl.decayCHAIN("neutrino")
except:
    print ("correctly trapped non existant unstable particle - neutrino")

print()
print("========  decayCHAIN: tests complete  ========")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for "Octupole" class
===============================

Octupole.py -- set "relative" path to code

"""

import numpy as np
import os

import BeamLineElement as BLE
import BeamLine        as BL

HOMEPATH = os.getenv('HOMEPATH')
filename = os.path.join(HOMEPATH, \
                        '11-Parameters/Dummy4Tests.csv')

##! Start:
print("========  Octupole: tests start  ========")

##! Test singleton class feature:
OctupoleTest = 1
print()
print("OctupoleTest:", OctupoleTest, \
      " check built-in methods.")

#.. __init__
print("    __init__:")
try:
    Oct = BLE.Octupole()
except:
    print('      ----> Correctly trapped no argument exception.')
rStrt = np.array([0.,0.,0.])
vStrt = np.array([[np.pi/2.,np.pi/2.],[0.,0.]])
drStrt = np.array([0.,0.,0.])
dvStrt = np.array([0.,0.,0.])
try:
    Oct = BLE.Octupole("Octopole", rStrt, vStrt, drStrt, dvStrt)
except:
    print('      ----> Correctly trapped no octupole length exception.')

#--------> Clean instances and restart (facility needed for transport):
BLE.BeamLineElement.cleaninstances()
BLI  = BL.BeamLine(filename)

#.. Create valid instance:
Oct = BLE.Octupole("Octopole", rStrt, vStrt, drStrt, dvStrt, 1.5)
iRefPrtcl   = BL.BeamLine.getcurrentReferenceParticle()
refPrtclSet = iRefPrtcl.setReferenceParticle(Oct)
    
#.. __repr__
print("    __repr__:")
print("      ---->", repr(Oct))
print("    <---- __repr__ done.")
#.. __str__
print("    __str__:")
print(str(Oct))
print("    <---- __str__ done.")

##! Check get methods:
OctupoleTest = 3
print()
print("OctupoleTest:", OctupoleTest, " check get methods.")
print("    ----> print() method; tests all get methods")
print(Oct)

##! Check set method:
OctupoleTest = 4
print()
print("OctupoleTest:", OctupoleTest, " check set method.")
BLE.Octupole.setDebug(True)
print(Oct)
BLE.Octupole.setDebug(False)

##! Check set method:
OctupoleTest += 1
print()
print("OctupoleTest:", OctupoleTest, " test transport through octupole.")
R      = np.array([0.005, 0.001, -0.003, -0.002, 0., 0.])
Rprime = Oct.Transport(R)
print("     ----> Input phase-space vector:", R)
print("     ----> Relevant portions of transfer matrix:")
print("         ", Oct.getTransferMatrix()[0,0], \
                   Oct.getTransferMatrix()[0,1], \
                   Oct.getTransferMatrix()[0,2], \
                   Oct.getTransferMatrix()[0,3])
print("         ", Oct.getTransferMatrix()[1,0], \
                   Oct.getTransferMatrix()[1,1], \
                   Oct.getTransferMatrix()[1,2], \
                   Oct.getTransferMatrix()[1,3])
print("         ", Oct.getTransferMatrix()[2,0], \
                   Oct.getTransferMatrix()[2,1], \
                   Oct.getTransferMatrix()[2,2], \
                   Oct.getTransferMatrix()[2,3])
print("         ", Oct.getTransferMatrix()[3,0], \
                   Oct.getTransferMatrix()[3,1], \
                   Oct.getTransferMatrix()[3,2], \
                   Oct.getTransferMatrix()[3,3])
print("     ----> Transported phase-space vector:", Rprime)
#.. Linear optics: octupole acts as a drift of length 1.5 m:
RprimeTest = np.array([ 0.0065,  0.001, -0.006, -0.002,  0. ,  0.  ])
Diff       = np.subtract(Rprime, RprimeTest)
Norm       = np.linalg.norm(Diff)
print("     ----> Difference Rprime - RprimeTest:", Diff)
print("     ----> Magnitude of Diff:", Norm)
if Norm > 1E-6:
    raise Exception(" !!!!----> FAILED: octupole transport result not as expected.")
else:
    print(" <---- Octupole transport test successful.")

##! Complete:
print()
print("========  Octupole: tests complete  ========")

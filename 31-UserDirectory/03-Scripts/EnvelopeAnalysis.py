#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

import numpy as np

import Particle as Prtcl
import Beam     as Bm
import BeamLine as BL
import BeamIO   as bmIO
import EnvelopeOptimisation as EO
import BeamLineElement      as BLE
import UserFramework        as UsrFw


def main(argv):

    Success, Debug, \
        beamspecfile, inputfile, bdsimfile, outputfile, nEvts = \
        UsrFw.startAnalysis(argv)
    if not Success:
        print(" <---- Failed at UsrFw.startAnalysis, exit")
        exit(1)

    Success, ibmIOr, ibmIOw = UsrFw.handleFILES(beamspecfile, \
                                                inputfile, \
                                                outputfile, \
                                                bdsimfile)
    if not Success:
        print(" <---- Failed at UsrFw.handleFILES, exit")
        exit(1)

    #.. ----> Instanciate user analysis:
    iEO = EO.UserAnal(Debug)
        
    #.. ----> Instanciate extrapolate beam class and extrapolate:
    iexBm = Bm.extrapolateBeam(ibmIOr, nEvts, None, None)
    iexBm.extrapolateBeam()
    Prtcl.Particle.cleanParticles()
    
    ibmIOwStrt = bmIO.BeamIO(None, '99-Scratch/ParticlesStrt.dat', True)
    BL.BeamLine.getinstances().writeBeamLine(ibmIOwStrt.getdataFILE())

    nEvtGen = BL.BeamLine.getinstances().trackBeam(nEvts, \
                                        ibmIOwStrt.getdataFILE(),
                                        None, None, False)
    iexBm.plotBeamProgression('99-Scratch/BeamProgressStrt.pdf')
    
    #.. ----> End of event loop, wrap up:
    if Debug:
        print(" UserAnalysis: calling UsrAnal.UserEnd after event loop:")
    iEO.UserEnd()
    if Debug:
        print(" <---- Done.")

#--------  Got "out of the box" size of beam, now try and optimise:

    Debug = False
    if Debug:
        print(" EnvelopeAnalysis: start iterations.")
    for i in range(2000):
        if Debug:
            print("     ----> Iteration:", i)
        
        ibmIOr = iEO.setupIteration(beamspecfile, ibmIOr, iEO)

        iexBm.extrapolateBeam()
        
        if Debug:
            print("     New extrapolated beam:")
            iBm = Bm.Beam.getinstances()[0]
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ---->   sigmaxy[0]:", iBm.getsigmaxy()[0])
                print("     ----> emittance[0]:", iBm.getemittance()[0])
                print("     ---->     Twiss[0]:", iBm.getTwiss()[0])
                print("     ----> Updated Beam instance:", iBm)
        
        iEO.UserEnd()

#---- Recalculate beam size:

    ibmIOr   = iEO.setupIteration(beamspecfile, ibmIOr, iEO, False)
    dataFILE = ibmIOw.getdataFILE()
    
    BL.BeamLine.getinstances().writeBeamLine(dataFILE)
    print(BL.BeamLine.getinstances())

    nEvtGen  = BL.BeamLine.getinstances().trackBeam(nEvts,     \
                                                   dataFILE, \
                                                   None, None, False)

    iexBm.extrapolateBeam()
    iexBm.plotBeamProgression('99-Scratch/BeamProgressEnd.pdf')
    ibmIOw.flushNclosedataFile(dataFILE)
        
"""
    BL.BeamLine.setDebug(True)
   Execute main"
"""
if __name__ == "__main__":
   main(sys.argv[1:])

sys.exit(0)

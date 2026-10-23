#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
import cProfile
"""

import os
import sys, getopt
import struct
import math as mth

import Particle as Prtcl
import BeamLine as BL
import Beam     as Bm

def main(argv):
    """
       Parse input arguments:
    """
    opts, args = getopt.getopt(argv,"hdi:o:b:n:l:",\
                               ["ifile=","nEvts", "ofile=","bfile", "iLoc"])

    beamlinefile = None
    inputfile    = None
    outputfile   = None
    strtloc      = None
    Debug        = False
    nEvts        = None
    for opt, arg in opts:
        if opt == '-h':
            print ( \
                    'plotextrapolateBeam.py '  + \
                    ' -i <inputfile> -n <nEvts> -o <outputfile>' + \
                    ' -l <startlocation> [-b <beamlinefile>]')
            sys.exit()
        if opt == '-d':
            Debug = True
        elif opt in ("-b", "--bfile"):
            beamlinefile = arg
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-n", "--nEvts"):
            nEvts = int(arg)
        elif opt in ("-l", "--iLoc"):
            strtloc = int(arg)

    if inputfile    == None:
        print ( \
                'plotextrapolateBeam.py '  + \
                ' -i <inputfile> -n <nEvts> -o <outputfile>' + \
                ' -l <startlocation> [-b <beamlinefile>]')
        sys.exit()

    print(" plotextrapolateBeam: start")
    
    print("     ----> Initialise:")
    
    HOMEPATH    = os.getenv('HOMEPATH')
    print("         ----> HOMEPATH:", HOMEPATH)
        
    #.. Create beam instance:
    if beamlinefile != None:
        print("         ----> Create beam instance:")
        filename     = os.path.join(HOMEPATH, beamlinefile)
        print("             ----> Beamline parameters will be read from:", \
              filename)
        
    particlefile  = os.path.join(HOMEPATH, inputfile)
    print("             ----> Particles will be read from:", \
          particlefile)
    
    CSVoutputFILE = None
    if outputfile != None:
        CSVoutputFILE = os.path.join(HOMEPATH, outputfile)
        print("         ----> Write beamline summary file to:", CSVoutputFILE)

    print("         ----> Start at location:", strtloc)
    iexBm = Bm.extrapolateBeam(particlefile, nEvts, CSVoutputFILE, strtloc)

    print("     <---- Beam instance initialised.")

    print("     ----> Create report:")

    iexBm.extrapolateBeam()
    iexBm.createReport()

    print("     <---- Create report:")
        
    print("     ----> Plot beam progression:")

    iexBm.plotBeamProgression()

    print("     <---- Beam progression plot done.")
        
    print(" plotextrapolateBeam: ends")
    
"""
   Execute main"
"""
if __name__ == "__main__":
    """
    cProfile.run('main(sys.argv[1:])', '99-Scratch/restats')
    """
    main(sys.argv[1:])

sys.exit(0)

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
import BeamIO   as bmIO

def main(argv):
    """
       Parse input arguments:
    """
    opts, args = getopt.getopt(argv,"hdi:o:b:n:",\
                               ["ifile=","ofile=","bfile", "nEvts"])

    inputfile    = None
    beamlinesummaryfile   = None
    Debug        = False
    nEvts        = None
    for opt, arg in opts:
        if opt == '-h':
            print ( \
                'readBEAMsim.py -i <inputfile> -o <beamlinesummaryfile>' + \
                ' -n <nEvts>')
            sys.exit()
        if opt == '-d':
            Debug = True
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            beamlinesummaryfile = arg
        elif opt in ("-n", "--nEvts"):
            nEvts = int(arg)

    if inputfile    == None:
        print ( \
                'readBEAMsim.py -i <inputfile> -o <beamlinesummaryfile>' + \
                ' -n <nEvts>')
        sys.exit()

    print(" readBEAMsim: start")
    
    print("     ----> Initialise:")
    
    HOMEPATH    = os.getenv('HOMEPATH')
    print("         ----> HOMEPATH:", HOMEPATH)
        
    #.. File handling:
    print("         ----> Check input and beamlinesummary files:")

    #.. ----> Input file, if specified:
    if inputfile != None and not os.path.isfile(inputfile):
        inputfile = os.path.join(HOMEPATH, inputfile)
    if not os.path.isfile(inputfile): 
        print("             ----> Input file does not exist.")
        print("                   Exit.")
        sys.exit(1)

    ibmIOr = bmIO.BeamIO(None, inputfile)
    #ibmIOr = bmIO.BeamIO("99-Scratch", "Data4Tests.dat")
    print("             ----> Input file:", inputfile)
    
    if beamlinesummaryfile != None and not os.path.isabs(beamlinesummaryfile): 
        beamlinesummaryfile = os.path.join(HOMEPATH, beamlinesummaryfile)
    if beamlinesummaryfile != None and \
       not os.path.isdir(os.path.dirname(beamlinesummaryfile)):
        print("                 ----> Directory for beamlinesummary file", \
              os.path.dirname(beamlinesummaryfile), "does not exist.")
    else:
        print("             ----> Write beamline summary file to:",
              beamlinesummaryfile)
    #    print("                   Exit.")
    #    sys.exit(1)
        print("                   Beamlinesummary file not implemented.")
    
    print("     <---- Initialisation complete.")

    print("     ----> Read data file:")

    EndOfFile = False
    iEvt = 0
    iCnt = 0
    Scl  = 10
    print("         ----> Read data file:")
    while not EndOfFile:
        EndOfFile = ibmIOr.readBeamDataRecord()
    
        if not EndOfFile:
            iEvt += 1
            if (iEvt % Scl) == 0:
                print("         ----> Read event ", iEvt)
                iCnt += 1
                if iCnt == 10:
                    iCnt = 1
                    Scl  = Scl * 10
        if iEvt < 1:
            print(BL.BeamLine.getinstance())
            print(Prtcl.Particle.getinstances()[iEvt-1])
        if nEvts != None and iEvt == nEvts:
            break

    print("     <----", iEvt, "events read")

    print(" <---- Data-file reading done.")
        
    print(" ----> Plot progression:")
    
    Prtcl.Particle.plotTraceSpaceProgression()
    Prtcl.Particle.plotLongitudinalTraceSpaceProgression()
    
    print(" <---- Done.")

    print(" readBEAMsim: ends")
    
"""
   Execute main"
"""
if __name__ == "__main__":
    """
    cProfile.run('main(sys.argv[1:])', '99-Scratch/restats')
    """
    main(sys.argv[1:])

sys.exit(0)

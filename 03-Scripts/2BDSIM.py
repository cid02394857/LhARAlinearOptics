#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys, getopt
import struct
import math as mth

import Particle as Prtcl
import BeamLine as BL
import Beam     as Bm
import BeamIO   as bmIO

def main(argv):
    """
       Parse input arguments:
    """
    opts, args = getopt.getopt(argv,"hd:i:o:l:n:",\
                               ["ifile=","ofile=","loc", "nEvts"])

    inputfile    = None
    outputfile   = None
    Debug        = False
    iLoc         = 1
    nEvts        = None
    for opt, arg in opts:
        if opt == '-h':
            print ( \
                    '2BDSIM.py -i <inputfile> -o <outputfile>' + \
                    ' -l <start location> -n <nEvts>')
            sys.exit()
        if opt == '-d':
            Debug = True
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-l", "--loc"):
            iLoc = int(arg)
        elif opt in ("-n", "--nEvts"):
            nEvts = int(arg)

    if inputfile    == None:
        print ( \
                '2BDSIM.py -i <inputfile> -o <outputfile>' + \
                ' -l <start location> -n <nEvts>')
        sys.exit()

    print(" 2BDSIM: start")
    
    print("     ----> Initialise:")
    
    HOMEPATH    = os.getenv('HOMEPATH')
    print("         ----> HOMEPATH:", HOMEPATH)
        
    #.. File handling:
    print("         ----> Check input and output files:")

    #.. ----> Input file, if specified:
    if inputfile != None and not os.path.isfile(inputfile):
        inputfile = os.path.join(HOMEPATH, inputfile)
    if not os.path.isfile(inputfile): 
        print("             ----> Input file does not exist.")
        print("                   Exit.")
        sys.exit(1)

    ibmIOr = bmIO.BeamIO(None, inputfile)
    print("             ----> Input file:", inputfile)

    ibmIOw = bmIO.BeamIO(None, outputfile, True, True)
    print("             ----> Output file:", outputfile)
    
    print("     <---- Initialisation complete.")

    print("     ----> Read data file:")

    EndOfFile = False
    iEvt = -1
    iCnt = 0
    Scl  = 10
    print("         ----> Read data file:")
    while not EndOfFile:
        EndOfFile = ibmIOr.readBeamDataRecord()
        iEvt += 1
        if Debug and iEvt < 1:
            print(BL.BeamLine.getinstances()[0])
        if len(Prtcl.Particle.getinstances()) <= 1:
            continue
        if not EndOfFile:
            if (iEvt % Scl) == 0:
                print("         ----> Read event ", iEvt)
                iCnt += 1
                if iCnt == 10:
                    iCnt = 1
                    Scl  = Scl * 10

        iPrtcl = Prtcl.Particle.getinstances()[-1]
        iPrtcl.writeParticleBDSIM(ibmIOw.getdataFILE(), 1, True)
                    
        if nEvts != None and iEvt == nEvts:
            break

    print("     <----", iEvt, "events read")

    print(" <---- Data-file reading done.")
        
    print(" <---- Done.")

    print(" 2BDSIM: ends")
    
"""
   Execute main"
"""
if __name__ == "__main__":
   main(sys.argv[1:])

sys.exit(0)

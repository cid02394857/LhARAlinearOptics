#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
import cProfile
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import os
import sys, getopt
import numpy as np

import visualise as vis
import BeamLine  as BL
import Beam      as Bm
import BeamIO    as bmIO
import Particle  as Prtcl
import BeamLineElement as BLE

def main(argv):
    """
       Parse input arguments:
    """
    opts, args = getopt.getopt(argv,"hdi:o:b:n:",\
                               ["ifile=","ofile=","bfile", "nEvts"])

    beamlinefile = None
    inputfile    = None
    plotfile     = None
    Debug        = False
    nEvts        = 100
    for opt, arg in opts:
        if opt == '-h':
            print ( \
                    'visualiseBeam.py '  + \
                    ' -i <inputfile> -o <plotfile>' + \
                    ' -n <nEvts> [-b <beamlinefile>]')
            sys.exit()
        if opt == '-d':
            Debug = True
        elif opt in ("-b", "--bfile"):
            beamlinefile = arg
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            plotfile = arg
        elif opt in ("-n", "--nEvts"):
            nEvts = int(arg)

    if inputfile == None or \
       plotfile    == None:
        print ( \
                'visualiseBeam.py '  + \
                ' -i <inputfile> -o <plotfile>' + \
                ' -n <nEvts> [-b <beamlinefile>]')
        sys.exit()

    print(" visualiseBeam: start")
    
    print("     ----> Initialise:")
    
    HOMEPATH    = os.getenv('HOMEPATH')
    print("         ----> HOMEPATH:", HOMEPATH)
        
    #.. File handling:
    print("         ----> Check input and output files:")

    #.. ----> Input file, if specified:
    if inputfile != None and not os.path.isfile(inputfile):
        inputfile = os.path.join(HOMEPATH, inputfile)
    if inputfile != None and not os.path.isfile(inputfile): 
        print("             ----> Input file does not exist.")
        print("                   Exit.")
        sys.exit(1)

    ibmIOr = bmIO.BeamIO(None, inputfile)
    print("             ----> Input file:", inputfile)
    
    #.. ----> Input file, if specified:
    if not os.path.isabs(plotfile): 
        plotfile = os.path.join(HOMEPATH, plotfile)
    if not os.path.isdir(os.path.dirname(plotfile)):
        print("                 ----> Directory for plot file", \
              os.path.dirname(plotfile), "does not exist.")
        print("                   Exit.")
        sys.exit(1)

    print("             ----> Write plots to:", plotfile)
    
    print("     <---- Initialisation complete.")

    print("     ----> Visualisation:")

    #.. Read first record to set up geometry
    EndOfFile = False
    EndOfFile = ibmIOr.readBeamDataRecord()

    #.. Read particles to plot:
    for i in range(nEvts):
        ibmIOr.readBeamDataRecord()
    
    Prtcl.Particle.fillPhaseSpaceAll()
        
    ivisRPLCx = vis.visualise("RPLC", "xs")
    ivisRPLCy = vis.visualise("RPLC", "ys")

    ivisLabx = vis.visualise("Lab", "xz")
    ivisLaby = vis.visualise("Lab", "yz")
    
    font = {'family': 'serif', \
            'color':  'darkred', \
            }

    plt.rcParams["figure.figsize"] = (10., 7.5)

    with PdfPages(plotfile) as pdf:
        fig, axs = plt.subplots(nrows=2, ncols=1, \
                                layout="constrained")
        # add an artist, in this case a nice label in the middle...
        Facility = BL.BeamLine.getElement()[0].getName()
        Ttl = Facility + " (reference particle local coordinates)"
        fig.suptitle(Ttl, fontdict=font)
        axs[0].set_xlim(-0.5, 21.5)
        axs[0].set_ylim(-0.05,  0.05)
        axs[1].set_xlim(-0.5, 21.5)
        axs[1].set_ylim(-0.05,  0.05)

        ivisRPLCx.Particles(axs[0], nEvts)
        ivisRPLCx.BeamLine(axs[0])
        ivisRPLCy.Particles(axs[1], nEvts)    
        ivisRPLCy.BeamLine(axs[1])
    
        pdf.savefig()
        plt.close()
    
        fig, axs = plt.subplots(nrows=3, ncols=1, \
                                layout="constrained")
        gs     = axs[2].get_gridspec()
        axs[2].remove()
        gs     = axs[1].get_gridspec()
        axs[1].remove()
        axs[1] = fig.add_subplot(gs[1:])

        # add an artist, in this case a nice label in the middle...
        Ttl = Facility + " (laboratory reference frame)"
        fig.suptitle(Ttl, fontdict=font)

        axs[0].set_xlim(-0.5, 18.1)
        axs[0].set_ylim(-0.075,  0.075)
        axs[1].set_xlim(-0.5, 18.1)
        axs[1].set_ylim(-1.5,   5.9)

        ivisLabx.Particles(axs[0], nEvts)
        ivisLabx.BeamLine(axs[0])
        ivisLaby.Particles(axs[1], nEvts)
        ivisLaby.BeamLine(axs[1])
    
        pdf.savefig()
        plt.close()
    
    print("     <---- Visualisation done.")
        
    print(" visualiseBeam: ends")
    
"""
   Execute main"
"""
if __name__ == "__main__":
    """
    cProfile.run('main(sys.argv[1:])', '99-Scratch/restats')
    """
    main(sys.argv[1:])

sys.exit(0)

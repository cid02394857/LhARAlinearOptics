#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class Beam:
===========

  In some sense a "sister" class to Particle.  Whereas an instance of
  Particle records the passage of a particle travelling through the
  beam line, an instance of Beam records the collective properties of
  the beam such as emittance, etc., as the beam progresses through the
  beam line. 


  Class attributes:
  -----------------
    instances : List of instances of Particle class
  __Debug     : Debug flag

      
  Instance attributes:
  --------------------
   All instance attributes are initialised to Null

   Input arguments:
  _InputDataFile  : Path to BeamIO data file containing events to be read.
 _BeamLineInstance: Instance of BeamLine class to which the this instance of
                    the Beam class refers.
  _nEvtMax        : Maximum number of events to read, if not set, read 'em all
  _outputCSVfile  : Path to csv file in which summary of beam processing will
                    be written
  _startlocation  : Location at which extrapolation should start. Integer.
_beamlineSpecificationCSVfile : Kept for backward compatibility, optional.
                                Path to csv file containing specification of
                                beam line.

   _Location[] :   str : Name of location where parameters are recorded
   _s[]        : float : s coordinate at which parameters are recorded
   _nParticles :  list : Number of particles arriving at location
   _CovMtrx    :  list : Covariance matrix by location
   _sigmaxy    :  list : RMS x and y by location.  [0] sigmax, [1] sigmay
   _emittance  :  list : emittance by location; calculated from CovMtrx
                         [0] e_x, [1] e_y, [2] e_L, [3] e_{xy}, [4] e_{6D}
   _Twiss      :  list : Twiss parameters by location
                         [0][0] alpha_x, [0][1] beta_x, [0][2] gamma_x
                         [1][0] alpha_y, [1][1] beta_y, [1][2] gamma_y


  Note on indexing:
  -----------------
  When addressing the lists of BeamLineElements, PhaseSpaceRecords, and
  Bean parameters, the intention is to stick to the following local
  variable-name convention:

    iLoc        -- Address of location of BeamLineElement in
                   BLE.BeamLineElement.getinstances()
    iPhsSpcRcrd -- Address of locationof phase-space record in
                   Prtcl.Partiicle.getTraceSpace()
    iAddr       -- Address in list of CovMtrx etc.

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Set methods:
     setDebug: set class debug flag
           Input: bool, True/False
          Return: None

  resetBeamInsances:
               Sets cls.instances []

  setAll2None: Set all instance attributes to None.
        No input or return.

  setLocation: str : Set location string at which phase-space is stored.

         sets: float : Set s coordinate at which phase-space is stored.

  Other set methods believed to be self documenting:
   setbeamlineSpecificationCSVfile, setInputDataFile, setoutputCSVfile,
   setnEvtMax, setLocation, sets, setsigmaxy, setEmittance, setTwiss
   

  Get methods:
      getDebug, getbeamlineSpecificationCSVfile,getInputDataFile, 
      getoutputCSVfile, getBeamInstances(cls), getLocation, 
      getnEvtMax, getCovSums, getnParticles, getCovarianceMatrix,
      getsigmaxy, getemittance, getTwiss
          -- thought to be self documenting!

  Processing methods:
    cleanBeams : Deletes all Beam instances and resets list of
                 Beams.
         No input; Returns bool flag, True means all good.

    initialiseSums: initiaise sums used to calculate covariance matrix

     incrementSums: increment sums used to calculate covariance matrix
         Input: Instance of particle class

printProgression : prints evolution of beam parameters by location.

    createReport : creates csv file with evolutio of beam paramters by
                   location.

   initaliseSums : Initials sums needed to calculate covariance matrix

   incrementSums : Incrememnt sums
                input : iPrtcl : instance of particle class by which sums are
                                 to be updated

 calcCovarianceMatrix : Calculate covariance matrix given sums

   evaluateBeam : Read data file and increment sums.  Then, calculate
                  covariance matrix, RMS x and y, emittance, and Twiss
                  paramters
        Input : TrackBeam (Bool) : If reading particles from file,
                                       TrackBeam = False
                                   If needing to read particle at source
                                   and propagate through beam line,
                                       TrackBeam = True
                                  

plotBeamProgression : create standard plots
                      (in 99-Scratch/BeamProgressionPlot.pdf!) of beam
                      parameters.


  I/o methods:
    CSV file written out using Report module.

  Exceptions:
    badBeam, badParameter, badParameter, badTraceSpace


Created on Mon 28Feb24: Version history:
----------------------------------------
 1.0: 28Feb24: First implementation

@author: kennethlong
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from   copy   import deepcopy
import math   as     mth
import numpy  as     np
import os

import visualise         as vis
import Particle          as Prtcl
import BeamLine          as BL
import BeamLineElement   as BLE
import Report            as Rprt
import BeamIO            as bmIO


class Beam:
    instances  = []
    __Debug    = False


#--------  "Built-in methods":
    def __init__(self, _InputDataFile=None, _nEvtMax=None, \
                       _outputCSVfile=None, _startlocation=None, \
                       _beamlineSpecificationCSVfile=None):

        if self.getDebug():
            print(' Beam.__init__: \n', \
                  'creating the Beam object; input:', \
                  '     ----> InputDataFile, nEvtMax,', \
                  'outputCSVfile, startlocation,', \
                  'beamlineSpecificationCSVfile:', \
                  _InputDataFile, _nEvtMax, _outputCSVfile, _startlocation, \
                  _beamlineSpecificationCSVfile)
                  
        Beam.instances.append(self)

        self.setAll2None()

#--------  Check and initialise all inputs:  --------  --------  --------

        #.. Load parameter file
        #.. Check and open input data file
        if _InputDataFile                == None and \
           _beamlineSpecificationCSVfile == None:
            raise Exception( \
                " Beam.__init__: no input data file given ... and ...", \
                 "no csv file.")

        if _nEvtMax == None:
            pass
        elif isinstance(_nEvtMax, int):
            self.setnEvtMax(_nEvtMax)
        else:
            raise Exception(" Bad maximum number of events to read")

        #.. Check and output data file
        if _outputCSVfile == None:
            pass
        else:
            self.setoutputCSVfile(_outputCSVfile)
            dirname, filename = os.path.split(self.getoutputCSVfile())
            if not os.path.isdir(dirname):
                raise Exception( \
                        " Beam.__init__: output data frame invalid.")

        #.. Check and set start location:
        if _startlocation == None:
            pass
        else:
            self.setstartlocation(_startlocation)

#--------  <---- Check and initialise all inputs done.:  --------  --------

#--------  Open input data file, read first record, initialise sums:  -----

        if _InputDataFile != None:
            if isinstance(_InputDataFile, bmIO.BeamIO):
                _ibmIOr = _InputDataFile
            else:
                _ibmIOr = bmIO.BeamIO(None, _InputDataFile)
                
            self.setBeamIOread(_ibmIOr)
            ParticleFILE = self.getBeamIOread().getdataFILE()
            self.setInputDataFile(ParticleFILE)

            EndOfFile = False
            if BL.BeamLine.getinstances() == None:
                EndOfFile = self.getBeamIOread().readBeamDataRecord()

        iBm = BL.BeamLine.getinstances()
        
        if iBm == None:
            self.setbeamlineSpecificationCSVfile( \
                                        _beamlineSpecificationCSVfile)
            if self.getbeamlineSpecificationCSVfile() == None:
                raise badBeam(" No beam line read from i/p data file and " + \
                              "no beam-line specification csv file given.")
            iBm = BL.BeamLine(self.getbeamlineSpecificationCSVfile())

        #.. Set beam line instance:
        self.setBeamLineInstance(iBm)
        
        #.. Set locations and initialise sums:
        for iBLE in BLE.BeamLineElement.getinstances():
            if not isinstance(iBLE, BLE.Facility):
                self.setLocation(iBLE.getName())
        self.initialiseSums()
        
        #.. Must have reference particle:
        if not isinstance(BL.BeamLine.getcurrentReferenceParticle(), \
                          Prtcl.ReferenceParticle):
            raise noReferenceBeam(" Reference particle, ", \
                                      "not first in particle list.")

#--------  <---- Open input data file and load reference particle: done.  --
        if self.getDebug():
            print("     ----> New Beam instance: \n", \
                  Beam.__str__(self))
            print(" Beam.__init__: no maximum number of events requested,", \
                  " will read 'em all!")
            print(" <---- Beam instance created.")
            
    def __repr__(self):
        return "Beam(<InputDataFile>, " + \
               "nEvtMax=None, <OutputFile>=None, " + \
               "[<BeamLineSpecCSV>=None] )"

    def __str__(self):
        self.print()
        return " Beam __str__ done."

    def print(self):
        print(" Beam:")
        print(" -----")
        print("     ----> Debug flag:", self.getDebug())
        print("     ----> Beam specification file:", \
              self.getbeamlineSpecificationCSVfile())
        print("     ----> Input data file:", \
              self.getInputDataFile())
        print("     ----> Beam line:", \
              self.getBeamLineInstance().getElement()[0].getName())
        print("     ----> Number of events to read:", \
              self.getnEvtMax())
        print("     ----> Output data file:", \
              self.getoutputCSVfile())
        print("     ----> Start location:", \
              self.getstartlocation())
        print("     ----> Beam parameters by location:")

        print("         ----> start location, locations:",
              self.getstartlocation(), len(self.getLocation()))
        for iLoc in range(self.getstartlocation(), len(self.getLocation())+1):
            iAddr = iLoc - 1
            print("         ----> iLoc, iAddr:", iLoc, iAddr)
            print("         ----> iLoc:", iLoc, self.getLocation()[iAddr])
            if len(self.getnParticles()) > iAddr:
                print("             ----> Number of particles:", \
                      " nParticles:", int(self.getnParticles()[iAddr]))
            if len(self.getsigmaxy()) > iAddr:
                print("             ---->   sigma_x,   sigma_y:", \
                      self.getsigmaxy()[iAddr][0], \
                      self.getsigmaxy()[iAddr][1])
            if len(self.getemittance()) > iAddr:
                print("             ----> epsilon_x, epsilon_y:", \
                      self.getemittance()[iAddr][0], \
                      self.getemittance()[iAddr][1])
                print("             ----> epsilon_4, epsilon_l:", \
                      self.getemittance()[iAddr][2], \
                      self.getemittance()[iAddr][3])
                print("             ---->            epsilon_6:", \
                      self.getemittance()[iAddr][4])
            if len(self.getTwiss()) > iAddr:
                print("             ---->      Twiss paramters:", \
                      self.getTwiss()[iAddr])

        return " <---- Beam parameter dump complete."

    
#--------  "Set method" only Debug
#.. Method believed to be self documenting(!)

    @classmethod
    def setDebug(cls, Debug=False):
        if cls.__Debug:
            print(" Beam.setDebug: ", Debug)
        cls.__Debug = Debug

    @classmethod
    def resetBeamInstances(cls):
        if len(cls.instances) > 0:
            cls.instances = []
        
    def setAll2None(self):
        self._bmIOr                        = None
        self._InputDataFile                = None
        self._nEvtMax                      = None
        self._outputCSVfile                = None
        self._startlocation                = None
        self._beamlineSpecificationCSVfile = None
        self._BeamLineInstance                 = None

        self._Location   = []
        self._CovSums    = []
        self._nParticles = []
        self._CovMtrx    = []
        self._sigmaxy    = []
        self._emittance  = []
        self._Twiss      = []

    def setBeamLineInstance(self, _BeamLineInstance):
        self._BeamLineInstance = _BeamLineInstance
        
    def setbeamlineSpecificationCSVfile(self, _beamlineSpecificationCSVfile):
        self._beamlineSpecificationCSVfile = \
                        _beamlineSpecificationCSVfile

    def setBeamIOread(self, _bmIOr):
        self._bmIOr = _bmIOr

    def setInputDataFile(self, _InputDataFile):
        self._InputDataFile = _InputDataFile

    def setoutputCSVfile(self, _outputCSVfile):
        self._outputCSVfile = _outputCSVfile

    def setnEvtMax(self, _nEvtMax):
        if isinstance(_nEvtMax, int):
            self._nEvtMax = _nEvtMax
        else:
            raise Exception(" Bad maximum number of events to read")
        
    def setstartlocation(self, _startlocation):
        if isinstance(_startlocation, int):
            self._startlocation = _startlocation
        else:
            raise Exception(" Bad start location")
        
    def setLocation(self, Location):
        Success = False
        if isinstance(Location, str):
            self._Location.append(Location)
            Success = True
        return Success

    def sets(self, s):
        Success = False
        if isinstance(s, float):
            self._s.append(s)
            Success = True
        return Success

    def setsigmaxy(self):
        if self.getDebug():
            print(" Beam.setsigmaxy: start")

        for iAddr in range(len(self.getCovarianceMatrix())):
            if self.getDebug():
                print("     ----> sx2, sy2:", \
                      self.getCovarianceMatrix()[iAddr][0,0], \
                      self.getCovarianceMatrix()[iAddr][2,2])
            sx1  = mth.sqrt(self.getCovarianceMatrix()[iAddr][0,0])
            sy1  = mth.sqrt(self.getCovarianceMatrix()[iAddr][2,2])
                
            if self.getDebug():
                print("     <---- iAddr, sx1, sy1:", \
                      iAddr, sx1, sy1)
                
            self._sigmaxy.append([sx1, sy1])

        if len(self._sigmaxy) >= len(BLE.BeamLineElement.getinstances()):
            print(len(self._sigmaxy), len(BLE.BeamLineElement.getinstances()))
            raise BadCovMtrx(" too many sigmaxy!")
        if self.getDebug():
            print(" <---- Beam.setsigmaxy: done.")

    def setEmittance(self):
        if self.getDebug():
            print(" Beam.setEmittance: start")

        CovX = None
        CovY = None
        CovL = None
        Cov4 = None
        
        if len(self.getnParticles()) < 1 or \
           self.getnParticles()[0] < 10:
            return
        else:
            nPrtcls = self.getnParticles()[0] 
            
        for iAddr in range(len(self.getCovarianceMatrix())):
            if iAddr < len(self.getnParticles()):
                nPrtcls = self.getnParticles()[iAddr]
                
            if nPrtcls < 10: break
            
            if self.getDebug():
                print("     ----> iAddr:", iAddr)
                with np.printoptions(linewidth=500,precision=10,\
                                     suppress=True):
                    print("         CovMtrx: \n", \
                          self.getCovarianceMatrix()[iAddr]) 

            CovX = self.getCovarianceMatrix()[iAddr][0:2,0:2]
            CovY = self.getCovarianceMatrix()[iAddr][2:4,2:4]
            CovL = self.getCovarianceMatrix()[iAddr][4:6,4:6]
            Cov4 = self.getCovarianceMatrix()[iAddr][0:4,0:4]
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=10,\
                                     suppress=True):
                    print("                CovMtrx_X: \n", CovX)
                    print("                CovMtrx_Y: \n", CovY)
                    print("                CovMtrx_L: \n", CovL)
                    print("                CovMtrx_4: \n", Cov4)

            e2X  = np.linalg.det(CovX)
            e2Y  = np.linalg.det(CovY)
            e2L  = np.linalg.det(CovL)
            e24  = np.linalg.det(Cov4)
            e26  = np.linalg.det(self.getCovarianceMatrix()[iAddr])
            if self.getDebug():
                print("                e2X:", e2X)
                print("                e2Y:", e2Y)
                print("                e2L:", e2L)
                print("                e24:", e24)
                print("                e24:", e26)

            if e2X < 0.: e2X = 0.
            if e2Y < 0.: e2Y = 0.
            if e2L < 0.: e2L = 0.
            if e24 < 0.: e24 = 0.
            if e26 < 0.: e26 = 0.
            self._emittance.append([                \
                                     mth.sqrt(e2X), \
                                     mth.sqrt(e2Y), \
                                     mth.sqrt(e2L), \
                                     mth.pow(e24, 0.25), \
                                     mth.pow(e26, 1./6.)])
            
            if self.getDebug():
                print(" Beam.getEmittance:")
                print("     ----> CovX: \n", CovX)
                print("     ----> CovY: \n", CovY)
                print("     ----> CovL: \n", CovL)
                print("     ----> Cov4: \n", Cov4)
                print("     <---- eX:", self.getemittance()[iAddr][0])
                print("     <---- eY:", self.getemittance()[iAddr][1])
                print("     <---- eL:", self.getemittance()[iAddr][2])
                print("     <---- e4:", self.getemittance()[iAddr][3])
                print("     <---- e6:", self.getemittance()[iAddr][4])

        if self.getDebug():
            print(" <---- Beam.setEmittance: done.")

    def setTwiss(self):
        if self.getDebug():
            print(" Beam.setTwiss: start")
 
        if len(self.getnParticles()) < 1 or \
           self.getnParticles()[0] < 10:
            return
        else:
            nPrtcls = self.getnParticles()[0] 
            
        for iAddr in range(len(self.getCovarianceMatrix())):
            if iAddr < len(self.getnParticles()):
                nPrtcls = self.getnParticles()[iAddr]
                
            if nPrtcls < 10: break
            
            if self.getDebug():
                print("     ----> iAddr:", iAddr)
                with np.printoptions(linewidth=500,precision=10,\
                                     suppress=True):
                    print("         CovMtrx: \n", \
                          self.getCovarianceMatrix()[iAddr])
                print("     ----> ex, ey:", self.getemittance()[iAddr][0], \
                                            self.getemittance()[iAddr][1])
            emX = None
            emY = None
            if self.getemittance()[iAddr][0] > 0.:
                emX = self.getCovarianceMatrix()[iAddr][0:2,0:2] / \
                                       self.getemittance()[iAddr][0]
            if self.getemittance()[iAddr][1] > 0.:
                emY = self.getCovarianceMatrix()[iAddr][2:4,2:4] / \
                                       self.getemittance()[iAddr][1]

            if self.getemittance()[iAddr][0] > 0.:
                ax = float(-emX[0,1])
                bx =  float(emX[0,0])
                gx =  float(emX[1,1])
            else:
                ax = None
                bx = None
                gx = None

            if self.getemittance()[iAddr][1] > 0.:
                ay = float(-emY[0,1])
                by =  float(emY[0,0])
                gy =  float(emY[1,1])
            else:
                ay = None
                by = None
                gy = None

            self._Twiss.append([[ax, bx, gx], [ay, by, gy]])

            if self.getDebug():
                print("         ----> Twiss paramters, x, y:", \
                      self._Twiss[iAddr])

        if self.getDebug():
            print(" <---- Beam.setTwiss: done.")

                      
#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    @classmethod
    def getDebug(cls):
        return cls.__Debug

    def getBeamLineInstance(self):
        return self._BeamLineInstance
        
    def getbeamlineSpecificationCSVfile(self):
        return self._beamlineSpecificationCSVfile

    def getBeamIOread(self):
        return self._bmIOr
        
    def getInputDataFile(self):
        return self._InputDataFile
        
    def getoutputCSVfile(self):
        return self._outputCSVfile
        
    @classmethod
    def getinstances(cls):
        return cls.instances

    def getLocation(self):
        return self._Location
    
    def getnEvtMax(self):
        return self._nEvtMax

    def getstartlocation(self):
        if self._startlocation == None:
            startlocation = 1
        else:
            startlocation = self._startlocation
        return startlocation

    def getCovSums(self):
        return self._CovSums
    
    def getCovMtrx(self):
        return self._CovMtrx
    
    def getnParticles(self):
        return self._nParticles
    
    def getCovarianceMatrix(self):
        return self._CovMtrx

    def getsigmaxy(self):
        return self._sigmaxy

    def getemittance(self):
        return self._emittance

    def getTwiss(self):
        return self._Twiss

    
#--------  Processing methods:
    def initialiseSums(self):
        if self.getDebug():
            print(" Beam.intialiseSums start:")
            
        CovSums = np.zeros((6,6))

        self._CovSums    = []
        self._nParticles = []

        iLocMin = self.getstartlocation()
            
        for iLoc in range(iLocMin, \
                          len(BLE.BeamLineElement.getinstances())):
            
            iBLE = BLE.BeamLineElement.getinstances()[iLoc]
            if self.getDebug():
                print("     ----> iLoc, BLE name, type:", \
                      iLoc, iBLE.getName(), type(iBLE))

            self._CovSums.append(deepcopy(CovSums))
            self._nParticles.append(0.)
                
        if self.getDebug():
            print(" Beam.initialiseSums: n, CovSums:")
            for i in range(len(self.getnParticles())):
                with np.printoptions(linewidth=500,precision=7,
                                     suppress=True):
                    print("     ----> ", i, "\n", self._CovSums[i])
            
    def incrementSums(self, iPrtcl):
        startlocation = self.getstartlocation()
            
        if self.getDebug():
            print(" Beam.incrementSums start:")
            print("     ----> Start location:", startlocation)
            print("     ----> Number of locations:", \
                  len(self.getCovSums()))
            print("     ----> Particle trace space:")
            for iAddr in range(len(iPrtcl.getTraceSpace())):
                with np.printoptions(linewidth=500,precision=7,\
                                     suppress=True):
                    print("         ----> iAddr, trace space:", \
                          iAddr, \
                          iPrtcl.getTraceSpace()[iAddr])
            
        for iPhsSpcRcrd in range(startlocation-1, \
                                 len(iPrtcl.getTraceSpace())):
            iLoc  = iPhsSpcRcrd + 1
            iAddr = iLoc - startlocation
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=10,\
                                     suppress=True):
                    print("         ----> Phs, iLoc, iAddr:", \
                      iPhsSpcRcrd, iLoc, iAddr, \
                      iPrtcl.getTraceSpace()[iPhsSpcRcrd])
                    print("         ----> Location:", iLoc)
                    print("                CovSums_(i): \n", \
                          self.getCovSums()[iAddr])
                    print("         ----> iAddr, startlocation, ", \
                          "len(getnPartcles):", iAddr, startlocation, \
                          len(self.getnParticles()))
                    
            self._nParticles[iAddr] += 1

            for i in range(6):
                for j in range(i,6):
                    self._CovSums[iAddr][i,j]     = \
                        self._CovSums[iAddr][i,j] + \
                        iPrtcl.getTraceSpace()[iPhsSpcRcrd][i] * \
                        iPrtcl.getTraceSpace()[iPhsSpcRcrd][j]
                    if i != j:
                        self._CovSums[iAddr][j,i] = \
                            deepcopy(self._CovSums[iAddr][i,j])

            if self.getDebug():
                print("         ----> Location:", iLoc)
                with np.printoptions(linewidth=500,precision=7,\
                                     suppress=True):
                    print("                CovSums_(i+1): \n", \
                          self.getCovSums()[iAddr]) 
                                
        if self.getDebug():
            print(" <---- Beam.incrementSums: Done")

    def calcCovarianceMatrix(self):
        if self.getDebug():
            print(" Beam.calcCovarianceMatrix start:")
            print("     ----> Number of locations:", \
                  len(self.getCovSums()))

        for iAddr in range(len(self.getCovSums())):
            if self.getDebug():
                print("         ----> Location:", iAddr)
                print("               Number of particles:", \
                      self.getnParticles()[iAddr])
                with np.printoptions(linewidth=500,precision=7,\
                                     suppress=True):
                    print("                CovSums: \n", \
                          self.getCovSums()[iAddr])
                    
            if self.getnParticles()[iAddr] < 1:
                break
            
            self._CovMtrx.append(                             \
                                self.getCovSums()[iAddr] /         \
                                float(self.getnParticles()[iAddr]) \
                                 )
            
            if self.getDebug():
                print("     <---- Covariance matrix:")
                with np.printoptions(linewidth=500,precision=7, \
                                     suppress=True):
                    print("                CovMtrx: \n", \
                          self.getCovarianceMatrix()[iAddr]) 

    def evaluateBeam(self, TrackBeam=False):
        #.. Need to track beam if not reading particles from file.
        if self.getDebug():
            print(" Beam.evaluateBeam: ", \
                  "perform` sums to get covariance matrices")
            print("     ----> TrackBeam:", TrackBeam)
        
        EndOfFile = False
        iEvt = 0
        iCnt = 0
        Scl  = 1

        ParticleFILE = self.getInputDataFile()
        if self.getDebug():
            print("     ----: event loop")

        nEvtMax = self.getnEvtMax()
        if nEvtMax == None:
            nEvtMax = 1000
        iEvtStopClean = max(0, nEvtMax-1000)
        while not EndOfFile:
            EndOfFile = self.getBeamIOread().readBeamDataRecord()
            if not EndOfFile:
                iEvt += 1
                if (iEvt % Scl) == 0:
                    if self.getDebug():
                        print("         ----> Read event ", iEvt)
                    iCnt += 1
                    if iCnt == 10:
                        iCnt = 1
                        Scl  = Scl * 10

                iPrtcl = Prtcl.Particle.getinstances()[-1]
                if TrackBeam:
                    nEvtGen = BL.BeamLine.getinstances().trackBeam( \
                            1, None, iPrtcl, self.getstartlocation(), False)
                    
                self.incrementSums(iPrtcl)
                #.. Keep a few particles for plotting:
                Cleaned = False
                if iEvt < iEvtStopClean:
                    Cleaned = Prtcl.Particle.cleanParticles()
            
            if self.getDebug():
                print("     ----> Cleaned:", Cleaned)

            if self.getnEvtMax() != None and iEvt >= self.getnEvtMax():
                break
            
        if self.getDebug():
            print("     <----", iEvt, "events read")

        if self.getDebug():
            print("     ----> calculate covariance matrix:")
        self.calcCovarianceMatrix()
        if self.getDebug():
            print("     <---- done.")
            print("     ----> calculate sigma x, y:")
        self.setsigmaxy()
        if self.getDebug():
            print("     <---- done.")
            print("     ----> calculate emittances:")
        self.setEmittance()
        if self.getDebug():
            print("     <---- done.")
            print("     ----> Twiss paramters:")
        self.setTwiss()
        if self.getDebug():
            print("     <---- done.")
        
        
#--------  Utilities:
    @classmethod
    def cleanBeams(cls):
        DoneOK = False
        
        for iBm in cls.getinstances():
            del iBm
            
        cls.resetBeamInstances()
        DoneOK = True

        return DoneOK
    
    def printProgression(self):
        for iAddr in range(len(self.getLocation())):
            with np.printoptions(linewidth=500,precision=5, \
                                 suppress=True):
                print(self.getLocation()[iAddr], \
                      ": z, s, trace space:", \
                      self.getz()[iAddr], self.gets()[iAddr])

    def getHeader(self):
        HeaderList = ["Name", "nPrtcls", "sigmaxy", "emittance", "Twiss"]
        return HeaderList

    def getLines(self):
        DataList = []
        iLocMin = self.getstartlocation()

        if self.getDebug():
            print(" Beam.getLines:")
            print("     ----> Start location:", self.getstartlocation())
            print("     ----> Number of locations:", \
                  len(self.getLocation()))
            print("     ----> Start location:", iLocMin)
        
        for iLoc in range(iLocMin, \
                          len(BLE.BeamLineElement.getinstances())):
            iAddr = iLoc - iLocMin
            if self.getDebug():
                print("         ----> iLoc, Name, iAddr:", \
                      iLoc, self.getLocation()[iLoc-1], iAddr)
                print("         ----> len: location, nParticles, ", \
                      "sigmaxy, emittance, Twiss:", \
                      len(self.getLocation()), \
                      len(self.getnParticles()), \
                      len(self.getsigmaxy()), \
                      len(self.getemittance()), \
                      len(self.getTwiss()))

            nPrtcls = None
            if iAddr < len(self.getnParticles()):
                nPrtcls = self.getnParticles()[iAddr]
            else:
                nPrtcls = None
                
            if iAddr < len(self.getemittance()):
                DataList.append([  \
                                   self.getLocation()[iLoc-1], \
                                   nPrtcls, \
                                   self.getsigmaxy()[iAddr], \
                                   self.getemittance()[iAddr], \
                                   self.getTwiss()[iAddr] \
                                 ])
            elif iAddr < len(self.getsigmaxy()):
                DataList.append([  \
                                   self.getLocation()[iLoc-1], \
                                   nPrtcls, \
                                   self.getsigmaxy()[iAddr], \
                                   None, \
                                   None  \
                                 ])
            else:
                DataList.append([  \
                                   self.getLocation()[iLoc-1], \
                                   None, \
                                   None, \
                                   None, \
                                   None, \
                                 ])

        return DataList
    
    def createReport(self):
        if self.getoutputCSVfile() == None:
            print(" Beam.createReport: no data file given, skip.")
        else:
            Name   = BLE.BeamLineElement.getinstances()[0].getName()
            iRprt  = Rprt.Report(Name, None, self.getoutputCSVfile(),
                             self.getHeader(), self.getLines())
            if self.getDebug():
                print("Beam.createReport:")
                print(iRprt)

            iRprt.asCSV()

    def plotBeamProgression(self, \
                            plotFILE='99-Scratch/BeamProgressionPlot.pdf'):
        pathNAME = os.path.split(plotFILE)
        if not os.path.exists(pathNAME[0]):
            raise noPath4plotFILE( \
                          " Beam.plotBeamProgression:", \
                          " path for plotFILE does not exist!")

        if self.getDebug():
            print(" Beam.plotBeamProgression: start")

        font = {'family': 'serif', \
                'color':  'darkred', \
                }
        plt.rcParams["figure.figsize"] = (10., 10.0)

        iLocMin = self.getstartlocation()
        nLocs = len(BLE.BeamLineElement.getinstances()) - iLocMin
    
        s     = [None]*nLocs
        sx    = [None]*nLocs
        sy    = [None]*nLocs

        ex    = [None]*nLocs
        ey    = [None]*nLocs
        exy   = [None]*nLocs

        ax    = [None]*nLocs
        ay    = [None]*nLocs
        bx    = [None]*nLocs
        by    = [None]*nLocs
        gx    = [None]*nLocs
        gy    = [None]*nLocs
        
        """
        sx    = []
        sy    = []
        
        ex    = []
        ey    = []
        exy   = []

        ax    = []
        ay    = []
        bx    = []
        by    = []
        gx    = []
        gy    = []
        """

        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()

        if self.getDebug():
            print(BL.BeamLine.getinstances())
            print(self)
        
        if self.getDebug():
            print("     ----> Start location:", self.getstartlocation())
            print("     ----> Number of locations:", \
                  len(self.getLocation()))
        
        for iLoc in range(iLocMin, \
                          len(BLE.BeamLineElement.getinstances())):
            iAddr = iLoc - iLocMin
            if self.getDebug():
                print("         ----> start location, locations:",
                      self.getstartlocation(), len(self.getLocation()))
                print("         ----> iLoc, iAddr:", iLoc, iAddr)

                  
            s[iAddr] = iRefPrtcl.getsOut()[iLoc-1]
            if iAddr < len(self.getsigmaxy()):
                sx[iAddr] = self.getsigmaxy()[iAddr][0]
                sy[iAddr] = self.getsigmaxy()[iAddr][1]
                if iAddr < len(self.getemittance()):
                    ex[iAddr] = self.getemittance()[iAddr][0]
                    ey[iAddr] = self.getemittance()[iAddr][1]
                    exy[iAddr] = self.getemittance()[iAddr][3]

                    bx[iAddr] = self.getTwiss()[iAddr][0][1]
                    by[iAddr] = self.getTwiss()[iAddr][1][1]
                    ax[iAddr] = self.getTwiss()[iAddr][0][0]
                    ay[iAddr] = self.getTwiss()[iAddr][1][0]
                    gx[iAddr] = self.getTwiss()[iAddr][0][2]
                    gy[iAddr] = self.getTwiss()[iAddr][1][2]
            
            if self.getDebug():
                if iAddr < len(self.getsigmaxy()):
                    print("     ----> iLoc, s, sx, sy:", \
                          iLoc, self.getLocation()[iAddr], \
                          s[iAddr], sx[iAddr], sy[iAddr])

        with PdfPages(plotFILE) as pdf:
            fig, axs = plt.subplots(nrows=6, ncols=1, \
                                    layout="constrained")
            # add an artist, in this case a nice label in the middle...
            Ttl = "Test"
            fig.suptitle(Ttl, fontdict=font)

            ivisRPLCx = vis.visualise("RPLC", "xs")
            axs[0].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[0].set_ylim(-0.05, 0.05)
            ivisRPLCx.Particles(axs[0], 1000)
            ivisRPLCx.BeamLine(axs[0])
        
            ivisRPLCy = vis.visualise("RPLC", "ys")
            axs[1].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[1].set_ylim(-0.05, 0.05)
            ivisRPLCy.Particles(axs[1], 1000)
            ivisRPLCy.BeamLine(axs[1])
        
            gs     = axs[3].get_gridspec()
            axs[3].remove()
            gs     = axs[2].get_gridspec()
            axs[2].remove()
            
            axs[2] = fig.add_subplot(gs[2:4])

            gs     = axs[5].get_gridspec()
            axs[5].remove()
            gs     = axs[3].get_gridspec()
            axs[4].remove()
            axs[4] = fig.add_subplot(gs[4:])

            axs[2].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[2].plot(s[0:len(sx)], sx, \
                color='b', marker='o', markersize=4, label='s_x')
            axs[2].plot(s[0:len(sy)], sy, color='r', marker='s', markersize=4, \
                        linestyle='dashed', label='s_y')
            axs[2].legend()
            axs[2].set_xlabel('s (m)')
            axs[2].set_ylabel('s_{xy} (m)')

            axs[4].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[4].plot(s[0:len(ex)], ex, \
                          color='b', marker='o', markersize=4, label='e_x')
            axs[4].plot(s[0:len(ey)], ey, \
                          color='r', marker='s', markersize=4, 
                          linestyle='dashed', label='e_y')
            axs[4].legend(loc='upper left')
            axs[4].set_xlabel('s (m)')
            axs[4].set_ylabel('e_{xy} (m)')

            ax2 = axs[4].twinx()
            ax2.plot(s[0:len(exy)], exy, \
                       color='g', marker='d', markersize=4, \
                       linestyle='dotted', label='e_(xy)')
            ax2.legend()
            ax2.set_ylabel('e_(xy) (m)')

            pdf.savefig()
            plt.close()

            # Twiss parameters and dispersion
            
            fig, axs = plt.subplots(nrows=8, ncols=1, \
                                    layout="constrained")
            # add an artist, in this case a nice label in the middle...
            Ttl = "Test"
            fig.suptitle(Ttl, fontdict=font)

            ivisRPLCx = vis.visualise("RPLC", "xs")
            axs[0].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[0].set_ylim(-0.05, 0.05)
            ivisRPLCx.Particles(axs[0], 1000)
            ivisRPLCx.BeamLine(axs[0])
        
            ivisRPLCy = vis.visualise("RPLC", "ys")
            axs[1].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[1].set_ylim(-0.05, 0.05)
            ivisRPLCy.Particles(axs[1], 1000)
            ivisRPLCy.BeamLine(axs[1])
        
            gs     = axs[3].get_gridspec()
            axs[3].remove()
            gs     = axs[2].get_gridspec()
            axs[2].remove()
            
            axs[2] = fig.add_subplot(gs[2:4])

            gs     = axs[5].get_gridspec()
            axs[5].remove()
            gs     = axs[4].get_gridspec()
            axs[4].remove()
            
            axs[3] = fig.add_subplot(gs[4:6])
    
            gs     = axs[7].get_gridspec()
            axs[7].remove()
            gs     = axs[6].get_gridspec()
            axs[6].remove()
            
            axs[4] = fig.add_subplot(gs[6:])
    
            axs[2].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[2].plot(s[0:len(bx)], bx, \
                          color='b', marker='o', markersize=4, \
                          label='b_x')
            axs[2].plot(s[0:len(by)], by, \
                          color='r', marker='s', markersize=4, \
                          linestyle='dashed', label='b_y')
            axs[2].legend(loc='upper center')
            axs[2].set_xlabel('s (m)')
            axs[2].set_ylabel('b_{xy} (m)')

            ax1 = axs[2].twinx()
            ax1.plot(s[0:len(exy)], exy, \
                     color='g', marker='d', markersize=4, \
                     linestyle='dotted', label='e_(xy)')
            ax1.legend()
            ax1.set_ylabel('e_(xy) (m)')

            axs[3].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[3].plot(s[0:len(ax)], ax, \
                          color='b', marker='o', markersize=4, \
                          label='a_x')
            axs[3].plot(s[0:len(ay)], ay, \
                        color='r', marker='s', markersize=4, 
                        linestyle='dashed', label='a_y')
            axs[3].legend()
            axs[3].set_xlabel('s (m)')
            axs[3].set_ylabel('a_{xy}')

            axs[4].set_xlim(-0.5, \
                BL.BeamLine.getcurrentReferenceParticle().gets()[-1]+0.5)
            axs[4].plot(s[0:len(gx)], gx, \
                          color='b', marker='o', markersize=4, \
                          label='g_x')
            axs[4].plot(s[0:len(gy)], gy, \
                        color='r', marker='s', markersize=4, \
                        linestyle="dashed", label='g_y')
            axs[4].legend()
            axs[4].set_xlabel('s (m)')
            axs[4].set_ylabel('g_{xy}')

            pdf.savefig()
            plt.close()
            
        if self.getDebug():
            print(" <----  Beam.plotBeamProgression done.")

                
"""
Derived class extrapolateBeam(Beam):
=====================================

  The Beam class provides methods to calculate the collective
  properties of the beam such as emittance, etc., from the
  distributions of beam particles as the beam progresses through the 
  beam line.  The extrapolateBeam derived class uses the transfer
  matrices to extrapolate the initial phase-space distribution along
  the beam line.


  Class attributes:
  -----------------
    instances : List of instances of Particle class
  __Debug     : Debug flag

      
  Instance attributes:
  --------------------
   All instance attributes are initialised to Null

   Input arguments:
  _InputDataFile  : Path to BeamIO data file containing events to be read.
  _nEvtMax        : Maximum number of events to read, if not set, read 'em all
   _outputCSVfile : Path to csv file in which summary of beam processing will
                    be written

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Set methods:
     setDebug: set class debug flag
           Input: bool, True/False
          Return: None

  resetBeamInsances:
               Sets cls.instances []

  setAll2None: Set all instance attributes to None.
        No input or return.


  Get methods:
      getDebug, getbeamlineSpecificationCSVfile,getInputDataFile, 
      getoutputCSVfile, getextrapolateBeamInstances
          -- thought to be self documenting!

  Processing methods:
    cleanextrapolateBeams : Deletes all extrapolateBeam instances
                             and resets list of extrapolateBeam.
         No input; Returns bool flag, True means all good.


  I/o methods:
    CSV file written out using Report module.

  Exceptions:
    badBeam, badParameter, badParameter, badTraceSpace


Created on Mon 28Feb24: Version history:
----------------------------------------
 1.0: 09Apr24: First implementation

@author: kennethlong
"""

class extrapolateBeam(Beam):
    instances  = []
    __Debug    = False


#--------  "Built-in methods":
    def __init__(self, _InputDataFile=None, _nEvtMax=None, \
                       _outputCSVfile=None, _startlocation=None,\
                       _beamlineSpecificationCSVfile=None):
        
        if self.__Debug:
            print(' extrapolateBeam.__init__: ', \
                  'creating the Beam object:', \
                  'start location =', _startlocation)
            
        extrapolateBeam.instances.append(self)

#--------  Check and initialise all inputs:  --------  --------  --------

        Beam.__init__(self, _InputDataFile, _nEvtMax, _outputCSVfile, \
                      _startlocation, _beamlineSpecificationCSVfile)

#--------  <---- Check and initialise all inputs done.:  --------  --------

        if self.__Debug:
            print("     ----> New extrapolateBeam instance: \n", \
                  extrapolateBeam.__str__(self))
            print(" <---- Beam instance created.")
            
    def __repr__(self):
        return "extrapolateBeam(<BeamLineSpecCSV>, <InputDataFile>, " + \
               "nEvtMax=None, startlocation=None <OutputFile>=None)"

    def __str__(self):
        print(" extrapolateBeam:")
        print(" -----------------")
        print("     ----> Start location:", self.getstartlocation())
        self.print()
        return " extrapolateBeam __str__ done."


#--------  "Set method" only Debug
#.. Method believed to be self documenting(!)

    @classmethod
    def setDebug(cls, Debug=False):
        if cls.__Debug:
            print(" extrapolateBeam.setdebug: ", Debug)
        cls.__Debug = Debug

    @classmethod
    def resetextrapolateBeamInstances(cls):
        if len(cls.instances) > 0:
            cls.instances = []
        
    def setAll2None(self):
        Beam.setAll2None(self)
        
    
#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    @classmethod
    def getDebug(cls):
        return cls.__Debug

    @classmethod
    def getextrapolateBeamInstances(cls):
        return cls.instances

#--------  Processing methods:    
    def initialiseSums(self):
        if self.getDebug():
            print(" extrapolateBeam.intialiseSums start:")
            
        CovSums = np.zeros((6,6))

        iLoc = self.getstartlocation()
            
        iBLE = BLE.BeamLineElement.getinstances()[iLoc]
        if self.getDebug():
            print("     ----> iLoc, BLE name, type:", \
                  iLoc, iBLE.getName(), type(iBLE))

        self._CovSums.append(deepcopy(CovSums))
        self._nParticles.append(0.)
                
        if self.getDebug():
            print("     ----> n, CovSums:")
            for i in range(len(self.getnParticles())):
                with np.printoptions(linewidth=500,precision=7,
                                     suppress=True):
                    print("         ----> n =", i, "\n", self._CovSums[i])

    def incrementSums(self, iPrtcl):
        startlocation = self.getstartlocation()
            
        if self.getDebug():
            print(" extrapolateBeam.incrementSums start:")
            print("     ----> Number of locations:", \
                  len(self.getCovSums()))
            print("     ----> Start location", startlocation, \
                  BLE.BeamLineElement.getinstances()[startlocation].getName())
            print("     ----> Particle trace space:")
            print("         ----> Number of trace-space vectors:", \
                  len(iPrtcl.getTraceSpace()))
            if len(iPrtcl.getTraceSpace()) > startlocation:
                with np.printoptions(linewidth=500,precision=7,\
                                     suppress=True):
                    print("         ----> iLoc, trace space:", \
                          startlocation, \
                          iPrtcl.getTraceSpace()[startlocation-1])

        if len(iPrtcl.getTraceSpace()) >= startlocation:
            
            iAddr = 0
            if self.getDebug():
                print("     ----> Start location:", startlocation, \
                BLE.BeamLineElement.getinstances()[startlocation].getName())
            
            self._nParticles[iAddr] += 1

            if self.getDebug():
                print("         ----> Location:", startlocation)
                with np.printoptions(linewidth=500,precision=10,\
                                     suppress=True):
                    print("                CovSums_(i): \n", \
                          self.getCovSums()[iAddr])

            for i in range(6):
                for j in range(i,6):
                    self._CovSums[iAddr][i,j]     = \
                        self._CovSums[iAddr][i,j] + \
                        iPrtcl.getTraceSpace()[startlocation-1][i] * \
                        iPrtcl.getTraceSpace()[startlocation-1][j]
                    if i != j:
                        self._CovSums[iAddr][j,i] = \
                            deepcopy(self._CovSums[iAddr][i,j])
            if self.getDebug():
                print("         ----> Location:", startlocation)
                with np.printoptions(linewidth=500,precision=7,\
                                     suppress=True):
                    print("                CovSums_(i+1): \n", \
                          self.getCovSums()[iAddr]) 

        if self.getDebug():
            print(" <---- extrapolateBeam.incrementSums: Done")

    def extrapolateCovarianceMatrix(self):
        if self.getDebug():
            print(" extrapolateBeam.extrapolateCovarianceMatrix start:")

        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()

        iLocMin = self.getstartlocation()
        if self.getDebug():
            print("     ----> iLocMin:", iLocMin)
        
        for jLoc in range(iLocMin+1, \
                          len(BLE.BeamLineElement.getinstances())):
            
            jAddr       = jLoc - iLocMin 
            iLoc        = jLoc - 1
            iAddr       = iLoc - iLocMin 
            iPhsSpcRcrd = iLoc - 1
            if self.getDebug():
                print("         ---->           Location:", jLoc)
                print("         ---->  Previous location:", iLoc)
                print("         ----> Trace-space record:", iPhsSpcRcrd)
                print("         ---->              jAddr:", jAddr)
                print("         ---->              iAddr:", iAddr)
                          
            jBLE  = BLE.BeamLineElement.getinstances()[jLoc]
            iBLE  = BLE.BeamLineElement.getinstances()[iLoc]
            
            if self.getDebug():
                print("     ----> iBLE name, type:", \
                      iBLE.getName(), type(iBLE))
                print("     ----> jBLE name, type:", \
                      jBLE.getName(), type(jBLE))

            if isinstance(jBLE, BLE.Drift) or \
               isinstance(jBLE, BLE.Aperture) or \
               isinstance(jBLE, BLE.Octupole) or \
               isinstance(jBLE, BLE.CylindricalRFCavity) or \
               isinstance(jBLE, BLE.RPLCswitch):
                pass
            else:
                TrcSpc       = iRefPrtcl.getTraceSpace()[iPhsSpcRcrd]
                TrnsfrMtrx   = jBLE.setTransferMatrix(TrcSpc)
                
            TrnsfrMtrx       = jBLE.getTransferMatrix()
            TrnspsTrnsfrMtrx = np.transpose(TrnsfrMtrx)
                
            if self.getDebug():
                print("         ----> Name:", jBLE.getName())
                with np.printoptions(linewidth=500,precision=10,\
                                     suppress=True):
                    print("         ----> Transfer matrix: \n", \
                          TrnsfrMtrx)
                    print("         ----> Transpose of transfer matrix: \n",\
                          TrnspsTrnsfrMtrx)
                    print("         ----> iLoc, covariance matrix:", \
                          iLoc, " \n", \
                          self.getCovMtrx()[iAddr])

            CovInv  = np.matmul(self.getCovMtrx()[iAddr], TrnspsTrnsfrMtrx)
            CovMtrx = np.matmul(TrnsfrMtrx, CovInv)

            self._CovMtrx.append(CovMtrx)
            
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=10, \
                                     suppress=True):
                    print("         ----> jLoc, CovMtrx: \n", \
                          self.getCovMtrx()[jAddr])

    def extrapolateBeam(self):
        if self.getDebug():
            print(" Beam.extrapolateBeam: transport", \
                  "beam envelope")
            print("     ----> BeamLine: id(BL):", \
                  id(BL.BeamLine.getinstances()))
            print(BL.BeamLine.getinstances())
        
        EndOfFile = False
        iEvt = 0
        iCnt = 0
        Scl  = 1

        ParticleFILE = self.getInputDataFile()

        Generate  = False
        GotCvMtrx = False
            
        #.. if ParticleFILE is None, assume will generate starting
        #   covariance matrix
        if self.getDebug():
            print("     ----> Choice of initial covariance matrix:")
            print("         ----> len(self.getCovMtrx():", \
                  len(self.getCovMtrx()))
            print("         ----> ParticleFILE.closed:", ParticleFILE.closed)
            
        if len(self.getCovMtrx()) == 1:
            if self.getDebug():
                print("     ----> Got covariance matrix.")
            GotCvMtrx = True
        
        elif ParticleFILE == None:
            if self.getDebug():
                print("     ----> No particle file closed, so generate", \
                      "initial covariance matrix as specified in csv file.")
            Generate = True
            
        #.. if ParticleFILE is closed, assume dont need to make initial
        #   covariance matrix
        elif ParticleFILE.closed:
            if self.getDebug():
                print("     ----> Particle file closed, so continue from", \
                      " stored source covariance matrix.")
            GotCvMtrx = True

        if not GotCvMtrx:
            if self.getDebug():
                print("     ----> Event loop")

            nEvtMax = self.getnEvtMax()
            if nEvtMax == None:
                nEvtMax = 1000
            iEvtStopClean = max(0, nEvtMax-1000)
            Cleaned       = None
            while not EndOfFile:
                if not Generate:
                    EndOfFile = self.getBeamIOread().readBeamDataRecord()

                if not EndOfFile:
                    iPrtcl = Prtcl.Particle.getinstances()[-1]
                    rEvt = BL.BeamLine.trackBeam(1, None, iPrtcl)

                    iEvt += 1
                    if (iEvt % Scl) == 0:
                        if self.getDebug():
                            print("         ----> Read event ", iEvt)
                        iCnt += 1
                        if iCnt == 10:
                            iCnt = 1
                            Scl  = Scl * 10

                    iPrtcl = Prtcl.Particle.getinstances()[-1]
                    self.incrementSums(iPrtcl)

                    #.. Keep a few particles for plotting:
                    if iEvt < iEvtStopClean:
                        Cleaned = Prtcl.Particle.cleanParticles()
            
                    if self.getDebug():
                        print("     ----> Cleaned:", Cleaned)

                    if self.getnEvtMax() != None and iEvt >= \
                       self.getnEvtMax():
                        break
            if self.getDebug():
                print("     <----", iEvt, "events read")

            if self.getDebug():
                print("     ----> calculate covariance matrix:")
            self.calcCovarianceMatrix()
                
        if self.getDebug():
            print("     <---- done.")
            print("     ----> extrapolate covariance matrix:")
        self.extrapolateCovarianceMatrix()
        if self.getDebug():
            print("     <---- done.")
            print("     ----> calculate sigma x, y:")
        self.setsigmaxy()
        if self.getDebug():
            print("     <---- done.")
            print("     ----> calculate emittances:")
        self.setEmittance()
        if self.getDebug():
            print("     <---- done.")
            print("     ----> Twiss paramters:")
        self.setTwiss()
        if self.getDebug():
            print("     <---- done.")
        

#--------  Utilities:
    @classmethod
    def cleanextrapolateBeams(cls):
        DoneOK = False
        
        for iexBm in cls.getextrapolateBeamInstances():
            del iexBm
            
        cls.resetextrapolateBeamInstances()
        DoneOK = True

        return DoneOK

    
#--------  Exceptions:
class noReferenceBeam(Exception):
    pass

class badBeam(Exception):
    pass

class badParameter(Exception):
    pass

class badTraceSpace(Exception):
    pass

class noPath4plotFILE(Exception):
    pass

class BadCovMtrx(Exception):
    pass

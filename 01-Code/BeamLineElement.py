#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class BeamLineElement:
======================

  Class to contain beam line elements.  Parent class defines principal
  coordinates of the beam-line element and transfer matrix.  Properties
  of the beam-line elements are derived classes.

  Classes derived from BeamLineElement:
             Drift, Aperture, FocusQuadrupole, DeFocusQuadrupole,
             SectorDipole, Octupole, Solenoid, CylindricalRFCavity, Source


  Class attributes:
  -----------------
        instances : List of instances of BeamLineElement class
      __Debug     : Debug flag
constants_instance: Instance of PhysicalConstants class
    speed_of_light: Speed of light from PhysicalConstants

      
  Instance attributes:
  --------------------
  Calling arguments:
       _Name : string; name of element, should be identifiable, e.g.,
               Drift1 
       _rStrt : numpy array; x, y, z position (in m) of start of
                element. 
       _vStrt : numpy array; theta, phi of principal axis of element.
      _drStrt : "error", displacement of start from nominal position
                (rad). 
      _dvStrt : "error", deviation in theta and phy from nominal axis
                in form of Euler angles (rad).  Convention is an
                "intrinsic" rotation in the "z-y-z" convention.

  Calculated and set in derived classes:
    _Strt2End : Calculated; vector that translates from start to end of
                element in lab coordinates.  Set in derived class.
  _Rot2LbStrt : Calculated; rotation matrix that takes RPLC axes to Lab
                axes.  3x3 np.ndarray
   _Rot2LbEnd : Calculated; otation matrix that takes RPLC axes to Lab
                axes.  3x3 np.ndarray.  Set in derived class.
    _TrnsMtrx : Calculated; transfer matrix (6x6).  Set to Null in
                __init__, initialised (to Null) in
                BeamLineElement.__init__, set in derived classes.

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

      SummaryStr: returns summary string, presently just position of
                  element.
            Input: None
           Return: str="Pos: [x, y, z] = " + str(self.getrStrt())

  Set methods:
         setDebug  : Set debug flag
       setAll2None : Set all attributes to Null
          setName  : Set name of element
         setrStrt  : Set start of element, x, y, z (m)
         setvStrt  : Set orientation of element, theta, phi (rad)
        setdrStrt  : Set offset of start of element, x, y, z (m)
        setdvStrt  : Set offset orientation of element, Euler angles (rad)
     setRot2LbStrt : set rotation matrix totransform from RLBC to lab at
                     start.

  Get methods:
         getDebug  : get debug flag
     getinstances  : get list of instances
          getName  : Get name of element
          getrStrt : Get start of element, x, y, z (m)
          getvStrt : Get orientation of element, theta, phi (rad)
         getdrStrt : Get offset from nominal start of element, x, y, z (m)
         getdvStrt : Get offset of orientation of element Euler angles (rad)
       getStrt2End : Get vector to translate from start to end in lab
     getRot2LbStrt : Get rotation matrix totransform from RLBC to lab at
                     start.
      getRot2LbEnd : Get rotation matrix totransform from RLBC to lab at end.
 getTransferMatrix : Get transfer matrix.


  Processing methods:
OutsideBeamPipe : Returns true of  particle outside beam pipe defined in
                  Facility
             Input: R: np.ndarray trace-space vector.

      Transport : Applies transfer matrix to phase-space vector.
             Input: 6D phase-space vector, np.array.
            Return: 6D phase-space vector after element

    Shift2Local : Shift (translate) to local coordinates of element
             Input: 6D phase-space vector, np.array.
            Return: Transformed 6D phase-space vector

     Shift2RPLC : Shift (translate) back to RPLC coordinates
             Input: 6D phase-space vector, np.array.
            Return: Transformed 6D phase-space vector

   
  I/o methods:
      writeElement : Class method; write element data to "dataFILE"
          i/p: dataFILE; i.o writer

       readElement : Class method; read element fromd to "dataFILE"
          i/p: dataFILE; i.o reader


  Utilities:
    cleaninstances : Class method: cleans (using "del") instances and resets
                     list.

    removeInstance : Class method: remove (using remove) instance, inst, from
                     memory and list of instances.
           input: inst; instance of BLE class


Created on Mon 12Jun23: Version history:
---------------------------------------- 
 2.0: 09Apr25: Include electron temperature update from Sadur and Zakhir.
               Also, slim down input arguments required for laser-driven
               source.
 1.1: 10Jan24: Update setTraceSpaceAtSource to make it match documented
               definititions.
 1.0: 12Jun23: First implementation

@author: kennethlong
"""

import warnings as wrnngs

from copy import deepcopy
import matplotlib.patches as patches
import scipy  as sp
import numpy  as np
import math   as mth
import random as rnd
import scipy
from scipy.optimize import fsolve
import struct as strct
import math
import pandas   as pnds
import warnings as wrngs

import BeamLine          as BL
import PhysicalConstants as PhysCnst
import Particle          as Prtcl
import LaTeX             as LTX
import BeamIO            as bmIO

#--------  Physical Constants
constants_instance = PhysCnst.PhysicalConstants()

#.. SI units:
alpha              = constants_instance.alpha()
epsilon0SI         = constants_instance.epsilon0SI()
electronCHARGESI   = constants_instance.electricCHARGE()
eps0SI               = constants_instance.epsilon0()

electronMASSSI     = constants_instance.meSI()
protonMASSSI       = constants_instance.mpSI()

speed_of_light     = constants_instance.SoL()

#.. Natural units:
electricCHARGE     = mth.sqrt(4.*mth.pi*alpha)
epsilon0           = 1.

electronMASS       = constants_instance.me()
protonMASS         = constants_instance.mp()

#.. Conversion factors and units:
Joule2MeV          = constants_instance.Joule2MeV()
m2InvMeV           = constants_instance.m2InvMeV()


class BeamLineElement:

#-------- Class attributes  --------  --------

#.. List of instances and debug flag
    instances      = []
    __Debug        = False
    __NOdispersion = False

#--------  "Built-in methods":
    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None):
        if self.getDebug():
            print(' BeamLineElement.__init__: ', \
                  'creating the BeamLineElement object')
            print("     ---->               Name:", _Name)
            print("     ---->           Position:", _rStrt)
            print("     ---->        Orientation:", _vStrt)
            print("     ---->    Position offset:", _drStrt)
            print("     ----> Orientation offset:", _dvStrt)

        BeamLineElement.instances.append(self)

        BeamLineElement.setAll2None(self)
        
        if  not isinstance( _Name, str)        or \
            not isinstance( _rStrt, np.ndarray) or \
            not isinstance( _vStrt, np.ndarray) or \
            not isinstance(_drStrt, np.ndarray) or \
            not isinstance(_dvStrt, np.ndarray):
            raise badBeamLineElement( \
                  " BeamLineElement: no default beamline element!"
                                      )
        self.setName(_Name)
        self.setrStrt(_rStrt)
        self.setvStrt(_vStrt)
        self.setdrStrt(_drStrt)
        self.setdvStrt(_dvStrt)
        self.setRot2LbStrt()
        
        if self.getDebug():
            print("     ----> New BeamLineElement instance: \n", \
                  BeamLineElement.__str__(self))
            
    def __repr__(self):
        return "BeamLineElement()"

    def __str__(self):
        print(" BeamLineElement:")
        print(" ----------------")
        print("     ---->             Debug flag:", \
              BeamLineElement.getDebug())
        if BeamLineElement.getNOdispersion():
            print("     ---->            No dispersion flag:", \
                  BeamLineElement.getNOdispersion())
        print("     ---->                   Name:", self.getName())
        print("     ---->               Position:", self.getrStrt())
        print("     ---->            Orientation: \n", self.getvStrt())
        print("     ---->        Position offset:", self.getdrStrt())
        print("     ---->     Orientation offset:", self.getdvStrt())
        print("     ---->    Magnitude of dvStrt:", \
              np.linalg.norm(self.getdvStrt()))
        if np.linalg.norm(self.getdvStrt()) != 0.:
            print("          ---->    drRotStrt:", self.getdRotStrt())
            print("          ----> drRotStrtINV:", self.getdRotStrtINV())
        print("     ---->    Start to end vector:", self.getStrt2End())
        print("     ----> Rotate to lab at start: \n", self.getRot2LbStrt())
        print("     ---->   Rotate to lab at end: \n", self.getRot2LbEnd())
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ---->        Transfer matrix: \n", \
                  self.getTransferMatrix())
        return " <---- BeamLineElement parameter dump complete."

    def SummaryStr(self):
        Str = "Pos = " + str(self.getrStrt()) + \
              " dr = " + str(self.getdrStrt()) + \
              " dv = " + str(self.getdvStrt())
        if self.getDebug():
            if np.linalg.norm(self.getdvStrt()) != 0.:
                Str += \
                    "\n dR = " + str(self.getdRotStrt()) + \
                    "\n dRINV = " + str(self.getdRotStrtINV())
        return Str

    
#--------  "Set method" only Debug
#.. Method believed to be self documenting(!)
    @classmethod
    def setDebug(self, Debug=False):
        if self.getDebug():
            print(" BeamLineElement.setdebug: ", Debug)
        self.__Debug = Debug

    @classmethod
    def setNOdispersion(cls, _NOdispersion):
        cls.__NOdispersion = _NOdispersion

    def setAll2None(self):
        self._Name       = None
        self._rStrt      = None
        self._vStrt      = None
        self._drStrt     = None
        self._dvStrt     = None
        self._Length     = None
        self._Strt2End   = None
        self._Rot2LbStrt = None
        self._Rot2LbEnd  = None
        self._TrnsMtrx   = None
    
    def setName(self, _Name):
        if not isinstance(_Name, str):
            raise badParameter(" BeamLineElement.setName: bad name:", \
                               _Name)
        self._Name = _Name
        
    def setrStrt(self, _rStrt):
        if not isinstance(_rStrt, np.ndarray):
            raise badParameter(" BeamLineElement.setrStrt: bad start:", \
                               _rStrt)
        self._rStrt = _rStrt
        
    def setvStrt(self, _vStrt):
        if not isinstance(_vStrt, np.ndarray):
            raise badParameter( \
                        " BeamLineElement.setvStrt: bad orienttion:", \
                               _vStrt)
        if len(np.shape(_vStrt)) == 1:
            raise badParameter( \
                        " BeamLineElement.setvStrt: legacy _vStrt, ", \
                               "np.shape(_vStrt)=", np.shape(_vStrt))
        self._vStrt = _vStrt
        
    def setdrStrt(self, _drStrt):
        if not isinstance(_drStrt, np.ndarray):
            raise badParameter(" BeamLineElement.setdrStrt:", \
                               " bad start offset:", \
                               _drStrt)
        self._drStrt = _drStrt
        
    def setdvStrt(self, _dvStrt):
        if not isinstance(_dvStrt, np.ndarray):
            raise badParameter(" BeamLineElement.setdvStrt:", \
                               " bad orienttion offset:", \
                               _dvStrt)

        if self.getDebug():
            print(" BeamLineElement.setdvStrt: dvStrt:", _dvStrt, \
                  _dvStrt.shape)

        if str(_dvStrt.shape) == "(2, 2)":
            if _dvStrt[0][0] == 0 and \
               _dvStrt[0][1] == 0 and \
               _dvStrt[1][0] == 0 and \
               _dvStrt[1][1] == 0:
                _dvStrt = [0., 0., 0.]
            else:
                raise need2convertdvStrt()
     
        self._dvStrt = _dvStrt

        #.. Also calculate rotation matrix:
        R1 = np.array( [ [ mth.cos(_dvStrt[0]), -mth.sin(_dvStrt[0]), 0.], \
                         [ mth.sin(_dvStrt[0]),  mth.cos(_dvStrt[0]), 0.], \
                         [                  0.,                   0., 1.] ] )
        R2 = np.array( [ [ mth.cos(_dvStrt[1]), 0., -mth.sin(_dvStrt[1])], \
                         [                  0., 1.,                   0.], \
                         [ mth.sin(_dvStrt[1]), 0.,  mth.cos(_dvStrt[1])] ] )
        R3 = np.array( [ [ mth.cos(_dvStrt[2]), -mth.sin(_dvStrt[2]), 0.], \
                         [ mth.sin(_dvStrt[2]),  mth.cos(_dvStrt[2]), 0.], \
                         [                  0.,                   0., 1.] ] )

        self._dRotStrt    = np.matmul( R3, np.matmul(R2, R1) )
        self._dRotStrtINV = np.linalg.inv(self._dRotStrt)

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> dv:", self.getdvStrt())
                print("         ----> R1:", R1)
                print("         ----> R2:", R2)
                print("         ----> R3:", R3)
                print("     ---->    dRot:", self.getdRotStrt())
                print("     ----> dRotINV:", self.getdRotStrtINV())
                print(" <---- Done.")

    def setLength(self, _Length):
        self._Length = _Length

    def setRot2LbStrt(self):
        if not isinstance(self.getvStrt(), np.ndarray):
            raise badParameter(" BeamLineElement.setRot2LbStrt:" + \
                               " vStrt not set.")
        if len(np.shape(self.getvStrt())) != 2:
            raise badParameter(" BeamLineElement.setRot2LbStrt:" + \
                               " bad vStrt=", np.shape(self.getvStrt()))

        jr = np.array( [ \
              mth.sin(self.getvStrt()[0][0]) * \
                         mth.cos(self.getvStrt()[0][1]), \
              mth.sin(self.getvStrt()[0][0]) * \
                         mth.sin(self.getvStrt()[0][1]), \
              mth.cos(self.getvStrt()[0][0]) \
                         ])
        kr = np.array( [ \
              mth.sin(self.getvStrt()[1][0]) * \
                         mth.cos(self.getvStrt()[1][1]), \
              mth.sin(self.getvStrt()[1][0]) * \
                         mth.sin(self.getvStrt()[1][1]), \
              mth.cos(self.getvStrt()[1][0]) \
                         ])
        ir = np.cross(jr, kr)

        Rot2LbStrt = np.array( [ \
                                 [ir[0], jr[0], kr[0]], \
                                 [ir[1], jr[1], kr[1]], \
                                 [ir[2], jr[2], kr[2]]  \
                                 ] )
        
        self._Rot2LbStrt = Rot2LbStrt
        
        if self.getDebug():
            print(" BeamLineElement.setRot2LbStrt:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> vStrt: \n", self.getvStrt())
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Rot2LbStrt: \n", Rot2LbStrt)
            print(" <---- BeamLineElement.setRot2LbStrt: done.")

    def setStrt2End(self, t):
        if not isinstance(t, np.ndarray):
            raise badParameter(" BeamLineElement.setStrt2End:" + \
                               " bad translation (RLPC)")

        self._Strt2End = np.matmul(self.getRot2LbStrt(), t)
        
        if self.getDebug():
            print(" BeamLineElement.setStrt2End:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> t: \n", t)
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Rot2LbStrt: \n", self.getRot2LbStrt())
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Strt2End: \n", self.getStrt2End())
            print(" <---- BeamLineElement.setStrt2End: done.")
        
    def setRot2LbEnd(self, _Rot2LbEnd):
        if not isinstance(_Rot2LbEnd, np.ndarray):
            raise badParameter(" BeamLineElement.setRot2LbEnd:", \
                               " bad argument")
        self._Rot2LbEnd = _Rot2LbEnd

        if self.getDebug():
            print(" BeamLineElement.setRot2LbEnd:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Rot2LbEnd: \n", self.getRot2LbEnd())
            print(" <---- BeamLineElement.setRot2LbEnd: done.")

            
#--------  "Get methods" only; version, reference, and constants
#.. Methods believed to be self documenting(!)

    @classmethod
    def getDebug(cls):
        return cls.__Debug

    @classmethod
    def getNOdispersion(cls):
        return cls.__NOdispersion
        
    @classmethod
    def getinstances(self):
        return self.instances

    def getName(self):
        return self._Name

    def getrStrt(self):
        return self._rStrt

    def getvStrt(self):
        return self._vStrt

    def getdrStrt(self):
        return self._drStrt

    def getdvStrt(self):
        return self._dvStrt

    def getdRotStrt(self):
        return self._dRotStrt
    
    def getdRotStrtINV(self):
        return self._dRotStrtINV
    
    def getLength(self):
        return self._Length

    def getStrt2End(self):
        return self._Strt2End

    def getRot2LbStrt(self):
        return self._Rot2LbStrt

    def getRot2LbEnd(self):
        return self._Rot2LbEnd

    def getvEnd(self):
        if self.getDebug():
            print(" BeamLineElement.getvEnd; start.")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> self.getRot2LbEnd() \n", \
                      self.getRot2LbEnd())
            
        kx = self.getRot2LbEnd()[0][2]
        ky = self.getRot2LbEnd()[1][2]
        kz = self.getRot2LbEnd()[2][2]
        
        jx = self.getRot2LbEnd()[0][1]
        jy = self.getRot2LbEnd()[1][1]
        jz = self.getRot2LbEnd()[2][1]

        ktheta = mth.acos(kz)
        sktheta = mth.sqrt(1. - kz**2)
        
        jtheta = mth.acos(jz)
        sjtheta = mth.sqrt(1. - jz**2)

        if abs(sktheta) > 1.E-12:
            kcphi   = kx/sktheta
            ksphi   = ky/sktheta
        else:
            kcphi = 1.
            ksphi = 0.

        if self.getDebug():
            print("     ----> jx, jy, sjtheta:", jx, jy, sjtheta)
        jcphi   = jx/sjtheta
        jsphi   = jy/sjtheta

        kphi = mth.atan2(ksphi, kcphi)
        if kphi < 0.: kphi = 2.*mth.pi + kphi

        jphi = mth.atan2(jsphi, jcphi)
        if jphi < 0.: jphi = 2.*mth.pi + jphi

        vEnd = np.array( [ [jtheta, jphi], [ktheta, kphi] ])

        return vEnd

    def getMode(self):
        return self._Mode
    
    def getModeText(self):
        return self._ModeText
    
    def getParameterText(self):
        return self._ParameterText
    
    def getParameterUnit(self):
        return self._ParameterUnit
    
    def getParameters(self):
        return self._Params
    
    def getTransferMatrix(self):
        return self._TrnsMtrx

    def getLines(self):
        Lines = []
        return Lines


#--------  Processing methods:
    def OutsideBeamPipe(self, _R):
        Outside = False
        Rad = np.sqrt(_R[0]**2 + _R[2]**2)
        if Rad >= Facility.getinstances().getVCMVr():
            Outside = True
        return Outside

    def ExpansionParameterFail(self, _R):
        iLctn = BeamLineElement.getinstances().index(self)
        iAddr = iLctn - 1
        
        if self.getDebug():
            print(" Particle.ExpansionParameterFail: start")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> TraceSpace:", _R)
                print("     ----> iLctn:", iLctn, self.getName())
                print("     ----> iAddr:", iAddr)
                
        Fail = False
        
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()

        p0    = mth.sqrt(np.dot(iRefPrtcl.getPrIn()[iAddr][:3], \
                                iRefPrtcl.getPrIn()[iAddr][:3]))
        E0    = iRefPrtcl.getPrOut()[iAddr][3]
        b0    = p0/E0
        D     = mth.sqrt(1. + \
                         2.*_R[5]/b0 +
                         _R[5]**2)
        eps   = ( _R[1]**2 + _R[3]**2  ) / (2.*D**2)
        if self.getDebug():
            print("     ----> p0, E0, b0, D:", p0, E0, b0, D)
            print("     <---- Epsilon:", eps)
            
        if eps > 1.0:
            Fail = True

        if self.getDebug():
            print(" <---- Return, Fail:", Fail)

        return Fail

    def Transport(self, __R):        #<---- class BeamLineElement:
        #.. Protect input vector; planning to transform it to coordinate
        #   system referred to beam-line element
        _R = deepcopy(__R)
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
                        " BeamLineElement.Transport: bad input vector:", \
                                _R)

        if self.getDebug():
            print(" BeamLineElement.Transport:", \
                  Facility.getinstances().getName(), \
                  Facility.getinstances().getVCMVr())
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> _R:", _R)

        _R = self.Shift2Local(_R)
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Shift to element-centred coordinates:")
                print("           _R:", _R)
                
        _R = self.Tilt2Local(_R)
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Tilt to element-centred coordinates:")
                print("           _R:", _R)
                
        if self.getDebug():
            print("     ----> Outside:", self.OutsideBeamPipe(_R))
            print("     ----> Expansion parameter fail:", \
                  self.ExpansionParameterFail(_R))

        NotCut = True
        if isinstance(self, Aperture):
            NotCut = self.passTHROUGH(_R)
            if not NotCut: return None
        if self.OutsideBeamPipe(_R) or \
           self.ExpansionParameterFail(_R) or \
           abs(_R[4]) > 5.:
            return None
        else:
            if isinstance(self, DefocusQuadrupole) or \
               isinstance(self, FocusQuadrupole)   or \
               isinstance(self, Solenoid)          or \
               isinstance(self, SectorDipole)      or \
               isinstance(self, GaborLens)         or \
               isinstance(self, QuadDoublet)       or \
               isinstance(self, QuadTriplet):         \
                self.setTransferMatrix(_R)

            detTrnsfrMtrx = np.linalg.det(self.getTransferMatrix())
            error         = abs(1. - abs(detTrnsfrMtrx))
            if error > 1.E-6:
                print(" BeamLineElement.Transport: detTrnsfrMtrx:", \
                      detTrnsfrMtrx)

            _Rprime = self.getTransferMatrix().dot(_R)

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Rprime:", _Rprime)

        if isinstance(_Rprime, np.ndarray):
            _Rprime = self.Tilt2RPLC(_Rprime)
            _Rprime = self.Shift2RPLC(_Rprime)
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print( \
                "     ----> Shift back from element-centred coordinates:")
                print("           _Rprime:", _Rprime)
                
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" <---- Rprime:", _Rprime)

        if not isinstance(_Rprime, np.ndarray):
            pass

        return _Rprime

    def Shift2Local(self, _R):
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
                        " BeamLineElement.Shift2Local: bad input vector:", \
                                _R)
        
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" Shift2Local: R:", _R)

        _Rprime    = deepcopy(_R)
        _Rprime[0] -= self._drStrt[0]
        _Rprime[2] -= self._drStrt[1]

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" <---- Rprime:", _Rprime)
        
        return _Rprime

    def Tilt2Local(self, _R):
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
                        " BeamLineElement.Tilt2Local: bad input vector:", \
                                _R)
        
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" Tilt2Local: R:", _R)
                print("     ----> Norm(dv):", np.linalg.norm(self.getdvStrt()))

        _Rprime = deepcopy(_R)

        if np.linalg.norm(self.getdvStrt()) != 0:
            r    = np.array([_Rprime[0], _Rprime[2], _Rprime[4]])
            rPRM = np.matmul(self.getdRotStrtINV(), r)
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=7,suppress=True):
                    print("     ---->       r:", r)
                    print("     ----> dRotINV: \n", self.getdRotStrtINV())
                    print("     <----    rPRM:", rPRM)
                    
            dzds   = mth.sqrt(1. - _Rprime[1]**2 - _Rprime[3]**2)
            vec    = np.array([_Rprime[1], _Rprime[3], dzds])
            vecPRM = np.matmul(self.getdRotStrtINV(), vec)
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=7,suppress=True):
                    print("     ---->    vec:", vec)
                    print("     ----> dRotINV: \n", self.getdRotStrtINV())
                    print("     <---- vecPRM:", vecPRM)
                
            _Rprime[0] = rPRM[0]
            _Rprime[1] = vecPRM[0]
            _Rprime[2] = rPRM[1]
            _Rprime[3] = vecPRM[1]
            _Rprime[4] = rPRM[2]

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" <---- Rprime:", _Rprime)
        
        return _Rprime

    def Shift2RPLC(self, _R):
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
                        " BeamLineElement.Shift2RPLC: bad input vector:", \
                                _R)
        
        _Rprime    = deepcopy(_R)
        _Rprime[0] += self._drStrt[0]
        _Rprime[2] += self._drStrt[1]
        
        return _Rprime

    def Tilt2RPLC(self, _R):
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
                        " BeamLineElement.Tilt2RPLC: bad input vector:", \
                                _R)
        
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" Tilt2RPLC: R:", _R)
                print("     ----> Norm(dv):", np.linalg.norm(self.getdvStrt()))

        _Rprime = deepcopy(_R)

        if np.linalg.norm(self.getdvStrt()) != 0:
            r    = np.array([_Rprime[0], _Rprime[2], _Rprime[4]])
            rPRM = np.matmul(self.getdRotStrt(), r)
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=7,suppress=True):
                    print("     ---->       r:", r)
                    print("     ----> dRot: \n", self.getdRotStrt())
                    print("     <----    rPRM:", rPRM)
                    
            dzds   = mth.sqrt(1. - _Rprime[1]**2 - _Rprime[3]**2)
            vec    = np.array([_Rprime[1], _Rprime[3], dzds])
            vecPRM = np.matmul(self.getdRotStrt(), vec)
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=7,suppress=True):
                    print("     ---->    vec:", vec)
                    print("     ---->   dRot: \n", self.getdRotStrt())
                    print("     <---- vecPRM:", vecPRM)
                
            _Rprime[0] = rPRM[0]
            _Rprime[1] = vecPRM[0]
            _Rprime[2] = rPRM[1]
            _Rprime[3] = vecPRM[1]
            _Rprime[4] = rPRM[2]

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" <---- Rprime:", _Rprime)
        
        return _Rprime

    
#--------  I/o methods:
    def writeElement(self, dataFILE):        #<---- class BeamLineElement:
        if self.getDebug():
            print(" BeamLineElement.writeElement starts.")

        bLocation = bytes(self.getName(), 'utf-8')
        record    = strct.pack(">i", len(bLocation))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Length of element name:", \
                  strct.unpack(">i", record))
        
        record = bLocation
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Location:", bLocation.decode('utf-8'))

        record = strct.pack(">13d",                 \
                            self.getrStrt()[0],     \
                            self.getrStrt()[1],     \
                            self.getrStrt()[2],     \
                            self.getvStrt()[0][0],  \
                            self.getvStrt()[0][1],  \
                            self.getvStrt()[1][0],  \
                            self.getvStrt()[1][1],  \
                            self.getdrStrt()[0],    \
                            self.getdrStrt()[1],    \
                            self.getdrStrt()[2],    \
                            self.getdvStrt()[0], \
                            self.getdvStrt()[1], \
                            self.getdvStrt()[2]
                            )
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> rStrt, vStrt:", \
                      strct.unpack(">13d",record))
        
        if self.getDebug():
            print(" <---- BeamLineElement.writeElement done.")
        
    @classmethod
    def readElement(cls, dataFILEinst):        #<---- class BeamLineElement:
        if cls.getDebug():
            print(" BeamLineElement.readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        nFls = 0
        for ibmIOr in filter(lambda iFL : \
                    iFL.getdataFILE() == dataFILE,\
                    bmIO.BeamIO.getinstances()):
            nFls += 1

        if nFls != 1:
            raise cantFINDfile()

        brecord = dataFILE.read(4)
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True
            
        record = strct.unpack(">i", brecord)
        nChr   = record[0]
        if cls.getDebug():
            print("     ----> Number of characters:", nChr)
                
        brecord  = dataFILE.read(nChr)
        Location = brecord.decode('utf-8')
        if cls.getDebug():
            print("                   Location:", Location)

        if ibmIOr.getdataFILEversion() < 4:
            brecord = dataFILE.read((7*8))
        elif ibmIOr.getdataFILEversion() == 4:
            brecord = dataFILE.read((14*8))
        else:
            brecord = dataFILE.read((13*8))
        if brecord == b'':
            return True, None, None, None, None
        
        if ibmIOr.getdataFILEversion() < 4:
            record = strct.unpack(">7d", brecord)
        elif ibmIOr.getdataFILEversion() == 4:
            record = strct.unpack(">14d", brecord)
        else:
            record = strct.unpack(">13d", brecord)
        r      = np.array([float(record[0]), float(record[1]),       \
                           float(record[2])])
        v      = np.array([[float(record[3]), float(record[4])],     \
                           [float(record[5]), float(record[6])]]  )
        
        if ibmIOr.getdataFILEversion() < 4:
            dr     = np.array([0., 0., 0.])
            dv     = np.array([0., 0., 0.])
        elif ibmIOr.getdataFILEversion() == 4:
            dr     = np.array([float(record[7]), float(record[8]),       \
                               float(record[9])])
            dv     = np.array([[0., 0.], [0., 0.]])
        else:
            dr     = np.array([float(record[7]), float(record[8]),       \
                               float(record[9])])
            dv     = np.array([ float(record[10]), float(record[11]),  \
                                float(record[12]) ])
        
        if cls.getDebug():
            print("     ----> r, v, dr, dv:", r, v, dr, dv)

        return EoF, Location, r, v, dr, dv

    
#--------  Utilities:
    @classmethod
    def cleaninstances(cls):
        if cls.getDebug():
            print(' BeamLineElement.cleaninstance: clean instances.')
        for inst in cls.instances:
            if cls.getDebug():
                print("     ----> Kill:", inst.getName())
            if isinstance(inst, Facility):
                Facility.instance = None
            del inst
        cls.instances = []
        if cls.getDebug():
            print(' <---- BeamLineElement.cleaninstance: instances removed.')

    @classmethod
    def removeInstance(cls, inst):
        if inst in cls.getinstances():
            if cls.getDebug():
                print(" BeamLineElement.removeInstance: remove", \
                      inst.getName(), "from beamline instances list")
            cls.getinstances().remove(inst)
        else:
            if cls.getDebug():
                print(" BeamLineElement.removeInstance: instance", \
                      inst.Name, "not in BeamLineElement.Instances!")
                
    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" BeamLineElement.visualise: start")
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        xlim  = axs.get_xlim()
        ylim  = axs.get_ylim()

        wdth = self.getLength()
        if wdth == 0.: wdth = (xlim[1] - xlim[0]) / 100.
        hght = (ylim[1] - ylim[0]) / 10.
        angl = 0.
        abt  = 'center'
            
        if CoordSys == "RPLC":
            if self.getDebug():
                print("     ----> RPLC:")
                
            iRefPrtcl = Prtcl.ReferenceParticle.getinstances("All")
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)

            sxy   = [ sStrt, -hght/2. ]
                
        elif CoordSys == "Lab":
            if self.getDebug():
                print("     ----> Lab:")
                
            strt  = self.getrStrt()
            if self.getDebug():
                print("         ----> Centre:", strt)

            if Proj == "xz":
                sxy   = [ strt[2], strt[0]-hght/2. ]
            elif Proj == "yz":
                sxy   = [ strt[2], strt[1]-hght/2. ]
        

        if self.getDebug():
            print("     ---->    sxy:", sxy)
            print("     ---->   wdth:", wdth)
            print("     ---->   hght:", hght)
            print("     ---->   angl:", angl)
            print("     ---->    abt:", abt)
            
        Patch = patches.Rectangle(sxy, wdth, hght, \
                                  angle=angl, \
                                  rotation_point=abt, \
                                  facecolor=('green', 0.5), \
                                  zorder=2)

        axs.add_patch(Patch)
                                   
        if self.getDebug():
            print(" <---- BeamLineElement.visualise: ends.")

#--------  Derived classes  --------  --------  --------  --------  --------
"""
Derived class Facility:
=======================

  Facility class derived from BeamLineElement to contain paramters for a
  Facility space.  __init__ sets parameters of Facility.  Methods to
  apply transfer matix.


  Class attributes:
  -----------------
    instance : List of instances of Facility(BeamLineElement) class
  __Debug    : Debug flag


  Parent class attributes:
  ------------------------
        _Name : str : Name of facility
       _rStrt : numpy array; x, y, z position (in m) of start of
                element. 
       _vStrt : numpy array; theta, phi of principal axis of element.
      _drStrt : "error", displacement of start from nominal position
                (rad). 
      _dvStrt : "error", deviation in theta and phy from nominal axis
                (rad).

  Calculated and set in derived classes:
    _Strt2End, _Rot2LbStrt, _Rot2LbEnd, _TrnsMtrx not needed, so not
    set.  Will remain "None" as set in parent class.


  Instance attributes to define Facility:
  ---------------------------------------
  _p0   : float : Reference particle momentum; unit MeV
 _VCMVr : float : Vacuum chamber mother volume; unit m.  Idea here is to
                  define the radius at which a particle trajectory is
                  terminated.  It may be necessary to introduce a beam
                  pipe later.
 _species0in: str of list of strings :
                  Gives list of species of reference particles to be
                  propagated through facility.
    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Set methods:
          setp0 : float : e.g. 15 MeV
       setVCMVr : float : e.g. 0.5 m
    setspecies0 : list of strings : e.g. [pion, muon]

  Get methods:
     getp0, getVCMVr, getspecies0

  I/o methods:
      writeElement : Class method; write element data to "dataFILE"
          i/p: dataFILE; i.o writer

       readElement : Class method; read element fromd to "dataFILE"
          i/p: dataFILE; i.o reader


"""
class Facility(BeamLineElement):
    instance  = None
    __Debug   = False

    
#--------  "Built-in methods":
    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _p0in=None, _VCMVr=None, _species0in=None):

        if self.__Debug:
            print(' Facility.__init__: ', \
                      'entering the Facility object creation method:')

        if Facility.instance == None:
            Facility.instance = self
            
            if self.__Debug:
                print('     ----> Creating the Facility object:')
                print("         ----> Name:", _Name)
                print("         ----> Reference particle momentum:", _p0in)
                print("         ----> Vacuum chamber mother volume radius", \
                      _VCMVr)
                print("         ----> Species of reference particle(s):", \
                      _species0in)

            Facility.instance = self

            #.. BeamLineElement class initialisation:
            BeamLineElement.__init__(self, \
                                     _Name, _rStrt, _vStrt, _drStrt, _dvStrt)

            if _p0in == None:
                raise badBeamLineElement( \
            " Facility: bad specification for reference particle mometnum!"
                                         )
            if isinstance(_p0in, float):
                _p0 = [_p0in]
            elif isinstance(_p0in, list):
                _p0 = []
                for p0 in _p0in:
                    if isinstance(p0, float):
                        _p0.append(p0)
                    else:
                        raise badBeamLineElement( \
            " Facility: bad specification for reference particle mometnum!"
                                                 )
            
            if not isinstance(_VCMVr, float):
                raise badBeamLineElement( \
            " Facility: bad specification for vacuum chamber mother volume!"
                                         )

            if _species0in == None or _species0in == []:
                wrngs.warn("Facility(BeamLineElement): no species given," + \
                           " use default, proton")
                _species0 = ["proton"]
            elif isinstance(_species0in, str):
                _species0 = [_species0in]
            elif isinstance(_species0in, list):
                _species0 = _species0in
            if not self.validspecies0(_species0):
                raise badBeamLineElement( \
                        " Facility: bad species0:" + str(_species0in) \
                                         )

            if self.getDebug():
                print("         ---->       _p0:", _p0, "\n", \
                      "              _species0:", _species0)

            if len(_p0) != len(_species0):
                raise badBeamLineElement( \
            " Facility: length of p0 and species0 not the same!" \
                                         )
            self.setp0(_p0)
            self.setVCMVr(_VCMVr)
            self.setspecies0(_species0)
                
            self.setLength(0.)
            self.setStrt2End(np.array([0., 0., self.getLength()]))
            self.setRot2LbEnd(self.getRot2LbStrt())
        
            if self.__Debug:
                print("     ----> New Facility instance: \n", \
                      self)
        else:
            print(' Facility(BeamLineElement).__init__: ',       \
                  " attempt to create facility.", \
                  " Abort!")
            raise secondFacility(" Second call not allowed.")
            
    def __repr__(self):
        return "Facility()"

    def __str__(self):
        print(" Facility:")
        print(" ------")
        print("     ----> Debug flag:", Facility.getDebug())
        print("     ----> Name      :", self.getName())
        print("     ----> p0 (MeV/c):", self.getp0())
        print("     ----> Vacuum chamber mother volume radius (m):", \
              self.getVCMVr())
        print("     ----> List of species:", self.getspecies0())
        BeamLineElement.__str__(self)
        return " <---- Facility parameter dump complete."

    def SummaryStr(self):
        Str  = "Facility         : " + BeamLineElement.SummaryStr(self) + \
            "; Name = " + self.getName() + "; p0 = " + str(self.getp0()) + \
            "; VCMVr = " + str(self.getVCMVr()) + \
            "; species0 = " + str(self.getspecies0()) 
        return Str

    def validspecies0(self, spcs0):
        if self.getDebug():
            print(" Facility.validspecies0: start:")
            print("     ---->   Stable species:", \
                  Prtcl.Particle.stable_species)
            print("     ----> Unstable species:", \
                  Prtcl.Particle.unstable_species)
            
        success = True
        for spcs in spcs0:
            if self.getDebug():
                print("     ----> Test:", spcs)
            if not(spcs.lower() in Prtcl.Particle.stable_species) and \
               not(spcs.lower() in Prtcl.Particle.unstable_species):
                success = False
        if self.getDebug():
            print(" <---- Success:", success)
        
        return success

#--------  "Set methods"
#.. Methods believed to be self documenting(!)
    def setp0(self, _p0=None):
        if not isinstance(_p0, list):
            raise badParameter( \
                     " BeamLineElement.Facility.setp0: bad p0",
                                _p0)
        self._p0 = _p0

    def setVCMVr(self, _VCMVr=None):
        if not isinstance(_VCMVr, float):
            raise badParameter( \
                     " BeamLineElement.Facility.setVCMVr: bad VCMVr",
                                _VCMVr)
        self._VCMVr = _VCMVr

    def setspecies0(self, _species0=None):
        if not isinstance(_species0, list):
            raise badParameter( \
                     " BeamLineElement.Facility.setspecies0: bad species0",
                                _species0)
        self._species0 = _species0

        
#--------  "get methods"
#.. Methods believed to be self documenting(!)
    @classmethod
    def getinstances(cls):
        return cls.instance
    
    def getp0(self):
        return self._p0
    
    def getVCMVr(self):
        return self._VCMVr
    
    def getspecies0(self):
        return self._species0
    

#--------  I/o methods:
    def getLines(self):
        Lines = []

        print(" TimeBomb!  This routine was modified to deal with", \
              "list of reference particles; need to check it!")
        
        Stage   = 0
        Section = "Facility"
        Element = "Global"
        Type    = "Name"
        Param   = "Name"
        Value   = self.getName()
        Unit    = ""
        Comment = ""
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        Type    = "Reference particle"
        Param   = "Kinetic energy"
        Unit    = "MeV"
        idRfPrtcl = 0
        for iRefPrtcl in Prtcl.ReferenceParticle.getinstances("All"):
            p0      = self.getp0()[idRfPrtcl]
            particleMASS = PhysCnst.PhysicalConstants().getparticleMASS( \
                                        self.getspecies0()[idRfPrtcl] )
            E       = mth.sqrt( particleMASS**2 + p0**2)
            K       = E - particleMASS
            Value   = K
            Lines.append([Stage, Section, Element, Type, \
                          Param, Value, Unit, Comment])
            
            idRfPrtcl += 1
            Param = Param + str(idRfPrtcl)

        Type    = "Reference particle"
        Param   = "Species"
        Unit    = ""
        idRfPrtcl = 0
        for iRefPrtcl in Prtcl.ReferenceParticle.getinstances("All"):
            Value = self.getspecies0()[idRfPrtcl]

            Value   = self.getspecies0()[idRfPrtcl]
            Lines.append([Stage, Section, Element, Type, \
                          Param, Value, Unit, Comment])
            
            idRfPrtcl += 1
            Param = Param + str(idRfPrtcl)
            
        Type    = "Vacuum chamber"
        Param   = "Mother volume radius"
        VCMVr   = self.getVCMVr()
        Value   = VCMVr
        Unit    = "m"
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        return Lines

    def writeElement(self, dataFILE):
        if self.getDebug():
            print(" Facility(BeamLineElement).writeElement starts.")

        derivedCLASS = "Facility"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Derived class:", bversion.decode('utf-8'))
            
        #.. Write reference particle 0 and VCMV:
        if self.getDebug():
            print("     ----> Reference particle 0 and VCMV:")
        record = strct.pack(">2d", self.getp0()[0], self.getVCMVr())
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> p0, VCMVr:", strct.unpack(">2d",record))

        species = self.getspecies0()[0]
        bversion = bytes(species, 'utf-8')
        record   = strct.pack(">i", len(species))
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length of species 0 record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Species:", bversion.decode('utf-8'))
        
        #.. Write reference particles for decay products:
        if self.getDebug():
            print("     ----> Reference particle(s) 1:",
                  len(self.getspecies0())-1)

        record     = strct.pack(">i", len(self.getspecies0())-1 )
        dataFILE.write(record)
        if self.getDebug():
            print( \
            "         ----> Number of decay-produce reference particles:", \
                  strct.unpack(">i", record))
        
        for idREFprtcl in range(1,len(self.getspecies0())):
            if idREFprtcl ==0: continue
            
            if self.getDebug():
                print("         ----> idREFprtcl, species, p0:", \
                      idREFprtcl, \
                      self.getspecies0()[idREFprtcl], \
                      self.getp0()[idREFprtcl])

            record = strct.pack(">1d", self.getp0()[idREFprtcl])
            dataFILE.write(record)
            if self.getDebug():
                print("             ----> idREFprtcl, p0:", idREFprtcl,\
                      strct.unpack(">1d",record))
                
            species = self.getspecies0()[idREFprtcl]
            bversion = bytes(species, 'utf-8')
            record   = strct.pack(">i", len(species))
            dataFILE.write(record)
            if self.getDebug():
                print("             ----> Length of species record:", \
                      strct.unpack(">i", record))
            record   = bversion
            dataFILE.write(record)
            if self.getDebug():
                print("             ----> Species:", bversion.decode('utf-8'))

        if self.getDebug():
            print("     <---- Done with Facility specific stuff;", \
                  "Facility(BeamLineElement).writeElement starts.")
            
        BeamLineElement.writeElement(self, dataFILE)

        if self.getDebug():
            print(" <---- Facility(BeamLineElement).writeElement done.")

    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" Facility(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        p0list       = []
        species0list = []
        
        EoF      = False
        species0 = "proton"

        brecord = dataFILE.read((2*8))
        if brecord == b'':
            return True, None, None, None
        
        record  = strct.unpack(">2d", brecord)
        p0      = float(record[0])
        VCMVr   = float(record[1])
        
        p0list.append(p0)

        if cls.getDebug():
            print("     ----> Data file version:", \
                  dataFILEinst.getdataFILEversion())
        if dataFILEinst.getdataFILEversion() < 7:
            if cls.getDebug():
                print("     ----> p0, VCMVr, species0:", \
                      p0list, VCMVr, species0list)
            return EoF, p0list, VCMVr, species0list
        
        brecord = dataFILE.read(4)
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True, p0list, VCMVr, species0list

        record = strct.unpack(">i", brecord)
        nChr   = record[0]
        if cls.getDebug():
            print("     ----> Number of characters:", nChr)
        
        brecord  = dataFILE.read(nChr)
        species0 = brecord.decode('utf-8')

        species0list.append(species0)

        if dataFILEinst.getdataFILEversion() < 8:
            if cls.getDebug():
                print(" <---- p0, VCMVr:", p0list, VCMVr, species0list)

            return EoF, p0list, VCMVr, species0list

        #.. Reference particles for decay products from version 8:
        brecord = dataFILE.read(4)
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True, p0list, VCMVr, species0list
        record = strct.unpack(">i", brecord)
        nDCYprdctREFprtcls = record[0]
        if cls.getDebug():
            print("     ----> Number of decay produce reference particles:", \
                  nDCYprdctREFprtcls)

        for iRefPrtcl in range(nDCYprdctREFprtcls):
            brecord = dataFILE.read((1*8))
            if brecord == b'':
                return True, p0list, VCMVr, species0list
            record  = strct.unpack(">1d", brecord)
            p0      = float(record[0])
            if cls.getDebug():
                print("     ----> p0:", p0)

            p0list.append(p0)    

            brecord = dataFILE.read(4)
            if brecord == b'':
                if cls.getDebug():
                    print(" <---- end of file, return.")
                    return True, p0list, VCMVr, species0list
            record = strct.unpack(">i", brecord)
            nChr   = record[0]
            if cls.getDebug():
                print("     ----> Number of characters:", nChr)
            brecord = dataFILE.read(nChr)
            species = brecord.decode('utf-8')
            if cls.getDebug():
                print("     ----> Species:", species)

            species0list.append(species)
            
        if cls.getDebug():
            print(" <---- p0list, species0list:", p0list, species0list)

        return EoF, p0list, VCMVr, species0list
        
"""
Derived class Drift:
====================

  Drift class derived from BeamLineElement to contain paramters for a drift
  space.  __init__ sets parameters of drift.  Methods to apply transfer
  matix.


  Class attributes:
  -----------------
    instances : List of instances of Drift(BeamLineElement) class
  __Debug     : Debug flag


  Parent class attributes:
  ------------------------
        _Name : str : Name of facility
       _rStrt : numpy array; x, y, z position (in m) of start of
                element. 
       _vStrt : numpy array; theta, phi of principal axis of element.
      _drStrt : "error", displacement of start from nominal position
                (rad). 
      _dvStrt : "error", deviation in theta and phy from nominal axis
                (rad).

  Calculated and set in derived classes:
    _Strt2End, _Rot2LbStrt, _Rot2LbEnd, _TrnsMtrx not needed, so not
    set.  


  Instance attributes to define drift:
  ------------------------------------
  _Length : Length (in m) of drift length
  
    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Set methods:
        setLength  : Set _Length
 setTransferMatrix : "Calculate" amd set transfer matrix.
                Input: _Length (m)

  Get methods:
     getLength  : get debug flag

"""
class Drift(BeamLineElement):
    instances  = []
    __Debug    = False

    
#--------  "Built-in methods":
    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Length=None):
        if self.__Debug:
            print(' Drift.__init__: ', \
                  'creating the Drift object: Length=', _Length)

        Drift.instances.append(self)

        #.. BeamLineElement class initialisation:
        BeamLineElement.__init__(self, _Name, _rStrt, _vStrt, _drStrt, _dvStrt)

        if not isinstance(_Length, float):
            raise badBeamLineElement( \
                  " Drift: bad specification for length of drift!"
                                      )
        self.setLength(_Length)

        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        self.setTransferMatrix()
                
        if self.__Debug:
            print("     ----> New drift instance: \n", \
                  self)
            
    def __repr__(self):
        return "Drift()"

    def __str__(self):
        print(" Drift:")
        print(" ------")
        print("     ----> Debug flag:", Drift.getDebug())
        print("     ----> Length (m):", self.getLength())
        print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- Drift parameter dump complete."

    def SummaryStr(self):
        Str  = "Drift            : " + BeamLineElement.SummaryStr(self) + \
            "; length = " + str(self.getLength())
        return Str


#--------  "Set methods"
#.. Methods believed to be self documenting(!)

    def setLength(self, _Length):
        if not isinstance(_Length, float):
            raise badParameter(" BeamLineElement.Drift.setLength: bad length:",
                               _Length)
        self._Length = _Length

    def setTransferMatrix(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        E0        = iRefPrtcl.getPrOut()[iPrev][3]
        b02       = (p0/E0)**2
        g02       = 1./(1.-b02)
        
        if self.getDebug():
            print(" Drift(BeamLineElement).setTransferMatrix:")
            print("     ----> Reference particle 4-mmtm:", \
                  iRefPrtcl.getPrIn()[0])
            print("         ----> p0, E0:", p0, E0)
            print("     <---- b02, g02:", b02, g02)

        l = self.getLength()

        TrnsMtrx = np.array( [ \
                              [1.,  l, 0., 0., 0.,        0.], \
                              [0., 1., 0., 0., 0.,        0.], \
                              [0., 0., 1.,  l, 0.,        0.], \
                              [0., 0., 0., 1., 0.,        0.], \
                              [0., 0., 0., 0., 1., l/b02/g02], \
                              [0., 0., 0., 0., 0.,        1.]  \
                             ] )

        if self.getDebug():
            print("     ----> Drift transfer matrix:")
            print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx
        
        if self.getDebug():
            print(" <---- Done.")
        
#--------  "get methods"
#.. Methods believed to be self documenting(!)

    def getLength(self):
        return self._Length
    
        
#--------  I/o methods:
    def getLines(self):
        Lines = []

        Fields  = self.getName().split(":")
        
        Stage   = Fields[1]
        Section = Fields[2]
        Element = Fields[3]
        Type    = ""
        
        Param   = "Length"
        Value   = self.getLength()
        Unit    = "m"
        Comment = ""
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        return Lines
    
    def writeElement(self, dataFILE):
        if self.getDebug():
            print(" Drift(BeamLineElement).writeElement starts.")

        derivedCLASS = "Drift"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Derived class:", bversion.decode('utf-8'))

        if self.getDebug():
            print("     ----> Write parameter:")
        record = strct.pack(">d", self.getLength())
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length:", strct.unpack(">d", record))
        if self.getDebug():
            print("     <---- Done.")
            
        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print(" <---- Drift(BeamLineElement).writeElement done.")
    
    
    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" Drift(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read((1*8))
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True
            
        record  = strct.unpack(">d", brecord)
        Length  = float(record[0])
        if cls.getDebug():
            print("     ----> Length:", Length)
            
        return EoF, Length


"""
Derived class Aperture:
=======================

  Aperture class derived from BeamLineElement to contain paramters for an
  aperture.  __init__ sets parameters of the aperture.  Methods to apply
  transfer cuts in position and the transfer matix are provided.


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
  Calling arguments:
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.
  _Param : Array, Type, then paramters for aperture.

  _TrnsMtrx : Transfer matrix.


  Instance attributes to define aperture:
  ---------------------------------------
  _Type      : Type of aperture:
               = 0: circular
               = 1: Elliptical
               = 2: Rectangular
 _Params     : Circular:
               _Params[0]: Radius of circular aperture
               Elliptical:
               _Params[0]: Radius of ellipse in horizontal (x) direction
               _Params[1]: Radius of ellipse in horizontal (y) direction
               Rectangular:
               _Params[0]: Radius of ellipse in horizontal (x) direction
               _Params[1]: Radius of ellipse in horizontal (y) direction

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising aperture
                parameterrs.

  Set methods:
setAoertureParameters: Sets aperture parameters:
                 Input: _Param: same i/p arguments as Aperture __init__
 setTransferMatrix : "Calculate" amd set transfer matrix.
                Input: _Length (m)

  Get methods:
     getDebug: get debug flag, bool
      getType: Get type, return Type, int
    getParameters: Return Params, list
   getLenbgth: Returns length of aperture (presently 0)

  Utilities:


"""
class Aperture(BeamLineElement):
    instances  = []
    __Debug    = False

#--------  "Built-in methods":    
    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Param=[]):
        if self.getDebug():
            print(' Aperture.__init__: ', \
                  'creating the Aperture object')
            print("     ----> Parameters:", _Param)

        Aperture.instances.append(self)
        
        #.. BeamLineElement class initialisation:
        BeamLineElement.__init__(self, _Name, \
                                 _rStrt, _vStrt, _drStrt, _dvStrt)

        if len(_Param) < 2:
            raise badBeamLineElement( \
                  " Aperture: bad specification for aperture!"
                                      )
        self.setApertureParameters(_Param)
        
        self.setStrt2End(np.array([0., 0., 0.]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        self.setTransferMatrix()
                
        if self.getDebug():
            print("     ----> New Aperture instance: \n", \
                  self)
            
    def __repr__(self):
        return "Aperture()"
    
    def __str__(self):
        print(" Aperture:")
        print(" ---------")
        print("     ----> Debug flag:", Aperture.getDebug())
        print("     ----> Type:", self.getType())
        if self.getType() == 0:
            print("     ----> Circular:")
            print("     ----> Radius (m)", self.getParameters()[0])
        elif self.getType() == 1:
            print("     ----> Elliptical:")
            print("     ----> Radius x, y (m)", self.getParameters())
        elif self.getType() == 2:
            print("     ----> Rectangular:")
            print("     ----> Half length x, y (m)", self.getParameters())
        BeamLineElement.__str__(self)
        return " <---- Aperture parameter dump complete."

    def SummaryStr(self):
        Str  = "Aperture         : " + BeamLineElement.SummaryStr(self) + \
            "; Type = " + str(self.getType()) + \
            "; Parameters = " + str(self.getParameters())
        return Str

    
#--------  "Set methods".
#.. Methods believed to be self documenting(!)
    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug
        
    def setApertureParameters(self, _Param):
        if self.getDebug():
            print(" Apperture.setApertureParamters; Parameters:", _Param)
        if not isinstance(_Param[0], int):
            raise badParameter( \
                        " BeamLineElement.Aperture.setApertureParameters:",\
                        " bad type:",
                        _Param[0])
        self._Type = _Param[0]

        self._Params = np.array([])
        if _Param[0] == 0:           #.. Circular aperture
            if not isinstance(_Param[1], float):
                raise badParameter( \
                        " BeamLineElement.Aperture.setApertureParameters:",\
                        " bad radius:",
                        _Param[1])
            self._Params = np.append(self._Params, _Param[1])
        elif _Param[0] == 1:           #.. Elliptical aperture
            if not isinstance(_Param[1], float) or \
               not isinstance(_Param[2], float):
                raise badParameter( \
                        " BeamLineElement.Aperture.setApertureParameters:",\
                        " bad radius:",
                        _Param[1], _Param[1])
            self._Params = np.append(self._Params, _Param[1])
            self._Params = np.append(self._Params, _Param[2])
        elif _Param[0] == 2:           #.. Rectanglar aperture
            if not isinstance(_Param[1], float) or \
               not isinstance(_Param[2], float):
                raise badParameter( \
                        " BeamLineElement.Aperture.setApertureParameters:",\
                        " half length:",
                        _Param[1], _Param[1])
            self._Params = np.append(self._Params, _Param[1])
            self._Params = np.append(self._Params, _Param[2])

    def setTransferMatrix(self):
        
        TrnsMtrx = np.array( [ \
                              [1., 0., 0., 0., 0., 0.], \
                              [0., 1., 0., 0., 0., 0.], \
                              [0., 0., 1., 0., 0., 0.], \
                              [0., 0., 0., 1., 0., 0.], \
                              [0., 0., 0., 0., 1., 0.], \
                              [0., 0., 0., 0., 0., 1.]  \
                             ] )
        self._TrnsMtrx = TrnsMtrx


#--------  "get methods"
#.. Methods believed to be self documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug
    
    def getType(self):
        return self._Type
    
    def getLength(self):
        return 0.

    
#--------  Processing method:
    def passTHROUGH(self, _R):
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
                        " BeamLineElement.passTHROUGH: bad input vector:", \
                                _R)

        if self.getDebug():
            print(" Aperture(BeamLineElement).passTHROUGH:", \
                  self.getType(), self.getParameters())
            
        NotCut = True
        if self.getType() == 0:
            Rad = np.sqrt(_R[0]**2 + _R[2]**2)
            if self.getDebug():
                print("     ----> Aperture cut: R, Raptr:", \
                      Rad, self.getParameters()[0])
            if Rad >= self.getParameters()[0]:
                NotCut = False
        elif self.getType() == 1:
            RadX2 = (_R[0]/self.getParameters()[0])**2
            RadY2 = (_R[2]/self.getParameters()[1])**2
            #print(" Aperture cut: RadX2, RapY2:", RadX2, RadY2)
            if (RadX2+RadY2) >= 1.:
                NotCut = False
        elif self.getType() == 2:
            lenX = abs(_R[0])
            lenY = abs(_R[2])
            #print(" Aperture cut: RadX2, RapY2:", RadX2, RadY2)
            if lenX > self.getParameters()[0] or \
               lenY > self.getParameters()[1]:
                NotCut = False

        if self.getDebug():
            print(" <---- Return, passsTHROUGH:", NotCut)
        return NotCut

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" Aperture(BeamLineElement).visualise: start")
            print("     ---->        CoordSys:", CoordSys)
            print("     ---->            Proj:", Proj)
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        xlim  = axs.get_xlim()
        ylim  = axs.get_ylim()

        ytot  = (ylim[1]-ylim[0])
        yzero = abs(ylim[0])/ytot
        if self.getDebug():
            print("     ----> xlim, ylim:", xlim, ylim)
            print("     ----> ytot, yzero:", ytot, yzero)

        Typ  = self.getType()
        xmin = self.getParameters()[0]
        ymin = self.getParameters()[0]
        if Typ == 1 or Typ == 2:
            ymin = self.getParameters()[1]
        if self.getDebug():
            print("     ----> xmin, ymin:", xmin, ymin)
            
        if CoordSys == "RPLC": 
            if self.getDebug():
                print("     ----> RPLC:")
                
            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)
                
            sz  = sStrt

        elif CoordSys == "Lab": 
            if self.getDebug():
                print("     ----> Lab:")
                
            cntr  = self.getrStrt()
            if self.getDebug():
                print("         ----> Centre:", cntr)
            sz = cntr[2]

        if Proj.find("x") >= 0:
            xyup  = yzero + xmin/ytot
            xydn  = yzero - xmin/ytot
        elif Proj.find("y") >= 0:
            xyup  = yzero + ymin/ytot
            xydn  = yzero - ymin/ytot
                
        amax = min(1., xyup + 1.0/ytot)
        amin = max(0., xydn - 1.0/ytot)
            
        if self.getDebug():
            print("     ---->     sz:", sz)
            print("     ---->   xyup:", xyup)
            print("     ---->   xydn:", xydn)
            
        axs.axvline(sz, xyup, amax, \
                    color="black", \
                    linewidth=1, \
                    zorder=2)
        axs.axvline(sz, xydn, amin, \
                    color="black", \
                    linewidth=1, \
                    zorder=2)
                                   
        if self.getDebug():
            print(" <---- Aperture(BeamLineElement).visualise: ends.")


#--------  I/o methods:
    def writeElement(self, dataFILE):
        if self.getDebug():
            print(" Aperture(BeamLineElement).writeElement starts.")

        derivedCLASS = "Aperture"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Derived class:", bversion.decode('utf-8'))

        record = strct.pack(">i", self.getType())
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Type:", strct.unpack(">i", record))

        record = strct.pack(">i", len(self.getParameters()))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Number of paramters:", \
                  strct.unpack(">i", record))

        if self.getDebug():
            print("     ----> Write parameters:")
        for iPrm in range(len(self.getParameters())):
            record = strct.pack(">d", self.getParameters()[iPrm])
            dataFILE.write(record)
            if self.getDebug():
                print("         ----> iPrm, value:", \
                      iPrm, strct.unpack(">d", record))
        if self.getDebug():
            print("     <---- Done.")
            
        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print(" <---- Aperture(BeamLineElement).writeElement done.")
        
    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" Aperture(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read(4)
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True
            
        record = strct.unpack(">i", brecord)
        Type = record[0]
        if cls.getDebug():
            print("     ----> Type:", Type)

        brecord = dataFILE.read(4)
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True
            
        record = strct.unpack(">i", brecord)
        nPrm = record[0]
        if cls.getDebug():
            print("     ----> Number of parameters:", nPrm)

        Params = []
        for iPrm in range(nPrm):
            brecord = dataFILE.read((1*8))
            if brecord == b'':
                return True
        
            record  = strct.unpack(">d", brecord)
            var     = float(record[0])
            Params.append(var)

        if cls.getDebug():
            print("     ----> Parameters:", Params)

        return EoF, Type, Params


"""
Derived class FocusQuadrupole:
==============================

  FocusQuadrupole class derived from BeamLineElement to contain paramters
  for an F-quad.


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
  Calling arguments:
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.


  _TrnsMtrx : Transfer matrix:
         Input: Brho: float ... B*rho (=3.3356E-3) * p (MeV)
                      at init in __init__ Brho=1. is used


  Instance attributes to define quadrupole:
  -----------------------------------------
  _Length  : Length of quad, m
  _Strength: Strength of quad in T/m

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising quad
                parameterrs.

  Set methods:
    setLength: set length:
          Input: _Length (float): length of quad (m)
  setStrength: _Strength (float): strength (gradient) of quad (T/m)
          Input: _Strength (float)
       setkFQ: _kFQ (float): quad focusing, k (m**-2)
          Input: _kFQ (float)

setTransferMatrix: Set transfer matrix; calculate using i/p brhop
          Input: Brho (T m)
         Return: np.array(6,6,) transfer matrix

  Get methods:
      getLength, getStrength

  Utilities:
    Transport: transport through focus quad.  Sets transfer matrix (using
               call to setTransferMatrix) given Brho
            Input:
                  R: 6D numpy array containing phase space at entry of
                      aperture
               Brho: Brho: float ... B*rho (~3.3356E-3) * p (MeV)
             Rprime: 6D phase space (numpy array) at exit of aperture
                     if the particle passes through the aperture ... or ...

"""
class FocusQuadrupole(BeamLineElement):
    instances = []
    __Debug   = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Length=None, _Strength=None, _kFQ=None):

        if self.getDebug():
            print(' FocusQuadrupole.__init__: ', \
                  'creating the FocusQuadrupole object: Length=', \
                  _Length, ', Strength=', _Strength, ', kFQ=', _kFQ)

        FocusQuadrupole.instances.append(self)

        self.setAll2None()

        #.. Hard wired:
        #   - FQmode = 0 ==> use particle momentum in calculation of k
        #            = 1 ==> use reference particle momentum and dispersion
        #                    calculation.
        self.setFQmode(0)

        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, _Name, _rStrt, _vStrt, _drStrt, _dvStrt)

        if not isinstance(_Length, float):
            raise badBeamLineElement("FocusQuadrupole: bad specification", \
                                     " for length!")
        if not isinstance(_Strength, float) and \
           not isinstance(_kFQ, float):
            raise badBeamLineElement("FocusQuadrupole: bad specification", \
                                     " for quadrupole strength!")

        self.setLength(_Length)
        if isinstance(_Strength, float):
            self.setStrength(_Strength)
            self.setkFQ(self.calckFQ())
        if isinstance(_kFQ, float):
            self.setkFQ(_kFQ)
            self.setStrength(self.calcStrength())

        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        if self.getDebug():
            print("     ----> New FocusQuadrupole instance: \n", self)
            print(" <---- Done.")

    def __repr__(self):
        return "FocusQuadrupole()"

    def __str__(self):
        print(" FocusQuadrupole:")
        print(" ----------------")
        print("     ---->     Debug flag:", FocusQuadrupole.getDebug())
        print("     ---->         FQmode:", self.getFQmode())
        print("     ---->     Length (m):", self.getLength())
        print("     ----> Strength (T/m):", self.getStrength())
        print("     ---->       kFQ (/m):", self.getkFQ())
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- FocusQuadrupole parameter dump complete."

    def SummaryStr(self):
        Str  = "FocusQuadrupole  : " + BeamLineElement.SummaryStr(self) + \
            "; Length = " + str(self.getLength()) + \
            "; Strength (gradient) = " + str(self.getStrength()) + \
            "; kFQ = " + str(self.getkFQ())

        return Str

    
# -------- "Set methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug
        
    def setAll2None(self):
        self._Length   = None
        self._Strength = None
        self._kFQ      = None
        self._TrnsMtrx = None
        self._FQmode   = None
        
    def setFQmode(self, _FQmode):
        if not isinstance(_FQmode, int):
            raise badParameter( \
                            "BeamLineElement.FocusQuadrupole.setFQmode:", \
                            " bad FQmode:", _FQmode)
        self._FQmode = _FQmode

    def setLength(self, _Length):
        if not isinstance(_Length, float):
            raise badParameter( \
                            "BeamLineElement.FocusQuadrupole.setLength:", \
                            " bad length:", _Length)
        self._Length = _Length

    def setStrength(self, _Strength):
        if not isinstance(_Strength, float):
            raise badParameter( \
                    "BeamLineElement.FocusQuadrupole.setStrength:", \
                    " bad quadrupole strength:", _Strength)
        self._Strength = _Strength

    def setkFQ(self, _kFQ):
        if not isinstance(_kFQ, float):
            raise badParameter( \
                    "BeamLineElement.FocusQuadrupole.setStrength:", \
                                " bad quadrupole k constant:", _kFQ)
        self._kFQ = _kFQ

    def setTransferMatrix(self, _R):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        E0        = iRefPrtcl.getPrOut()[iPrev][3]
        b0        = p0/E0
        b02       = b0**2
        g02       = 1./(1.-b02)
        
        if self.getDebug():
            print(" FocusQuadrupole(BeamLineElement).setTransferMatrix:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0, E0:", p0, E0)
            print("     <---- b0, b02, g02:", b0, b02, g02)

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Trace space:", _R)

        D   = 1.
        Scl = 1.
        if self.getFQmode() == 1:
            D = mth.sqrt(1. + 2.*_R[5]/b0 + _R[5]**2)
        else:
            particleMASS = PhysCnst.PhysicalConstants().getparticleMASS(   \
                                                   iRefPrtcl.getSpecies() )
            E   = E0 + _R[5]*p0
            p   = mth.sqrt(E**2 - particleMASS**2)
            if p > 0:
                Scl = p0 / p

        if self.getDebug():
            print("     ----> D:", D)
            print("     ----> E, p, Scl:", E, p, Scl)

        k = self.getkFQ() * Scl
        l = self.getLength()

        b = mth.sqrt(k/D)
        a = l * b
        b = b * D

        if self.getDebug():
            print("     ----> Length, Strength, kFQ:", \
                  self.getLength(), self.getStrength(), self.getkFQ())
            print("     ----> omegaPrime*L, omegaPrime*D:", a, b)

        TrnsMtrx = np.array([                                             \
            [   np.cos(a), np.sin(a)/b, 0., 0.,                   0., 0.],\
            [-b*np.sin(a),   np.cos(a), 0., 0.,                   0., 0.],\
            [          0.,          0., np.cosh(a), np.sinh(a)/b, 0., 0.],\
            [          0.,          0., b*np.sinh(a), np.cosh(a), 0., 0.],\
            [          0.,          0.,           0., 0.,  1., l/b02/g02],\
            [          0.,          0.,           0., 0.,  0.,        1.]\
                            ])

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

        
# -------- "Get methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug

    def getFQmode(self):
        return self._FQmode
    
    def getLength(self):
        return self._Length

    def getStrength(self):
        return self._Strength

    def getkFQ(self):
        return self._kFQ

    
# -------- Utilities:
    def calckFQ(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        
        if self.getDebug():
            print(" FocusQuadrupole(BeamLineElement).calckFQ:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0:", p0)

        Brho = (1./(speed_of_light*1.E-9))*p0/1000.
        kFQ  = self.getStrength() / Brho * \
            PhysCnst.PhysicalConstants().getparticleCHARGE(   \
                                                   iRefPrtcl.getSpecies() )

        if self.getDebug():
            print("     <---- kFQ:", kFQ)
            print(" <---- Done.")

        return kFQ

    def calcStrength(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        
        if self.getDebug():
            print(" FocusQuadrupole(BeamLineElement).calcStrength:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0:", p0)

        Brho = (1./(speed_of_light*1.E-9))*p0/1000.
        Strn = self.getkFQ() * Brho

        if self.getDebug():
            print("     <---- Strength:", Strn)
            print(" <---- Done.")

        return Strn
    
    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" FocusQuadrupole(BeamLineElement).visualise: start")
            print("     ----> CoordSys, Proj:", CoordSys, Proj)
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        xlim = axs.get_xlim()
        ylim = axs.get_ylim()

        wdth = self.getLength()
        hght = ylim[1]/2.
        
        angl = 0.
        abt  = 'xy'
        
        if CoordSys == "RPLC":
            if self.getDebug():
                print("     ----> RPLC:")
                
            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)
                
            if Proj == "xs":
                sxy   = [sStrt, 0.]
            elif Proj == "ys":
                sxy   = [sStrt, -hght]
                
        elif CoordSys == "Lab":
            if self.getDebug():
                print("     ----> RPLC:")
            
            BndPln = "xz"
            if abs(self.getStrt2End()[0]) < abs(self.getStrt2End()[1]):
                BndPln = "yz"
            if self.getDebug():
                print("         ----> Bending plane:", BndPln)
                
            strt = self.getrStrt()

            iCoord = 0
            if Proj.find("y") >=0: iCoord = 1
            
            if self.getDebug():
                print("         ----> strt:", strt)
                
            if Proj != BndPln:
                if Proj == "xz":
                    sxy   = [strt[2], 0.]
                elif Proj == "yz":
                    sxy   = [strt[2], -hght]
                
            elif Proj == BndPln:
                bbox = axs.get_window_extent()
                xax, yax = bbox.width, bbox.height
                xl = xlim[1] - xlim[0]
                yl = ylim[1] - ylim[0]
                xscl = xax*yl / (yax*xl)
                if self.getDebug():
                    print("         ----> xax, yax:", xax, yax)
                    print("         ---->   xy, yl:", xl, yl)
                    print("         ---->     xscl:", xscl)

                qVec1Lab  = self.getStrt2End()
                invRot    = np.linalg.inv(self.getRot2LbStrt())
                qVec1RPLC = np.matmul(invRot, qVec1Lab)
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("         ---->  qVec1Lab:", qVec1Lab)
                        print("     ----> InvRot: \n", invRot)
                        print("         ----> qVec1RPLC:", qVec1RPLC)
                
                hght = min(0.4, yl/4./xscl)
                if Proj == "xz":
                    qVec2RPLC = np.array([hght, 0., 0.])
                elif Proj == "yz":
                    qVec2RPLC = np.array([hght, 0., 0.])
                
                qVec2Lab = np.matmul(self.getRot2LbStrt(), qVec2RPLC)
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("         ----> qVec2RPLC:", qVec2RPLC)
                        print("     ----> Rot: \n", self.getRot2LbStrt())
                        print("         ---->  qVec2Lab:", qVec2Lab)

                qVec2Lab[iCoord] = qVec2Lab[iCoord] * xscl
                qVec2Lab[2]      = qVec2Lab[2]      / xscl
                
                crnr1 = self.getrStrt()
                crnr2 = crnr1 + qVec2Lab
                crnr3 = crnr2 + qVec1Lab
                crnr4 = crnr3 - qVec2Lab
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("         ----> crnr1:", crnr1)
                        print("         ----> crnr2:", crnr2)
                        print("         ----> crnr3:", crnr3)
                        print("         ----> crnr4:", crnr4)

                x = np.array([])
                y = np.array([])

                x = np.append(x, crnr1[2])
                x = np.append(x, crnr2[2])
                x = np.append(x, crnr3[2])
                x = np.append(x, crnr4[2])
                
                y = np.append(y, crnr1[iCoord])
                y = np.append(y, crnr2[iCoord])
                y = np.append(y, crnr3[iCoord])
                y = np.append(y, crnr4[iCoord])

                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("     ----> Rot: \n", self.getRot2LbStrt())
                        print("     ----> Inv: \n", invRot)
                        print("     ----> Angl:", angl)

                sxy  = [strt[2], strt[iCoord]]
                        
        if self.getDebug():
            print("     ---->   sxy:", sxy)
            print("     ---->  wdth:", wdth)
            print("     ---->  hght:", hght)
            print("     ---->  angl:", angl)
            print("     ---->   abt:", abt)

        if CoordSys == "RPLC" or \
           (CoordSys == "Lab" and Proj != BndPln):

            Patch = patches.Rectangle(sxy, wdth, hght, \
                                      angle=angl, \
                                      rotation_point=abt, \
                                      facecolor=('darkgreen'), \
                                      zorder=2)
            axs.add_patch(Patch)
        else:
            axs.fill(x, y, \
                     "darkgreen", \
                     zorder=2)
                                   
        if self.getDebug():
            print(" <---- FocusQuadrupole(BeamLineElement).visualise: ends.")

        
#--------  I/o methods:
    def getLines(self):
        Lines = []

        Fields  = self.getName().split(":")
        
        Stage   = Fields[1]
        Section = Fields[2]
        Element = "Fquad"
        Type    = ""
        
        Param   = "Length"
        Value   = self.getLength()
        Unit    = "m"
        Comment = ""
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        Param   = "Strength"
        Value   = self.getStrength()
        Unit    = "m"
        Comment = ""
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        return Lines
    
    def writeElement(self, dataFILE):
        if self.getDebug():
            print( \
                " FocusQuadrupole(BeamLineElement).writeElement starts.")

        derivedCLASS = "FocusQuadrupole"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Derived class:", bversion.decode('utf-8'))

        record = strct.pack(">3d", \
                            self.getLength(), \
                            self.getStrength(), \
                            self.getkFQ()
                            )

        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length, strength, kFQ, drStrt:", \
                  strct.unpack(">3d",record))

        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print( \
             " <---- FocusQuadrupole(BeamLineElement).writeElement done.")

    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" FocusQuadrupole(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read((3*8))
        record = strct.unpack(">3d", brecord)
       
        if brecord == b'':
            return True
        
        Ln  = None
        St  = None
        kFQ = None
        drStrt = np.array([0., 0., 0.])

        if float(record[0]) != -1.:
            Ln     = float(record[0])
        if float(record[1]) != -1.:
            St     = float(record[1])
        if float(record[2]) != -1.:
            kFQ   = float(record[2])

        if cls.getDebug():
            print("     ----> Ln, St, kDQ, drStrt:", Ln, St, kFQ, drStrt)
            
        return EoF, Ln, St, kFQ

    
"""
Derived class DefocusQuadrupole:
================================

  DefocusQuadrupole class derived from BeamLineElement to contain paramters
  for a D-quad.


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
  Calling arguments:
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.

  _TrnsMtrx : Transfer matrix:
         Input: Brho: float ... B*rho (=3.3356E-3) * p (MeV)
                      at init in __init__ Brho=1. is used


  Instance attributes to define quadrupole:
  -----------------------------------------
  _Length  : Length of quad, m
  _Strength: Strength of quad in T/m

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising quad
                parameterrs.

  Set methods:
    setLength: set length:
          Input: _Length (float): length of quad (m)
  setStrength: -Strength (float): strength of quad (T/m)
          Input: _Length (float): length of quad

setTransferMatrix: Set transfer matrix; calculate using i/p brhop
          Input: Brho (T m)
         Return: np.array(6,6,) transfer matrix

  Get methods:
      getLength, getStrength

  Utilities:
    Transport: transport through focus quad.  Sets transfer matrix (using
               call to setTransferMatrix) given Brho
            Input:
                  R: 6D numpy array containing phase space at entry of
                      aperture
               Brho: Brho: float ... B*rho (~3.3356E-3) * p (MeV
             Rprime: 6D phase space (numpy array) at exit of aperture
                     if the particle passes through the aperture ... or ...

"""
class DefocusQuadrupole(BeamLineElement):
    instances = []
    __Debug   = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Length=None, _Strength=None, _kDQ=None):

        if self.getDebug():
            print(' DefocusQuadrupole.__init__: ', \
                  'creating the DefocusQuadrupole object: Length=', \
                  _Length, ', Strength=', _Strength)

        DefocusQuadrupole.instances.append(self)

        self.setAll2None()

        #.. Hard wired:
        #   - DQmode = 0 ==> use particle momentum in calculation of k
        #            = 1 ==> use reference particle momentum and dispersion
        #                    calculation.
        self.setDQmode(0)
        
        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, _Name, \
                                 _rStrt, _vStrt, _drStrt, _dvStrt)

        if not isinstance(_Length, float):
            raise badBeamLineElement("DefocusQuadrupole:",\
                                     " bad specification for length!")
        if not isinstance(_Strength, float) and \
           not isinstance(_kDQ, float):
            raise badBeamLineElement("FocusQuadrupole: bad specification", \
                                     " for quadrupole strength!")

        self.setLength(_Length)
        if isinstance(_Strength, float):
            self.setStrength(_Strength)
            self.setkDQ(self.calckDQ())
        if isinstance(_kDQ, float):
            self.setkDQ(_kDQ)
            self.setStrength(self.calcStrength())
                
        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        if self.getDebug():
            print("     ----> New DefocusQuadrupole instance: \n", self)
            print(" <---- Done.")

    def __repr__(self):
        return "DefocusQuadrupole()"

    def __str__(self):
        print(" DefocusQuadrupole:")
        print(" -------------------")
        print("     ---->     Debug flag:", DefocusQuadrupole.getDebug())
        print("     ---->         DQmode:", self.getDQmode())
        print("     ---->     Length (m):", self.getLength())
        print("     ----> Strength (T/m):", self.getStrength())
        print("     ---->       kDQ (/m):", self.getkDQ())
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- DefocusQuadrupole parameter dump complete."

    def SummaryStr(self):
        Str  = "DefocusQuadrupole: " + BeamLineElement.SummaryStr(self) + \
            "; Length = " + str(self.getLength()) + \
            "; Strength = " + str(self.getStrength()) + \
            "; kDQ = " + str(self.getkDQ())

        return Str

    
# -------- "Set methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug
        
    def setAll2None(self):
        self._Length   = None
        self._Strength = None
        self._kDQ      = None
        self._TrnsMtrx = None
        self._DQmode   = None
        
    def setDQmode(self, _DQmode):
        if not isinstance(_DQmode, int):
            raise badParameter( \
                            "BeamLineElement.FocusQuadrupole.setDQmode:", \
                            " bad DQmode:", _DQmode)
        self._DQmode = _DQmode

    def setLength(self, _Length):
        if not isinstance(_Length, float):
            raise badParameter( \
                "BeamLineElement.DefocusQuadrupole.setLength:", \
                " bad length:", _Length)
        self._Length = _Length

    def setStrength(self, _Strength):
        if not isinstance(_Strength, float):
            raise badParameter( \
                "BeamLineElement.DefocusQuadrupole.setStrength:", \
                " bad quadrupole strength:", _Strength)
        self._Strength = _Strength

    def setkDQ(self, _kDQ):
        if not isinstance(_kDQ, float):
            raise badParameter( \
                    "BeamLineElement.DefocusQuadrupole.setkDQ:", \
                                " bad quadrupole k constant:", _kDQ)
        self._kDQ = _kDQ

    def setTransferMatrix(self, _R):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        E0        = iRefPrtcl.getPrOut()[iPrev][3]
        b0        = p0/E0
        b02       = b0**2
        g02       = 1./(1.-b02)
        
        if self.getDebug():
            print(" DefocusQuadrupole(BeamLineElement).setTransferMatrix:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0, E0:", p0, E0)
            print("     <---- b0, b02, g02:", b0, b02, g02)

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Trace space:", _R)

        D   = 1.
        Scl = 1.
        if self.getDQmode() == 1:
            D = mth.sqrt(1. + 2.*_R[5]/b0 + _R[5]**2)
        else:
            E   = E0 + _R[5]*p0
            particleMASS = PhysCnst.PhysicalConstants().getparticleMASS(   \
                                                   iRefPrtcl.getSpecies() )
            p   = mth.sqrt(E**2 - particleMASS**2)
            if p > 0:
                Scl = p0 / p
        
        if self.getDebug():
            print("     ----> D:", D)
            print("     ----> E, p, Scl:", E, p, Scl)

        k = self.getkDQ() * Scl
        l = self.getLength()

        b = mth.sqrt(k/D)
        a = l * b
        b = b * D

        if self.getDebug():
            print("     ----> Length, Strength, kDQ:", \
                  self.getLength(), self.getStrength(), self.getkDQ())
            print("     ----> omegaPrime*L, omegaPrime*D:", a, b)

        TrnsMtrx = np.array([                                               \
            [  np.cosh(a), np.sinh(a)/b,           0.,          0., 0., 0.],\
            [b*np.sinh(a),   np.cosh(a),           0.,          0., 0., 0.],\
            [          0.,           0.,    np.cos(a), np.sin(a)/b, 0., 0.],\
            [          0.,           0., -b*np.sin(a),   np.cos(a), 0., 0.],\
            [          0.,           0.,           0.,  0.,  1., l/b02/g02],\
            [          0.,           0.,           0.,  0.,  0.,        1.]\
        ])

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

        
# -------- "Get methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug
    
    def getDQmode(self):
        return self._DQmode
    
    def getLength(self):
        return self._Length

    def getStrength(self):
        return self._Strength

    def getkDQ(self):
        return self._kDQ

    
# -------- Utilities:
    def calckDQ(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        
        if self.getDebug():
            print(" DefocusQuadrupole(BeamLineElement).calckDQ:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0:", p0)

        Brho = (1./(speed_of_light*1.E-9))*p0/1000.
        kDQ  = self.getStrength() / Brho * \
            PhysCnst.PhysicalConstants().getparticleCHARGE(   \
                                                   iRefPrtcl.getSpecies() )

        if self.getDebug():
            print("     <---- kDQ:", kDQ)
            print(" <---- Done.")

        return kDQ

    def calcStrength(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        
        if self.getDebug():
            print(" DefocusQuadrupole(BeamLineElement).calcStrength:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0:", p0)

        Brho = (1./(speed_of_light*1.E-9))*p0/1000.
        Strn = self.getkDQ() * Brho  / \
            PhysCnst.PhysicalConstants().getparticleCHARGE(   \
                                                   iRefPrtcl.getSpecies() )

        if self.getDebug():
            print("     <---- Strength:", Strn)
            print(" <---- Done.")

        return Strn

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" DefocusQuadrupole(BeamLineElement).visualise: start")
            print("     ----> CoordSys, Proj:", CoordSys, Proj)
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        xlim = axs.get_xlim()
        ylim = axs.get_ylim()

        wdth = self.getLength()
        hght = min(0.4, abs(ylim[0])/2.)
        
        angl = 0.
        abt  = 'xy'

        if CoordSys == "RPLC":
            if self.getDebug():
                print("     ----> RPLC:")
                
            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)

            if Proj == "xs":
                sxy   = [sStrt, -hght]
            elif Proj == "ys":
                sxy   = [sStrt, 0.]

        elif CoordSys == "Lab":
            if self.getDebug():
                print("     ----> RPLC:")
                
            BndPln = "xz"
            if abs(self.getStrt2End()[0]) < abs(self.getStrt2End()[1]):
                BndPln = "yz"
            if self.getDebug():
                print("         ----> Bending plane:", BndPln)
                
            strt = self.getrStrt()
            
            iCoord = 0
            if Proj.find("y") >=0: iCoord = 1
            
            if self.getDebug():
                print("         ---->  strt:", strt)
                print("         ----> iCoord:", iCoord)
                
            if Proj != BndPln:
                if Proj == "xz":
                    sxy   = [strt[2], -hght]
                elif Proj == "yz":
                    sxy   = [strt[2], 0.]

            elif Proj ==BndPln:
                bbox = axs.get_window_extent()
                xax, yax = bbox.width, bbox.height
                xl = xlim[1] - xlim[0]
                yl = ylim[1] - ylim[0]
                xscl = xax*yl / (yax*xl)
                if self.getDebug():
                    print("         ----> xax, yax:", xax, yax)
                    print("         ---->   xy, yl:", xl, yl)
                    print("         ---->     xscl:", xscl)

                qVec1Lab  = self.getStrt2End()
                invRot    = np.linalg.inv(self.getRot2LbStrt())
                qVec1RPLC = np.matmul(invRot, qVec1Lab)
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("         ---->  qVec1Lab:", qVec1Lab)
                        print("     ----> InvRot: \n", invRot)
                        print("         ----> qVec1RPLC:", qVec1RPLC)
                
                hght = min(0.4, yl/4./xscl)
                if Proj == "xz":
                    qVec2RPLC = np.array([-hght, 0., 0.])
                elif Proj == "yz":
                    qVec2RPLC = np.array([-hght, 0., 0.])
                
                qVec2Lab = np.matmul(self.getRot2LbStrt(), qVec2RPLC)
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("         ----> qVec2RPLC:", qVec2RPLC)
                        print("     ----> Rot: \n", self.getRot2LbStrt())
                        print("         ---->  qVec2Lab:", qVec2Lab)

                qVec2Lab[iCoord] = qVec2Lab[iCoord] * xscl
                qVec2Lab[2]      = qVec2Lab[2]      / xscl
                
                crnr1 = self.getrStrt()
                crnr2 = crnr1 + qVec2Lab
                crnr3 = crnr2 + qVec1Lab
                crnr4 = crnr3 - qVec2Lab
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("         ----> crnr1:", crnr1)
                        print("         ----> crnr2:", crnr2)
                        print("         ----> crnr3:", crnr3)
                        print("         ----> crnr4:", crnr4)

                x = np.array([])
                y = np.array([])

                x = np.append(x, crnr1[2])
                x = np.append(x, crnr2[2])
                x = np.append(x, crnr3[2])
                x = np.append(x, crnr4[2])
                
                y = np.append(y, crnr1[iCoord])
                y = np.append(y, crnr2[iCoord])
                y = np.append(y, crnr3[iCoord])
                y = np.append(y, crnr4[iCoord])
                
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7,\
                                         suppress=True):
                        print("     ----> Rot: \n", self.getRot2LbStrt())
                        """
                        print("     ----> Rot: \n", invRot)
                        """
                        print("     ----> Angl:", angl)

                sxy  = [strt[2], strt[iCoord]]
                        
        if self.getDebug():
            if CoordSys == "Lab" and Proj == BndPln:
                print("     ---->    x, y:", x, y)
            else:
                print("     ---->    sxy:", sxy)
            print("     ---->   wdth:", wdth)
            print("     ---->   hght:", hght)
            print("     ---->   angl:", angl)
            print("     ---->    abt:", abt)

        if CoordSys == "RPLC" or \
           (CoordSys == "Lab" and Proj != BndPln):
    
            Patch = patches.Rectangle(sxy, wdth, hght, \
                                      angle=angl, \
                                      rotation_point=abt, \
                                      facecolor=('darkgreen'), \
                                      zorder=2)
            axs.add_patch(Patch)
        else:
            axs.fill(x, y, \
                     "darkgreen", \
                     zorder=2)
                                   
        if self.getDebug():
            print(" <---- DefocusQuadrupole(BeamLineElement).visualise: ends.")
        
#--------  I/o methods:
    def getLines(self):
        Lines = []

        Fields  = self.getName().split(":")
        
        Stage   = Fields[1]
        Section = Fields[2]
        Element = "Dquad"
        Type    = ""
        
        Param   = "Length"
        Value   = self.getLength()
        Unit    = "m"
        Comment = ""
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        Param   = "Strength"
        Value   = self.getStrength()
        Unit    = "m"
        Comment = ""
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        return Lines
    
    def writeElement(self, dataFILE):
        if self.getDebug():
            print( \
                " DefocusQuadrupole(BeamLineElement).writeElement starts.")

        derivedCLASS = "DefocusQuadrupole"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Derived class:", bversion.decode('utf-8'))

        record = strct.pack(">3d", \
                            self.getLength(), \
                            self.getStrength(), \
                            self.getkDQ()
                            )

        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length, strength, kDQ:", \
                  strct.unpack(">3d",record))

        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print( \
             " <---- DefocusQuadrupole(BeamLineElement).writeElement done.")

    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" DefocusQuadrupole(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read((3*8))
        record = strct.unpack(">3d", brecord)
       
        if brecord == b'':
            return True

        Ln  = None
        St  = None
        kDQ = None
        drStrt = np.array([0., 0., 0.])

        if float(record[0]) != -1.:
            Ln     = float(record[0])
        if float(record[1]) != -1.:
            St     = float(record[1])
        if float(record[2]) != -1.:
            kDQ   = float(record[2])

        if cls.getDebug():
            print("     ----> Ln, St, kDQ:", Ln, St, kDQ)
            
        return EoF, Ln, St, kDQ

    
"""
Derived class SectorDipole:
===========================

  SectorDipole class derived from BeamLineElement to contain parameters
  for a sector dipole


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
  Calling arguments:
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.

  _TrnsMtrx : Transfer matrix


  Instance attributes to define quadrupole:
  -----------------------------------------
  _Angle  : Bending angle, degrees

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising dipole
                parameterrs.

  Set methods:
    setAngle: set length:
          Input: _Angle (float): angle through which dipole bends reference
                                 particle

setTransferMatrix: Set transfer matrix; calculate using i/p brhop
          Input: Brho (T m)
         Return: np.array(6,6,) transfer matrix

  Get methods:
      getAngle

  Utilities:

"""
class SectorDipole(BeamLineElement):
    instances = []
    __Debug = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Angle=None, _B=None):
        if self.getDebug():
            print(' SectorDipole(BeamLineElement).__init__: ')
            print('     ----> Angle=', _Angle)
            print('     ---->     B=', _B)

        SectorDipole.instances.append(self)

        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, _Name, \
                                 _rStrt, _vStrt, _drStrt, _dvStrt)

        if not isinstance(_Angle, float):
            raise badBeamLineElement( \
                "SectorDipole: bad specification for bending angle (Angle)!")

        self.setAngle(_Angle)
        self.setB(_B)
        self.setLength()

        ChalfAngle = mth.cos(self.getAngle()/2.)
        ShalfAngle = mth.sin(self.getAngle()/2.)
        Radius     = self.getLength()/self.getAngle()
        Scale      = 2. * Radius * ShalfAngle
        deltaV     = np.array([-ShalfAngle, 0., ChalfAngle]) * Scale
        self.setStrt2End(deltaV)
        
        CAngle = mth.cos(self.getAngle())
        SAngle = mth.sin(self.getAngle())
        Rot2End    = np.array([ \
                                [CAngle, 0., -SAngle], \
                                [    0., 1.,      0.], \
                                [SAngle, 0.,  CAngle]  \
                                ])
        Rot2LbEnd  = np.matmul(self.getRot2LbStrt(), Rot2End)
        self.setRot2LbEnd(Rot2LbEnd)
        
        if self.getDebug():
            print("     ----> New SectorDipole instance: \n", self)

    def __repr__(self):
        return "SectorDipole()"

    def __str__(self):
        print(" SectorDipole:")
        print(" -------")
        print("     ----> Debug flag:", SectorDipole.getDebug())
        print("     ----> Bending Angle (Angle):", self.getAngle())
        print("     ----> Magnetic field:", self.getB())
        print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- SectorDipole parameter dump com/plete."

    def SummaryStr(self):
        Str  = "Sector dipole    : " + BeamLineElement.SummaryStr(self) + \
            "; Angle = " + str(self.getAngle()) + \
            "; B = " + str(self.getB())
        return Str

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" SectorDipole(BeamLineElement).visualise: start")
            print("     ----> CoordSys, Proj:", CoordSys, Proj)
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        strt = self.getrStrt()
        xlim = axs.get_xlim()
        ylim = axs.get_ylim()

        wdth = self.getLength()
        hght = 0.4
        angl = 0.
        abt  = 'center'
            
        if CoordSys == "RPLC":
            if self.getDebug():
                print("     ----> RPLC:")
                
            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)

            hght = (ylim[1] - ylim[0]) / 10.
            sxy   = [ sStrt, -hght/2. ]

        elif CoordSys == "Lab":
            if self.getDebug():
                print("     ----> Lab:")

            BndPln = "xz"
            if abs(self.getStrt2End()[0]) < abs(self.getStrt2End()[1]):
                BndPln = "yz"
            if self.getDebug():
                print("         ----> Bending plane:", BndPln)

            iCoord = 0
            if Proj == "yz": iCoord = 1
                
            if Proj != BndPln:
                wdth = self.getStrt2End()[2]
                hght = (ylim[1] - ylim[0]) / 10.
                sxy   = [ strt[2], strt[iCoord]-hght/2. ]
                
            elif Proj == BndPln:
                rRPLC = np.array([-self.getRadius(), 0., 0.])
                rLab   = np.matmul(self.getRot2LbStrt(), rRPLC)
                cntr   = strt + rLab
                vec    = strt - rLab
                theta1 = mth.atan2(vec[iCoord]-strt[iCoord], vec[2]-strt[2])
                if theta1 < 0.:
                    theta1 += 2.*mth.pi
                theta2 = theta1 + \
                    self.getAngle()*mth.copysign(1.,rLab[iCoord])
                if theta1 < theta2:
                    theta = np.array([theta1*180./mth.pi, theta2*180./mth.pi])
                else:
                    theta = np.array([theta2*180./mth.pi, theta1*180./mth.pi])
                
                sxy   = np.array([cntr[2], cntr[iCoord]])
                wdth = 0.1
                
                if self.getDebug():
                    print("     ----> Lab:")
                    with np.printoptions(linewidth=500, \
                                         precision=7,suppress=True):
                        print("         ---->  rRPLC:", rRPLC)
                        print("         ---->   rLab:", rLab)
                        print("         ---->   cntr:", cntr)
                        print("         ---->    vec:", vec)
                        print("         ---->  theta:", theta)

        if self.getDebug():
            print("     ---->  Start:", strt)
            print("     ---->    sxy:", sxy)
            print("     ---->   wdth:", wdth)
            print("     ---->   hght:", hght)
            print("     ---->   angl:", angl)
            print("     ---->    abt:", abt)

        if CoordSys == "Lab":
            if Proj == BndPln:
                Patch = patches.Wedge(sxy, self.getRadius()+wdth, \
                                      theta[0], theta[1], \
                                      width=2.5*wdth, \
                                      facecolor=('mediumvioletred', 1.), \
                                      zorder=2)
            elif Proj != BndPln:
                Patch = patches.Rectangle(sxy, wdth, hght, \
                                          angle=angl, \
                                          rotation_point=abt, \
                                          facecolor=('mediumvioletred', 1.), \
                                          zorder=2)
        else:
            Patch = patches.Rectangle(sxy, wdth, hght, \
                                      angle=angl, \
                                      rotation_point=abt, \
                                      facecolor=('mediumvioletred', 1.), \
                                      zorder=2)

        axs.add_patch(Patch)
                                   
        if self.getDebug():
            print(" <---- SectorDipole(BeamLineElement).visualise: ends.")


# -------- "Set methods"
#..  Methods believed to be self-documenting(!)
    def setAngle(self, _Angle):
        if not isinstance(_Angle, float):
            raise badParameter(\
                               "BeamLineElement.SectorDipole.setAngle:", \
                               "bad bending angle (Angle):", _Angle)
        self._Angle = _Angle

    def setB(self, _B):
        if not isinstance(_B, float):
            raise badParameter(\
                               "BeamLineElement.SectorDipole.setB:", \
                               "bad B:", _B)
        self._B = _B

    def setLength(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0   = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))

        if self.getDebug():
            print(" Dipole(BeamLineElement).setLength:")
            print("     ----> Reference particle 4-mmtm:", \
                  iRefPrtcl.getPrIn()[0])
            print("         ----> p0:", p0)
        
        Brho = (1./(speed_of_light*1.E-9))*p0/1000.
        r    = Brho / self.getB() / \
            PhysCnst.PhysicalConstants().getparticleCHARGE(   \
                                               iRefPrtcl.getSpecies())
        l    = r * self.getAngle()

        if self.getDebug():
            print("     ----> Brho, r, l:", Brho, r, l)

        self._Length = l

    def setTransferMatrix(self, _R):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        E0        = iRefPrtcl.getPrOut()[iPrev][3]
        b0        = p0/E0
        b02       = b0**2
        g02       = 1./(1.-b02)
        
        if self.getDebug():
            print(" Dipole(BeamLineElement).setTransferMatrix:")
            print("     ----> Reference particle 4-mmtm:", \
                  iRefPrtcl.getPrIn()[0])
            print("     ----> NOdispersion:", \
                  BeamLineElement.getNOdispersion())
            print("         ----> p0, E0:", p0, E0)
            print("     <---- b02, g02:", b02, g02)
        
        E    = E0 + p0*_R[5]
        particleMASS = PhysCnst.PhysicalConstants().getparticleMASS(   \
                                               iRefPrtcl.getSpecies() )
        p    = mth.sqrt(E**2 - particleMASS**2)
        
        if self.getDebug():
            print("     ----> Particle energy and mmtm:", E, p)

        Brho = (1./(speed_of_light*1.E-9))*p/1000.
        r    = Brho / self.getB()
        c    = np.cos(self.getAngle())
        s    = np.sin(self.getAngle())
        l    = self.getLength()

        if self.getDebug():
            print("     ----> Brho, r, c, s, l:", Brho, r, c, s, l)

        TrnsMtrx = np.array([
            [     c,            r*s, 0., 0., 0.,                r*(1-c)/b0],
            [  -s/r,              c, 0., 0., 0.,                      s/b0],
            [    0.,             0., 1.,  l, 0.,                        0.],
            [    0.,             0., 0., 1., 0.,                        0.],
            [ -s/b0, -(r/b0)*(1.-c), 0., 0., 1., l/b02/g02 - (l-r*s)/b0**2],
            [    0.,             0., 0., 0., 0.,                        1.]
        ])

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug

# -------- "Get methods"
#..  Methods believed to be self-documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug

    def getAngle(self):
        return self._Angle

    def getB(self):
        return self._B

    def getLength(self):
        return self._Length

    def getRadius(self):
        return self.getLength() / self.getAngle()


#--------  I/o methods:
    def writeElement(self, dataFILE):
        if self.getDebug():
            print( \
                " Dipole(BeamLineElement).writeElement starts.")

        derivedCLASS = "SectorDipole"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Derived class:", bversion.decode('utf-8'))

        record = strct.pack(">2d", \
                            self.getAngle(), \
                            self.getB())

        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Angle, B:", \
                  strct.unpack(">2d",record))

        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print( \
             " <---- Dipole(BeamLineElement).writeElement done.")

        return

    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" SectorDipole(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read((2*8))
        if brecord == b'':
            return True
        
        record = strct.unpack(">2d", brecord)

        Angl = float(record[0])
        B    = float(record[1])

        if cls.getDebug():
            print("     ----> Angl, B:", Angl, B)
            
        return EoF, Angl, B

    
class Octupole(BeamLineElement):
    instances = []
    __Debug = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Length=None):
        if self.__Debug:
            print(' Octupole.__init__: ', 'creating the Octupole object: Length=', _Length)

        Octupole.instances.append(self)

        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, _Name, _rStrt, _vStrt, _drStrt, _dvStrt)

        if not isinstance(_Length, float):
            raise badBeamLineElement("Octupole: bad specification for length!")

        self.setLength(_Length)

        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())

        self.setTransferMatrix()

        if self.__Debug:
            print("     ----> New Octupole instance: \n", self)

    def __repr__(self):
        return "Octupole()"

    def __str__(self):
        print(" Octupole:")
        print(" ---------")
        print("     ----> Debug flag:", Octupole.getDebug())
        print("     ----> Length (m):", self.getLength())
        print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- Octupole parameter dump complete."

    # -------- "Set methods"
    # Methods believed to be self-documenting(!)

    def setLength(self, _Length):
        if not isinstance(_Length, float):
            raise badParameter("BeamLineElement.Octupole.setLength: bad length:", _Length)
        self._Length = _Length

    def setTransferMatrix(self):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        l = self._Length

        TrnsMtrx = np.array([
            [1, l, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, l, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]
        ])

        self._TrnsMtrx = TrnsMtrx

    # -------- "Get methods"
    # Methods believed to be self-documenting(!)

    def getLength(self):
        return self._Length
    
"""
Derived class Solenoid:
=======================

  Solenoid class derived from BeamLineElement to contain paramters
  for a solenoid.


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Derived instance attributes:
  ----------------------------
  Calling arguments:
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.


  _TrnsMtrx : Transfer matrix (6x6).  Set to Null in __init__, initialised
              (to Null) in BeamLineElement.__init__, filled in derived
              classes.
                

  Instance attributes to define quadrupole:
  -----------------------------------------
  _Length  : Length of quad, m
  _Strength: Strength of solenoid (B0) in T

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising quad
                parameterrs.

  Set methods:
    setLength: set length:
          Input: _Length (float)   : length of quad (m)
  setStrength: set primary field strength
          Input: _Streng=th (float): strength of solenoid (B0) (T)

setTransferMatrix: Set transfer matrix; calculate using i/p kinetic
                   energy for test particle and reference particle
          Input: Brho (T m)
         Return: np.array(6,6,) transfer matrix
         Return: np.array(6,6,) transfer matrix

  Get methods:
      getLength, getStrength

  Utilities:
    Transport: transport through solenoid.  Sets transfer matrix (using
               call to setTransferMatrix) given Brho
            Input:
                  R: 6D numpy array containing phase space at entry of
                      solenoid
               Brho: Brho (T m)
             Rprime: 6D phase space (numpy array) at exit of aperture
                     if the particle passes through the solenoid.

"""
class Solenoid(BeamLineElement):
    instances = []
    __Debug = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Length=None, _Strength=None, _ksol=None):
        if self.getDebug():
            print(" Solenoid.__init__:", \
                  " creating the Solenoid object: Length=", _Length, " m,"\
                  " strength=", _Strength, " T")

        Solenoid.instances.append(self)

        self.setAll2None()

        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, \
                                 _Name, _rStrt, _vStrt, _drStrt, _dvStrt)

        if not isinstance(_Length, float):
            raise badBeamLineElement( \
                            "Solenoid: bad specification for length!")
        if not isinstance(_Strength, float) and \
           not isinstance(_ksol, float):
            raise badBeamLineElement( \
                            "Solenoid: bad specification for strength!")

        self.setLength(_Length)
        if isinstance(_Strength, float):
            self.setStrength(_Strength)
            self.setksol(self.calcksol())
        if isinstance(_ksol, float):
            self.setksol(_ksol)
            self.setStrength(self.calcStrength())
                
        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        if self.getDebug():
            print("     ----> New Solenoid instance: \n", self)
            print(" <---- Done.")

    def __repr__(self):
        return "Solenoid()"

    def __str__(self):
        print(" Solenoid:")
        print(" ---------")
        print("     ----> Debug flag:", Solenoid.getDebug())
        print("     ----> Length (m):", self.getLength())
        print("     ----> Strength:", self.getStrength())
        print("     ---->ksol (/m):", self.getksol())
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- Solenoid parameter dump complete."

    def SummaryStr(self):
        Str  = "Solenoid         : " + BeamLineElement.SummaryStr(self) + \
            "; Length = " + str(self.getLength()) + \
            "; Strength = " + str(self.getStrength()) + \
            "; ksol = " + str(self.getksol())
        return Str

    
# -------- "Set methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug
        
    def setAll2None(self):
        self._Length   = None
        self._Strength = None
        self._ksol     = None
        self._TrnsMtrx = None
        
    def setLength(self, _Length):
        if not isinstance(_Length, float):
            raise badParameter( \
                "BeamLineElement.Solenoid.setLength: bad length:", \
                                _Length)
        self._Length = _Length

    def setStrength(self, _Strength):
        if not isinstance(_Strength, float):
            raise badParameter("BeamLineElement.Solenoid.setStrength:" \
                               " bad strength value:", \
                               _Strength)
        self._Strength = _Strength

    def setksol(self, _ksol):
        if not isinstance(_ksol, float):
            raise badParameter( \
                    "BeamLineElement.Solenloid.setcsol:", \
                                " bad quadrupole k constant:", _kDQ)
        self._ksol = _ksol

    def setTransferMatrix(self, _R):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        E0        = iRefPrtcl.getPrOut()[iPrev][3]

        b02       = (p0/E0)**2
        g02       = 1./(1.-b02)
        
        if self.getDebug():
            print(" Solenoid(BeamLineElement).setTransferMatrix:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrIn()[0])
            print("         ----> p0, E0:", p0, E0)
            print("     <---- b02, g02:", b02, g02)

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Trace space:", _R)

        E    = E0 + p0*_R[5]
        particleMASS = PhysCnst.PhysicalConstants().getparticleMASS(   \
                                               iRefPrtcl.getSpecies() )
        p    = mth.sqrt(E**2 - particleMASS**2)
        
        if self.getDebug():
            print("     ----> Particle energy and mmtm:", E, p)

        Brho = (1./(speed_of_light*1.E-9))*p/1000.
        l  = self.getLength()
        k  = self.getksol() * p0 / p
        
        ckl  = mth.cos(k*l)
        skl  = mth.sin(k*l)
        sckl = ckl*skl

        if self.getDebug():
            print("     ----> Length, Strength:", \
                  self.getLength(), self.getStrength())
            print("     ----> Brho, k, l:", Brho, k, l)

        TrnsMtrx = np.array([                                      \
            [   ckl**2,    sckl/k,      sckl, (skl**2)/k, 0., 0.], \
            [  -k*sckl,    ckl**2, -k*skl**2,       sckl, 0., 0.], \
            [    -sckl, -skl**2/k,    ckl**2,     sckl/k, 0., 0.], \
            [ k*skl**2,     -sckl,   -k*sckl,     ckl**2, 0., 0.], \
            [       0.,        0.,        0.,  0.,  1., l/b02/g02],\
            [       0.,        0.,        0.,  0.,  0.,        1.] \
                             ])

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" Solenoid(BeamLineElement).visualise: start")
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        xlim = axs.get_xlim()
        ylim = axs.get_ylim()

        wdth = self.getLength()
        angl = 0.
        abt  = 'center'
            
        hght = min(0.4, (ylim[1] - ylim[0]) / 10.)

        if CoordSys == "RPLC":
            if self.getDebug():
                print("     ----> RPLC:")

            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)

            sxy   = [ sStrt, -hght/2. ]
            
        elif CoordSys == "Lab":
            if self.getDebug():
                print("     ----> Lab:")
            
            cntr = self.getrStrt()
            if self.getDebug():
                print("         ----> Cntr:", Cntr)

            if Proj == "xz":
                sxy   = [ cntr[2], cntr[0]-hght/2. ]
            elif Proj == "yz":
                sxy   = [ cntr[2], cntr[1]-hght/2. ]

        if self.getDebug():
            print("     ----> Centre:", cntr)
            print("     ---->    sxy:", sxy)
            print("     ---->   wdth:", wdth)
            print("     ---->   hght:", hght)
            print("     ---->   angl:", angl)
            print("     ---->    abt:", abt)
            
        Patch = patches.Rectangle(sxy, wdth, hght, \
                                  angle=angl, \
                                  rotation_point=abt, \
                                  facecolor=('orange', 1.), \
                                  zorder=2)

        axs.add_patch(Patch)
                                   
        if self.getDebug():
            print(" <---- Solenoid(BeamLineElement).visualise: ends.")

            
# -------- Get methods:
#..   Methods believed to be self-documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug
    
    def getLength(self):
        return self._Length

    def getStrength(self):
        return self._Strength
    
    def getksol(self):
        return self._ksol


# -------- Utilities:
    def calcksol(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        
        if self.getDebug():
            print(" Solenoid(BeamLineElement).calckDQ:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0:", p0)

        Brho = (1./(speed_of_light*1.E-9))*p0/1000.
        if self.getDebug():
            print("         ----> Brho, B:", Brho, self.getStrength())

        ksol  = self.getStrength() / (2. * Brho) * \
            PhysCnst.PhysicalConstants().getparticleCHARGE(   \
                                                   iRefPrtcl.getSpecies() )
        if self.getDebug():
            print("     <---- ksol:", ksol)
            print(" <---- Done.")

        return ksol

    def calcStrength(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0        = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))
        
        if self.getDebug():
            print(" Solenoid(BeamLineElement).calcStrength:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Reference particle 4-mmtm:", \
                      iRefPrtcl.getPrOut()[iPrev])
            print("         ----> p0:", p0)

        Brho = (1./(speed_of_light*1.E-9))*p0/1000.
        Strn = self.getksol() * 2. * Brho

        if self.getDebug():
            print("     <---- Strength:", Strn)
            print(" <---- Done.")

        return Strn

    
#--------  I/o methods:
    def writeElement(self, dataFILE):
        if self.getDebug():
            print(" Solenoid(BeamLineElement).writeElement starts.")

        derivedCLASS = "Solenoid"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Derived class:", bversion.decode('utf-8'))

        record = strct.pack(">3d", \
                            self.getLength(), \
                            self.getStrength(), \
                            self.getksol())

        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length, strength, k_sol:", \
                  strct.unpack(">3d",record))

        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print(" <---- Solenoid(BeamLineElement).writeElement done.")

        return self._ksol
    
    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" Solenoid(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read((3*8))
        if brecord == b'':
            return True
        
        record = strct.unpack(">3d", brecord)

        Ln     = None
        St     = None
        kSol   = None

        if float(record[0]) != -1.:
            Ln     = float(record[0])
        if float(record[1]) != -1.:
            St     = float(record[1])
        if float(record[2]) != -1.:
            kSol   = float(record[2])

        if cls.getDebug():
            print("     ----> Ln, St, kSol:", Ln, St, kSol)
            
        return EoF, Ln, St, kSol


"""
Derived class GaborLens:
========================

  GaborLens class derived from BeamLineElement to contain paramters
  for a Gabor lens.


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Derived instance attributes:
  ----------------------------
  Calling arguments:
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.


  _TrnsMtrx : Transfer matrix (6x6).  Set to Null in __init__, initialised
              (to Null) in BeamLineElement.__init__, filled in derived
              classes.
                

  Instance attributes to define quadrupole:
  -----------------------------------------
  _Length  : Length of quad, m
  _Strength: Strength of solenoid (B0) in T

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising quad
                parameterrs.

  Set methods:
           setLength: set length:
                  Input: _Length (float)   : length of lens (m)
  setElectronDensity: set primary field strength
                  Input: _ne (float): electron density (m^{-3))

setTransferMatrix: Set transfer matrix; calculate using i/p kinetic
                   energy for test particle and reference particle
          Input: Brho (T m)
         Return: np.array(6,6,) transfer matrix
         Return: np.array(6,6,) transfer matrix

  Get methods:
      getLength, getElectronDensity

  Utilities:
    Transport: transport through solenoid.  Sets transfer matrix (using
               call to setTransferMatrix) given Brho
            Input:
                  R: 6D numpy array containing phase space at entry of
                      solenoid
               Brho: Brho (T m)
             Rprime: 6D phase space (numpy array) at exit of aperture
                     if the particle passes through the solenoid.

"""
class GaborLens(BeamLineElement):
    instances = []
    __Debug = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Bz=None, _VA=None, _RA=None, _Rp=None, _Length=None, \
                 _Strength=None):
        
        if self.getDebug():
            print(" GaborLens.__init__:", \
                  " creating the GaborLens object:")
            print("     ---->       Bz:",      _Bz, " T")
            print("     ---->       VA:",       _VA, " V")
            print("     ---->       RA:",       _RA, " m")
            print("     ---->       RP:",       _Rp, " m")
            print("     ---->   Length:",   _Length, " m")
            print("     ----> Strength:", _Strength, " m^{-2}")

        GaborLens.instances.append(self)

        OK = self.setAll2None()
        
        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, \
                                 _Name, _rStrt, _vStrt, _drStrt, _dvStrt)

        if not isinstance(_Length, float):
            raise badBeamLineElement( \
                            "GaborLens: bad specification for length!")
        self.setLength(_Length)

        if _Strength == None:
            if not isinstance(_Bz, float):
                raise badBeamLineElement( \
                                "GaborLens: bad specification for Bz!")
            if not isinstance(_VA, float):
                raise badBeamLineElement( \
                            "GaborLens: bad specification for Bz!")
            if not isinstance(_RA, float):
                raise badBeamLineElement( \
                            "GaborLens: bad specification for Bz!")
            if not isinstance(_Rp, float):
                raise badBeamLineElement( \
                            "GaborLens: bad specification for Bz!")

            self.setBz(_Bz)
            self.setVA(_VA)
            self.setRA(_RA)
            self.setRp(_Rp)
            
        else:
            if not isinstance(_Strength, float):
                raise badBeamLineElement( \
                            "GaborLens: bad specification for Strength!")
            self.setStrength(_Strength)
            
        self.setElectronDensity()

        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        if self.getDebug():
            print("     ----> New GaborLens instance: \n", self)

    def __repr__(self):
        return "GaborLens()"

    def __str__(self):
        print(" GaborLens:")
        print(" ---------")
        print("     ----> Debug flag:", GaborLens.getDebug())
        print("     ----> Length (m):", self.getLength())
        print("     ----> n_e (/m^3):", self.getElectronDensity())
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- GaborLens parameter dump complete."

    def SummaryStr(self):
        Str  = "GaborLens        : " + BeamLineElement.SummaryStr(self) + \
            "; Length = " + str(self.getLength()) + \
            "; Electron density = " + str(self.getElectronDensity())
        return Str

    
# -------- "Set methods"
# Methods believed to be self-documenting(!)
    def setAll2None(self):
        self._Bz       = None
        self._VA       = None
        self._RA       = None
        self._Rp       = None
        self._Length   = None
        self._Strength = None
        self._TrnsMtrx = None
        
    def setBz(self, _Bz):
        if not isinstance(_Bz, float):
            raise badParameter( \
                "BeamLineElement.GaborLens.setBz: bad length:", _Bz)
        self._Bz = _Bz

    def setVA(self, _VA):
        if not isinstance(_VA, float):
            raise badParameter( \
                "BeamLineElement.GaborLens.setVA: bad length:", _VA)
        self._VA = _VA

    def setRA(self, _RA):
        if not isinstance(_RA, float):
            raise badParameter( \
                "BeamLineElement.GaborLens.setRA: bad length:", _RA)
        self._RA = _RA

    def setRp(self, _Rp):
        if not isinstance(_Rp, float):
            raise badParameter( \
                "BeamLineElement.GaborLens.setRp: bad length:", _Rp)
        self._Rp = _Rp

    def setLength(self, _Length):
        if not isinstance(_Length, float):
            raise badParameter( \
                "BeamLineElement.GaborLens.setLength: bad length:", _Length)
        self._Length = _Length

    def setStrength(self, _Strength):
        if not isinstance(_Strength, float):
            raise badParameter( \
                "BeamLineElement.GaborLens.setLength: bad strength:", \
                                _Strength)
        self._Strength = _Strength

    def setElectronDensity(self):
        if self.getDebug():
            print(" GaborLens(BeamLineElement).setElectronDensity:")
        if isinstance(self.getStrength(), float):
            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
                raise ReferenceParticleNotSpecified()

            iPrev = len(iRefPrtcl.getPrOut()) - 1
            p0    = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                                    iRefPrtcl.getPrOut()[iPrev][:3]))

            Brho = (1./(speed_of_light*1.E-9))*p0/1000.
            Rgty = p0 / constants_instance.getparticleCHARGE(
                                                iRefPrtcl.getSpecies())
            if self.getDebug():
                print("     ----> Brho, Rgty:", Brho, Rgty)

            E0    = mth.sqrt(Rgty**2 + protonMASS**2)
            b02   = (Rgty/E0)**2
            g02   = 1./(1.-b02)
            g0    = mth.sqrt(g02)
            
            B0 = self.getStrength()
            # Bug: Will Shields; remove * 2.*Brho <---- 06Mar24
            if self.getDebug():
                print("     ----> epsilon0SI, B0, protonMASSSI, g0:", \
                      epsilon0SI, B0, protonMASSSI, g0)
            ne = epsilon0SI * B0**2 * mth.sqrt(b02) / (2.*Rgty) * \
                constants_instance.kginMeV()

            ne_trans = ne
            ne_longi = ne
        else:
            if self.getBz() == None or              \
               self.getVA() == None or              \
               self.getRA() == None or              \
               self.getRp() == None:
                raise badParameter( \
                        "BeamLineElement.GaborLens.setElectronDensity:" \
                        " no parameters!")
            
            ne_trans = (epsilon0SI * self.getBz()**2) / (2. * electronMASSSI)
            ne_longi = 4. * epsilon0SI * self.getVA() / \
                (electronCHARGESI * self.getRp()**2 * \
                 (1. + 2.*mth.log(self.getRA()/self.getRp())) \
                                                 )
        if self.getDebug():
            print("     <---- ne_trans, ne_longi:", ne_trans, ne_longi)

        self._ElectronDensity = min(ne_trans, ne_longi)
        
        if self.getDebug():
            print(" <---- Electron density:", self.getElectronDensity())
            
    def setTransferMatrix(self, _R):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()

        iPrev = len(iRefPrtcl.getPrOut()) - 1

        p0  = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iPrev][:3], \
                              iRefPrtcl.getPrOut()[iPrev][:3]))
        E0  = iRefPrtcl.getPrOut()[iPrev][3]
        b02 = (p0/E0)**2
        g02 = 1./(1.-b02)
        g0  = mth.sqrt(g02)
        
        if self.getDebug():
            print(" GaborLens(BeamLineElement).setTransferMatrix:")
            print("     ----> Reference particle 4-mmtm:", \
                  iRefPrtcl.getPrIn()[0])
            b0 = mth.sqrt(b02)
            print("         ----> p0, E0:", p0, E0)
            print("     <---- b0, g0:", b0, g0)
            
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Trace space:", _R)

        E = E0 + p0*_R[5]
        particleMASS = PhysCnst.PhysicalConstants().getparticleMASS(   \
                                               iRefPrtcl.getSpecies() )
        p = mth.sqrt(E**2 - particleMASS**2)
        g = E / particleMASS  
        
        if self.getDebug():
            print("     ----> Particle energy and mmtm:", E, p)

        l      = self.getLength()
        ne     = self.getElectronDensity()
        if self.getDebug():
            print("     ----> l, ne:", l, ne)

        if self.getDebug():
            print("     ----> alpha, electricCHARGE, epsilon0:", \
                  alpha, electricCHARGE, epsilon0)
            print("     ----> Joule2MeV:", Joule2MeV)
            print("     ----> m2InvMeV:", m2InvMeV)

        k      = (electricCHARGE**2 * particleMASS * g) / \
                 (2.*epsilon0 * p**2) * \
                 ne /m2InvMeV * \
            PhysCnst.PhysicalConstants().getparticleCHARGE(   \
                                                   iRefPrtcl.getSpecies() )
        w      = mth.sqrt(k)
        if self.getDebug():
            print("     ----> k, w:", k, w)
        
        cwl  = mth.cos(w*l)
        swl  = mth.sin(w*l)

        TrnsMtrx = np.array([                               \
            [    cwl, swl/w,     0.,    0.,        0., 0.], \
            [ -w*swl,   cwl,     0.,    0.,        0., 0.], \
            [     0.,    0.,    cwl, swl/w,        0., 0.], \
            [     0.,    0., -w*swl,   cwl,        0., 0.], \
            [     0.,    0.,     0.,    0.,  1., l/b02/g02],\
            [     0.,    0.,     0.,    0.,  0.,        1.] \
                             ])
        
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" GaborLens(BeamLineElement).visualise: start")
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        xlim = axs.get_xlim()
        ylim = axs.get_ylim()

        wdth = self.getLength()
        hght = 0.4
        angl = 0.
        abt  = 'center'
            
        if CoordSys == "RPLC":
            if self.getDebug():
                print("     ----> RPLC:")
                
            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)

            hght = (ylim[1] - ylim[0]) / 10.
            sxy   = [ sStrt, -hght/2. ]
            
        elif CoordSys == "Lab":
            if self.getDebug():
                print("     ----> Lab:")

            iCoord = 0
            if Proj == "yz": iCoord = 1

            cntr = self.getrStrt()
            if self.getDebug():
                print("         ----> cntr:", cntr)
            
            hght = (ylim[1] - ylim[0]) / 10.
            sxy   = [ cntr[2], cntr[iCoord]-hght/2. ]        

        if self.getDebug():
            print("     ----> Centre:", cntr)
            print("     ---->    sxy:", sxy)
            print("     ---->   wdth:", wdth)
            print("     ---->   hght:", hght)
            print("     ---->   angl:", angl)
            print("     ---->    abt:", abt)
            
        Patch = patches.Rectangle(sxy, wdth, hght, \
                                  angle=angl, \
                                  rotation_point=abt, \
                                  facecolor=('orange', 1.), \
                                  zorder=2)

        axs.add_patch(Patch)
                                   
        if self.getDebug():
            print(" <---- GaborLens(BeamLineElement).visualise: ends.")

            
# -------- Get methods:
#..   Methods believed to be self-documenting(!)
    def getBz(self):
        return self._Bz

    def getVA(self):
        return self._VA

    def getRA(self):
        return self._RA

    def getRp(self):
        return self._Rp

    def getLength(self):
        return self._Length

    def getStrength(self):
        return self._Strength

    def getElectronDensity(self):
        return self._ElectronDensity
    

#--------  I/o methods:
    def writeElement(self, dataFILE):
        if self.getDebug():
            print( \
                " GaborLens(BeamLineElement).writeElement starts.")

        derivedCLASS = "GaborLens"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Derived class:", bversion.decode('utf-8'))

        if self.getBz() == None:
            Bz = -1.
        else:
            Bz = self.getBz()
            
        if self.getVA() == None:
            VA = -1.
        else:
            VA = self.getVA()

        if self.getRA() == None:
            RA = -1.
        else:
            RA = self.getRA()

        if self.getRp() == None:
            Rp = -1.
        else:
            Rp = self.getRp()

        if self.getLength() == None:
            Ln = -1.
        else:
            Ln = self.getLength()

        if self.getStrength() == None:
            St = -1.
        else:
            St = self.getStrength()
        
        record = strct.pack(">6d", Bz, VA, RA, Rp, Ln, St)

        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Bz, VA, RA, Rp, Ln, St:", \
                  strct.unpack(">6d",record))

        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print( \
             " <---- GaborLens(BeamLineElement).writeElement done.")

        return
    
    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" GaborLens(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read((6*8))
        if brecord == b'':
            return True
        
        record = strct.unpack(">6d", brecord)

        Bz     = None
        VA     = None
        RA     = None
        Rp     = None
        Ln     = None
        St     = None

        if float(record[0]) != -1.:
            Bz     = float(record[0])
        if float(record[1]) != -1.:
            VA     = float(record[1])
        if float(record[2]) != -1.:
            RA     = float(record[2])
        if float(record[3]) != -1.:
            Rp     = float(record[3])
        if float(record[4]) != -1.:
            Ln     = float(record[4])
        if float(record[5]) != -1.:
            St     = float(record[5])

        if cls.getDebug():
            print("     ----> Bz, VA, RA, Rp, Ln, St:", \
                  Bz, VA, RA, Rp, Ln, St)
            
        return EoF, Bz, VA, RA, Rp, Ln, St


"""
Derived class CylindricalRFCavity:
==================================

** Note; KL; 05Mar24: Linac convention is phase relative to crest.

  CylindricalRFCavity class derived from BeamLineElement to contain
  paramters for a cylindrical RF cavity operated in the TM(010) mode.


  Class attributes:
  -----------------
    instances : List of instances of CylindricalRFCavity(BeamLineElement)
                class
  __Debug     : Debug flag


  Parent class attributes:
  ------------------------
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.
  _TrnsMtrx : Transfer matrix.


  Instance attributes to define drift:
  ------------------------------------
  _       : 
  
    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

  Set methods:
 setTransferMatrix : "Calculate" and set transfer matrix.

  Get methods:


"""
class CylindricalRFCavity(BeamLineElement):
    instances = []
    __Debug = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Gradient=None, _Frequency=None, _Phase=None):
        if self.__Debug:
            print(" CylindricalRFCavity.__init__: ", \
                  "creating the CylindricalRFCavity object:")
            print("     ----> Gradient (MV/m):", _Gradient)
            print("     ----> Frequency (MHz):", _Frequency)
            print("     ----> Phase     (rad):", _Phase)

        CylindricalRFCavity.instances.append(self)

        OK = self.setAll2None()

        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, \
                                 _Name, _rStrt, _vStrt, _drStrt, _dvStrt)


        if not isinstance(_Gradient, float) or \
           not isinstance(_Frequency, float) or \
           not isinstance(_Phase, float):
            raise badBeamLineElement( \
                                      "CylindricalRFCavity:" + \
                                      " bad specification:"  + \
                                      " gradient, frequency, phase:" + \
                                      _Gradient, _Frequency, _Phase)


        self.setGradient(_Gradient)
        self.setFrequency(_Frequency)
        self.setPhase(_Phase)

        _AngularFrequency  = self.getFrequency()*2.*mth.pi * 10.**6
        self.setAngularFrequency(_AngularFrequency)
        _WaveNumber        = self.getAngularFrequency() / speed_of_light
        self.setWaveNumber(_WaveNumber)
        
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()
        iPrev = len(iRefPrtcl.getPrOut()) - 1
        b0        = iRefPrtcl.getb0(iPrev)
        g0b0      = iRefPrtcl.getg0b0(iPrev)
        
        _Length   = mth.pi*b0*speed_of_light / \
            self.getAngularFrequency()
        self.setLength(_Length)

        _Radius   = sp.special.jn_zeros(0, 1)[0]/self.getWaveNumber()
        self.setRadius(_Radius)

        """
          Eqn 3.141 in Wolsli seems to have an error, expression for
          transit time factor below rederived 16Feb24, need to check
          with Andy.
        """
        _TransitTimeFactor = (2.*b0)                                 / \
                             (self.getWaveNumber()*self.getLength()) * \
                  mth.sin(self.getWaveNumber()*self.getLength()/2./b0)
        self.setTransitTimeFactor(_TransitTimeFactor)
        _V0 = self.getLength()*self.getGradient()*self.getTransitTimeFactor()
        self.setV0(_V0)

        _alpha = self.getV0()/iRefPrtcl.getMomentumIn(iPrev)/1000. * \
            PhysCnst.PhysicalConstants().getparticleCHARGE(   \
                                                   iRefPrtcl.getSpecies() )
        self.setalpha(_alpha)

        _wperp = self.getWaveNumber()*mth.sqrt( \
                        self.getalpha()*mth.cos(self.getPhase())/2./mth.pi)
        self.setwperp(_wperp)
        _cperp = mth.cos(self.getwperp()*self.getLength())
        self.setcperp(_cperp)
        _sperp = mth.sin(self.getwperp()*self.getLength()) / self.getwperp()
        self.setsperp(_sperp)
        
        _wprll = self.getWaveNumber()*mth.sqrt( \
                        self.getalpha()*mth.cos(self.getPhase())/mth.pi) / \
                        g0b0
        self.setwprll(_wprll)
        _cprll = mth.cos(self.getwprll()*self.getLength())
        self.setcprll(_cprll)
        _sprll = mth.sin(self.getwprll()*self.getLength()) / self.getwprll()
        self.setsprll(_sprll)

        self.setTransferMatrix()
        self.setmrf()

        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        if self.__Debug:
            print("     ----> New CylindricalRFCavity instance: \n", self)

    def __repr__(self):
        return "CylindricalRFCavity()"

    def __str__(self):
        print(" CylindricalRFCavity:")
        print(" --------------------")
        print("     ---->        Debug flag:", CylindricalRFCavity.getDebug())
        print("     ---->   Gradient (MV/m):", self.getGradient())
        print("     ---->   Frequency (MHz):", self.getFrequency())
        print("     ---->       Phase (rad):", self.getPhase())
        print("     ----> Derived quantities:")
        print("         ---->  Angular frequency (rad/s):", \
              self.getAngularFrequency())
        print("         ---->           Wave number (/m):", \
              self.getWaveNumber())
        print("         ---->                 Length (m):", \
              self.getLength())
        print("         ---->                 Radius (m):", \
              self.getRadius())
        print("         ---->          TransitTimeFactor:", \
              self.getTransitTimeFactor())
        print("         ---->                    V0 (MV):", \
              self.getV0())
        print("         ---->                      alpha:", \
              self.getalpha())
        print("         ---->                      wperp:", \
              self.getwperp())
        print("         ---->                      cperp:", \
              self.getcperp())
        print("         ---->                      sperp:", \
              self.getsperp())
        print("         ---->                      wprll:", \
              self.getwprll())
        print("         ---->                      cprll:", \
              self.getcprll())
        print("         ---->                      sprll:", \
              self.getsprll())
        print("     <---- End derived quantities:")
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> m_rf: \n", self.getmrf())
        BeamLineElement.__str__(self)
        return " <---- CylindricalRFCavity parameter dump complete."

    def SummaryStr(self):
        Str  = "CylindricalCavity: " + BeamLineElement.SummaryStr(self) + \
            "; Gradient = " + str(self.getGradient())   + \
            "; Frequency = " + str(self.getFrequency()) + \
            "; Phase = " + str(self.getPhase())
        return Str


# -------- "Set methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def setDebug(cls, _Debug=False):
        if not isinstance(_Debug, bool):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setDebug:" + \
                                " bad flag:", _Debug)
        cls.__Debug = _Debug
        
    def setAll2None(self):
        self._Gradient          = None
        self._Frequency         = None
        self._Phase             = None
        self._TransitTimeFactor = None
        self._V0                = None
        self._alpha             = None
        self._wperp             = None
        self._cperp             = None
        self._sperp             = None
        self._wprll             = None
        self._cprll             = None
        self._sprll             = None

        self._TrnsMtrx  = None

    def setGradient(self, _Gradient):
        if not isinstance(_Gradient, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setVoltage:" + \
                    " bad gradient:", _Gradient)
        self._Gradient = _Gradient

    def setFrequency(self, _Frequency):
        if not isinstance(_Frequency, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setFrequency:" + \
                                " bad frequency:", _Frequency)
        self._Frequency = _Frequency

    def setAngularFrequency(self, _AngularFrequency):
        if not isinstance(_AngularFrequency, float):
            raise badParameter( \
             "BeamLineElement.CylindricalRFCavity.setAngularFrequency:" + \
                      " bad angular frequency:", _AngularFrequency)
        self._AngularFrequency = _AngularFrequency

    def setPhase(self, _Phase):
        if not isinstance(_Phase, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _Phase)
        self._Phase = _Phase

    def setWaveNumber(self, _WaveNumber):
        if not isinstance(_WaveNumber, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _WaveNumber)
        self._WaveNumber = _WaveNumber

    def setLength(self, _Length):
        if not isinstance(_Length, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _Length)
        self._Length = _Length        

    def setRadius(self, _Radius):
        if not isinstance(_Radius, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _Radius)
        self._Radius = _Radius        

    def setTransitTimeFactor(self, _TransitTimeFactor):
        if not isinstance(_TransitTimeFactor, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _TransitTimeFactor)
        self._TransitTimeFactor = _TransitTimeFactor        
        
    def setV0(self, _V0):
        if not isinstance(_V0, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _V0)
        self._V0 = _V0        
        
    def setalpha(self, _alpha):
        if not isinstance(_alpha, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _alpha)
        self._alpha = _alpha        
        
    def setwperp(self, _wperp):
        if not isinstance(_wperp, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _wperp)
        self._wperp = _wperp        
        
    def setcperp(self, _cperp):
        if not isinstance(_cperp, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _cperp)
        self._cperp = _cperp        
        
    def setsperp(self, _sperp):
        if not isinstance(_sperp, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _sperp)
        self._sperp = _sperp        
        
    def setwprll(self, _wprll):
        if not isinstance(_wprll, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _wprll)
        self._wprll = _wprll        
        
    def setcprll(self, _cprll):
        if not isinstance(_cprll, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _cprll)
        self._cprll = _cprll        
        
    def setsprll(self, _sprll):
        if not isinstance(_sprll, float):
            raise badParameter( \
                    "BeamLineElement.CylindricalRFCavity.setPhase:" + \
                                " bad phase:", _sprll)
        self._sprll = _sprll        
        
        
    def setmrf(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()
        iPrev  = len(iRefPrtcl.getPrOut()) - 1
        g02b02 = iRefPrtcl.getg0b0(iPrev)**2

        if self.getDebug():
            print(" CylindricalRFCavity(BeamLineElement).setmrf:")
            print("     ----> Reference particle g02b02:", g02b02)
        
        _mrf   = np.array([0., 0., 0., 0., \
        (1.-self.getcprll())*mth.tan(self.getPhase())/self.getWaveNumber(), \
        g02b02*self.getwprll()**2*self.getsprll()*                          \
                        mth.tan(self.getPhase())/self.getWaveNumber()])
            
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> mrf:", _mrf)

        self._mrf = _mrf

    def setTransferMatrix(self):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()
        iPrev  = len(iRefPrtcl.getPrOut()) - 1
        g02b02 = iRefPrtcl.getg0b0(iPrev)**2

        if self.getDebug():
            print(" CylindricalRFCavity(BeamLineElement).setTransferMatrix:")
            print("     ----> Reference particle g02b02:", g02b02)
        
        TrnsMtrx = np.array([
            [                    self.getcperp(), self.getsperp(), \
                                                         0., 0., 0., 0.], \
            [-self.getwperp()**2*self.getsperp(), self.getcperp(), \
                                                         0., 0., 0., 0.], \
            [0., 0.,                     self.getcperp(), self.getsperp(), \
                                                                 0., 0.], \
            [0., 0., -self.getwperp()**2*self.getsperp(), self.getcperp(), \
                                                                 0., 0.], \
            [0., 0., 0., 0., \
                 self.getcprll(),                 self.getsprll()/g02b02], \
            [0., 0., 0., 0., \
                 -g02b02*self.getwprll()**2*self.getsprll(), \
                                                         self.getcprll()] \
        ])

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" CylindricalRFCavity(BeamLineElement).visualise: start")
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        xlim = axs.get_xlim()
        ylim = axs.get_ylim()

        wdth = self.getLength()
        hght = 2.*self.getRadius()
        angl = 0.
        abt  = 'center'

        hght = min(2.*self.getRadius(), (ylim[1] - ylim[0]) / 7.5)
            
        if CoordSys == "RPLC":
            if self.getDebug():
                print("     ----> RPLC:")
                
            iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
            iAddr     = iRefPrtcl.getLocation().index(self.getName())
            sStrt = iRefPrtcl.gets()[iAddr-1]
            if self.getDebug():
                print("         ----> self.getName(), iAddr, sStrt:", \
                      self.getName(), iAddr, sStrt)

            sxy   = [ sStrt, -hght/2. ]
            
        elif CoordSys == "Lab":
            if self.getDebug():
                print("     ----> Lab:")
                
            cntr = self.getrStrt()
            if self.getDebug():
                print("         ----> cntr:", cntr)
                
            if Proj == "xz":
                sxy   = [ cntr[2], cntr[0]-hght/2. ]
            elif Proj == "yz":
                sxy   = [ cntr[2], cntr[1]-hght/2. ]
        if self.getDebug():
            print("     ----> Centre:", cntr)
            print("     ---->    sxy:", sxy)
            print("     ---->   wdth:", wdth)
            print("     ---->   hght:", hght)
            print("     ---->   angl:", angl)
            print("     ---->    abt:", abt)
            
        Patch = patches.Rectangle(sxy, wdth, hght, \
                                  angle=angl, \
                                  rotation_point=abt, \
                                  facecolor=('slategray', 1.), \
                                  zorder=2)

        axs.add_patch(Patch)
                                   
        if self.getDebug():
            print( \
                " <---- CylindricalRFCavity(BeamLineElement).visualise: ends.")

        
# -------- "Get methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug
    
    def getGradient(self):
        return self._Gradient

    def getFrequency(self):
        return self._Frequency

    def getAngularFrequency(self):
        return self._AngularFrequency

    def getPhase(self):
        return self._Phase

    def getWaveNumber(self):
        return self._WaveNumber

    def getLength(self):
        return self._Length

    def getRadius(self):
        return self._Radius

    def getTransitTimeFactor(self):
        return self._TransitTimeFactor

    def getV0(self):
        return self._V0

    def getalpha(self):
        return self._alpha

    def getwperp(self):
        return self._wperp

    def getcperp(self):
        return self._cperp

    def getsperp(self):
        return self._sperp

    def getwprll(self):
        return self._wprll

    def getcprll(self):
        return self._cprll

    def getsprll(self):
        return self._sprll

    def getmrf(self):
        return self._mrf

    
#--------  Utilities:
    def Transport(self, _R=None):
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        if not isinstance(iRefPrtcl, Prtcl.ReferenceParticle):
            raise ReferenceParticleNotSpecified()
        iPrev  = len(iRefPrtcl.getPrOut()) - 1
        g02b02 = iRefPrtcl.getg0b0(iPrev)**2
        
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
            " CylindricalRFCavity(BeamLineElement).Transport:" + \
                                "bad input vector:", \
                                _R)

        if self.getDebug():
            print(" CylindricalRFCavity(BeamLineElement).Transport:", \
                  "     ---->", self.SummaryStr())
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> _R:", _R)
            print("     ----> Outside:", self.OutsideBeamPipe(_R))
            
        if self.OutsideBeamPipe(_R) or \
           self.ExpansionParameterFail(_R) or \
           abs(_R[4]) > 5.:
            _Rprime = None
        else:
            detTrnsfrMtrx = np.linalg.det(self.getTransferMatrix())
            error         = abs(1. - abs(detTrnsfrMtrx))
            if error > 1.E-6:
                print(" CylindricalRFCavity(BeamLineElement).Transport:(3):" \
                      " detTrnsfrMtrx:", detTrnsfrMtrx)
            
            _Rprime = self.getTransferMatrix().dot(_R) + self.getmrf()

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" <---- Rprime:", _Rprime)

        if not isinstance(_Rprime, np.ndarray):
            pass

        return _Rprime

#--------  I/o methods:
    def writeElement(self, dataFILE):
        if self.getDebug():
            print( \
                " CylindricalRFCavity(BeamLineElement).writeElement starts.")

        derivedCLASS = "CylindricalRFCavity"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Derived class:", bversion.decode('utf-8'))

        record = strct.pack(">3d", \
                            self.getGradient(), \
                            self.getFrequency(), \
                            self.getPhase())

        dataFILE.write(record)
        if self.getDebug():
            print("         ----> Gradient, frequency, phase:", \
                  strct.unpack(">3d",record))

        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print( \
             " <---- CylindricalRFCavity(BeamLineElement).writeElement done.")

        return

    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" Source(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read((3*8))
        if brecord == b'':
            return True
        
        record  = strct.unpack(">3d", brecord)
        Grdnt   = float(record[0])
        Frqncy  = float(record[1])
        Phs     = float(record[2])

        if cls.getDebug():
            print("     ----> Grdnt, Frqncy, Phs:", Grdnt, Frqncy, Phs)
            
        return EoF, Grdnt, Frqncy, Phs
        
    
"""
Derived class Source:
=====================

  Source class derived from BeamLineElement to contain paramters for the
  source.  __init__ sets parameters of the source.  


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.


  Instance attributes to define source:
  -------------------------------------
  _Mode  : Int, Mode
  _Param : List of parameters,
           Parameterised laser-driven source (Mode=0):
             [0] - Wavelength - microns
             [1] - Power - W
             [2] - Strehl ratio
             [3] - Focal spot radius - m
             [4] - Laser-pulse duration - s
             [5] - Hot electron temperature - MeV
             [6] - Minimum proton kinetic energy - MeV
             [7] - Maximum proton kinetic energy - MeV
             [8] - Intercept of sigma_{theta_S} at K=0 [degrees]
             [9] - Scaled slope of sigma_{theta_S}     [degrees]
            [10] - rp max

           Gaussian (Mode=1):
             [0] - Sigma of x gaussian - m
             [1] - Sigma of y gaussian - m
             [2] - Minimum cos theta to generate
             [3] - Kinetic energy - MeV
             [4] - Sigma of gaussian - MeV

           Flat (Mode=2)

           Read from file (Mode=3)
             [0] = User should set to BeamIO instance of file to be read

           Uniform disc (Mode=4):
             [0] - Kinetic energy - MeV
             [1] - Sigma of gaussian - MeV
             [2] - Radius of disc

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising quad
                parameters.

  Set methods:
        setAll2None : Set all instance attributes to None
            setMode : Set mode (int)
        setModeText : Set mode text--string identifying source
      setParameters : Source parameters, fitting with Mode

  Get methods:
        getAll2None : Set all instance attributes to None
            getMode : Set mode (int)
        getModeText : Set mode text--string identifying source

  Utilities and processing methods:
     cleaninstances : Deletes source instance(s) and empties list of sources.

    CheckSoureParam : Checks source parameters are self consistent.  Calls
                      CheckMode and CheckParam
                  Input : Mode (Int), Parameters (list, as defined for Mode)
                 Return : True/False

          CheckMode : Checks mode
                  Input : Mode (Int)
                 Return : True/False

         CheckParam : Check parameters are consistent with Mode
                  Input : Mode (Int), Parameters (list, as defined for Mode)
                 Return : True/False

getParticleFromSource : Generate one particle at source.  Check that particle
                        meets cos(theta) requirement programmed in source
                        parameters.
                  Input : None
                 Return : np.ndarray : 6D phase space of particle at source.

        getParticle : Generate one particle at source using Mode.
                  Input : None
                 Return : np.ndarray : 6D phase space of particle at source.

    getFlatThetaPhi : Generate direction of particle at source flat in
                      cos(theta) and phi
                  Input : None
                 Return : cos(theta) [float], phi [Float]

getLaserDrivenProtonEnergy: Generate proton energy at source using
                            parameterised TNSA spectrum.
                  Input : None
                 Return : Energy [float]

      getTraceSpace : Convert x, y, energy, cos(theta), phi [input] to
                      trace space.

                  Input : x, y, energy, cos(theta), phi [floats]
                 Return : np.ndarray : 6D phase space of particle at source.


"""

class Source(BeamLineElement):
    instances  = []
    __Debug    = False

    ModeList   = [0, 1, 2, 3, 4]
    ModeText   = ["Parameterised laser driven", "Gaussian", "Flat", \
                  "Read from file", "UniformDisc"]
    
    ParamUnit  = [ ["$\\mu$m", "W", " ", "m", "s", "MeV", \
                    "MeV", "MeV", \
                    "degrees", "degrees", " "], \
                   ["m", "m", "MeV", "MeV", ""], \
                   [" "], [], \
                   ["MeV", "MeV", "m"] ]
    
    ParamText  = [ ["Wavelength", "Power", "Strehl ratio", "r0", \
                    "Duration", "Te", \
                    "Kmin", "Kmax", \
                    "SigmaThetaS0", "SlopeThetaS", "rpmax"], \
                   ["SigmaX", "SigmaY", "MinCTheta", \
                    "MeanEnergy", "SigmaEnergy"],    \
                   [], \
                   ["imput file (BeamIO instance)"], \
                   ["MeanEnergy", "SigmaEnergy", "Radius"] ]
    
    ParamList  = [ [float, float, float, float, float, float, \
                    float, float, float, float, float], \
                   [float, float, float, float, float],          \
                   [float, float, float, float, float],          \
                   [], \
                   [float, float, float] ]
    ParamLaTeX  = [ ["Wavelength", "Power", "Strehl ratio", \
                     "Focal spot radius",         \
                     "Electron Temperature", "Duration", \
                     "$K_{\\rm min}$", "$K_{\\rm max}$", \
                     "Intercept of $\\sigma_{\\theta_S}$",               \
                     "Scaled slope of $\\sigma_{\\theta_S}$", "rpmax"],  \
                    ["\\sigma_x", "\\sigma_y",                           \
                     "\\cos\\theta_S |_{\\rm min}",                      \
                     "Mean kinetic energy",                              \
                     "Kinetic energy standard deviation"],               \
                    [], [], \
                    ["Mean kinetic energy",                              \
                     "Kinetic energy standard deviation",                \
                     "Maximum radius"] ]

    Lsrdrvng_E = None
    LsrDrvnIni = False

#--------  Initialisation and built-in methods  --------  --------  --------
    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _Mode=None, _Param=[]):

        Source.LsrDrvnIni = False
        
        if self.getDebug():
            print(' Source.__init__: ', \
                  'creating the Source object')
            print("     ----> Parameters:", _Param)

        Source.instances.append(self)

        self.setAll2None()
        
        #.. BeamLineElement class initialisation:
        BeamLineElement.__init__(self, _Name, _rStrt, _vStrt, \
                                              _drStrt, _dvStrt)
        
        self.setLength(0.)
        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())
        
        #.. Check valid mode and parameters:
        if _Mode==None or _Param==None:
            print(" BeamLineElement(Source).__init__:", \
                  " bad source paramters:", \
                      _Mode, _Param)
            raise badBeamLineElement( \
                    " Source: bad specification for source!")
        
        if self.getDebug():
            print("     ----> Before parameter checks; Mode, Param:", \
                  _Mode, _Param)

        params = _Param
        try:
            ValidSourceParam = self.CheckSourceParam(_Mode, _Param)
        except:
            if len(_Param) == 16:
                params = self.paramsFROMlegacyPARAMS(_Param)

        ValidSourceParam = self.CheckSourceParam(_Mode, params)

        if not ValidSourceParam:
            print(" BeamLineElement(Source).__init__:", \
                  " bad source input; :", \
                  "_Mode, _Param:", _Mode, _Param)
            raise badSourceSpecification( \
                        " BeamLineElement(Source).__init__:", \
                        " bad source paramters. Exit", \
                        _Name
                                                            )
        self.setMode(_Mode)
        self.setModeText(Source.ModeText[_Mode])
        self.setParameterText(Source.ParamText[_Mode])
        self.setParameters(params)
        self.setParameterUnit(Source.ParamUnit[_Mode])
                
        if self.__Debug:
            print("     ----> New Source instance: \n", \
                  self)
            print(" <---- BeamLineElement(Source) instance created.")
            
    def __repr__(self):
        return "Source()"

    def __str__(self):
        print(" BeamLineElement(Source):")
        print(" ------------------------")
        print("     ----> Debug flag:", Source.getDebug())
        print("     ----> Mode      :", self.getMode())
        print("     ----> Mode text :", self.getModeText())
        print("     ----> Parameters:", self.getParameters())
        return " <---- Source parameter dump complete."

    def SummaryStr(self):
        Str  = "Source           : " + BeamLineElement.SummaryStr(self) + \
            "; Mode = " + str(self.getMode()) + "; paramters = " + \
            str(self.getParameters())
        return Str

#--------  "Set methods"  --------  --------  --------  --------  --------
#.. Methods believed to be self documenting(!)
    def setAll2None(self):
        self._Mode         = None
        self._ModeText     = None
        self._Param        = None
        self._derivedParam = []
        
    def setMode(self, _Mode):
        if self.getDebug():
            print(" Source.setMode; Mode:", _Mode)
        self._Mode = _Mode

    def setModeText(self, _ModeText):
        if self.getDebug():
            print(" Source.setParamters; Mode:", _ModeText)
        self._ModeText = _ModeText

    def setParameterText(self, _ParameterText):
        if self.getDebug():
            print(" Source.setParamters; Parameter:", _ParameterText)
        self._ParameterText = _ParameterText

    def setParameters(self, _Param):
        if self.getDebug():
            print(" Source.setParamters; Parameters:", _Param)
        self._Params = _Param
         
    def setParameterUnit(self, _ParameterUnit):
        if self.getDebug():
            print(" Source.setParamters; Parameter:", _ParameterUnit)
        self._ParameterUnit = _ParameterUnit

        
#--------  "get methods"  --------  --------  --------  --------  --------
#.. Methods believed to be self documenting(!)
    def getderivedParameters(self):
        """
            [0] - Gamma - Normalisation constant for cumulative
                          probability
        """
        
        return self._derivedParam

    
#--------  Processing methods:
    def getParticleFromSource(self):
        if self.getDebug():
            print(" BeamLineElement(Source).getParticleFromSource: start")

        #.. Generate initial particle:
        x, y, K, cTheta, Phi, xp, yp = self.getParticle()
        if self.__Debug:
            print("     ----> x, y, K, cTheta, Phi:", \
                  x, y, K, cTheta, Phi)

        #.. Convert to trace space:
        TrcSpc = self.getTraceSpace(x, y, K, cTheta, Phi, xp, yp)
        if self.__Debug:
            print("     ----> Trace space:", TrcSpc)

        if not isinstance(TrcSpc, np.ndarray):
            if self.__Debug:
                raise FailToCreateTraceSpaceAtSource()
                
        if self.getDebug():
            print(" <---- BeamLineElement(Source).getParticleFromSource,", \
                  " done.", \
                  '  --------  --------  --------  --------  --------')

        return TrcSpc

#--------  "Management" methods  --------  --------  --------  --------
    #.. clean all source instances:
    @classmethod
    def cleaninstances(cls):
        if cls.getDebug():
            print(" Source(BeamLineElement).cleaninstance:")
        for inst in cls.getinstances():
            if cls.getDebug():
                print("     ----> Kill:", inst.getName())
                
            iAddrBLE = BeamLineElement.getinstances().index(inst)
            BeamLineElement.getinstances().pop(iAddrBLE)

            del inst

        cls.instances = []
        if cls.getDebug():
            print(' <---- Instances removed.')

    @classmethod
    #.. Check through data frame and report legacy lines:
    def scanLEGACY(cls, pndsDF):
        if cls.getDebug():
            print(" Source(BeamLineElement).scanLEGACY starts.", \
                  "Source data frame:")
            print(pndsDF)
            
        cleanedPNDS = pndsDF
        
        legacyPARAMs = ["SigmaX", "SigmaY", "Emax", "nPnts", "MinCTheta", \
                        "Energy", "Intensity"]
        
        cleanedPNDS = pndsDF[~pndsDF['Parameter'].isin(legacyPARAMs)]

        #.. Minimum kinetic energy:
        if bool(cleanedPNDS[cleanedPNDS["Parameter"] == \
                           "Emin"].any().any()):
            print(" Source(BeamLineElement).scanLEGACY:", \
                  "Depricated parameter Emin,", \
                  "use Kmin.")

            cleanedPNDS.loc[ \
                cleanedPNDS["Parameter"]=="Emin", ["Parameter"]] = "Kmin"

        if cls.getDebug():
            print("     ----> Cleaned DF:")
            print(cleanedPNDS)

        if len(cleanedPNDS) > len(pndsDF):
            print(" Source(BeamLineElement).scanLEGACY:", \
              "The following legacy parameters will be ignored:")
            legacyPNDS   = pndsDF['Parameter'].isin(legacyPARAMs)
            legacyMASKED = pndsDF[legacyPNDS]
            print(legacyMASKED)
            print("     ----> Remaining parameters:")
            print(cleanedPNDS)
            print(" <---- BeamLineElement.Source.scanLEGACY: Done.")

        return cleanedPNDS

    @classmethod
    #.. Check through data frame and report legacy lines:
    def setDEFAULTparams(cls):
        if cls.getDebug():
            print(" BeamLineElement.Source.setDEFAULTparams starts.")
            print("     ----> Default to J-KAREN-P to compare with", \
                  "Dover et al,", \
                  "High Energy Density Physics 37 (2020) 100847.")

        wavelength   = 0.8                   #.. microns

        E_laser      = 10.                   #.. J
        Duration     = 40.E-15               #.. s
        power        = E_laser / Duration    #.. W
        strhlRATIO   = 3./5.

        r0           = 1.5E-6                #.. m
        Thickness    = 5.E-6                 #.. m
        
        Te           = 10.                   #.. MeV
        Kmin         = 1.                    #.. MeV
        Kmax         = None                  #.. MeV
        
        DivAngle     = 25.                   #.. degrees
        
        SigmaThetaS0 = 20.                   #.. degrees
        SlopeThetaS  = 15.                   #.. degrees
        rpmax        = -9999.

        if cls.getDebug():
            print("     ----> wavelength, power, strhlRATIO, r0, Duration,", \
                  "Te, Kmin, Kmax, Thickness, DivAngle, SigmaThetaS0,", \
                  "SlopeThetaS, rpmax:", wavelength, power, strhlRATIO, r0, \
                   Duration, Te, Kmin, Kmax, Thickness, DivAngle, \
                   SigmaThetaS0, SlopeThetaS, rpmax)
            print(" <---- Defaults set.")
        
        return wavelength, power, strhlRATIO, r0, Duration, Te, Kmin, Kmax, \
               Thickness, DivAngle, SigmaThetaS0, SlopeThetaS, rpmax

    @classmethod
    #.. Parse one parameter from Pandas data frame:
    def parseSINGLEparam(cls, pndsDF=None, Name=None, default=None):
        if cls.getDebug():
            print(" BeamLineElement.Source.parseSINGLEparam starts.", \
                  "id(pndsDF), Name, default:", id(pndsDF), Name, default)

        value = None

        if isinstance(pndsDF, pnds.core.frame.DataFrame) != None and \
           Name != None:
            if bool(pndsDF[pndsDF["Parameter"]==Name].any().any()):
                value = float( \
                 pndsDF[pndsDF["Parameter"]==Name]["Value"].iloc[0])
            else:
                if Name != "Te":
                    print("     ---->", Name, "is not defined;", \
                          "it will be set to", default)
                value = default

        else:
            raise badParameter(" BeamLineElement.Source.parseSINGLEparam:"+\
                              "id(pndsDF), Name, default:" + \
                              str(id(pndsDF)) + "," + \
                              str(Name) + "," + \
                              str(default))

        if cls.getDebug():
            print(" BeamLineElement.Source.parseSINGLEparam:", \
                  "Name=", Name, "; value:", value)
        
        return value

    @classmethod
    def calculateTe(cls, wavelength, r0, power, strhlRATIO):
        if cls.getDebug():
            print(" Source(BeamLineElement).calculateTe:", \
                  "wavelength, electronCHARGESI, electronMASS,", \
                  "electronMASSSI:", \
                  wavelength, electronCHARGESI, electronMASS, \
                  electronMASSSI)
            print("              ", \
                  "speed of light, r0, power, strhlRATIO, eps0 ", \
                  "                                  :", \
                  speed_of_light, r0, power, eps0SI       )

            I = power*strhlRATIO / (mth.pi * r0**2 * 10000.)
            print("              ", \
                  "I                                             ", \
                  "                              :", I)
            print(" ")

        a0 = (wavelength * electronCHARGESI) / \
            (2.*mth.pi*electronMASSSI*speed_of_light**2*r0) * \
            mth.sqrt(2.*power*strhlRATIO / (mth.pi*eps0SI*speed_of_light))

        if cls.getDebug():
            print("     ----> a0:", a0)
            print(" ")

        pt0 = a0 * electronMASSSI * speed_of_light
        y0  = a0 * wavelength / (2.*mth.pi)
                    
        if cls.getDebug():
            print("     ----> pt0, y0:", pt0, y0)
            print(" ")

        rLe = r0 / (2. * mth.log(2.))
                    
        if cls.getDebug():
            print("     ----> r0, rLe:", r0, rLe)
            print(" ")
                
        if rLe < y0:
            pt = pt0 * mth.sqrt(1. - (1. - rLe/y0)**2 )
        else:
            pt = pt0

        if cls.getDebug():
            print("     ----> pt:", pt)
            print(" ")

        Te = electronMASSSI * speed_of_light**2 * ( \
                mth.sqrt(1. + (pt/(electronMASSSI*speed_of_light))**2) \
                                                    -1. )

        Te /= electronCHARGESI * 1.E6
        
        if cls.getDebug():
            print(" <---- Te:", Te)

        return Te
    
    @classmethod
    def calculatet0(cls, power, strhlRATIO, r0, Thickness, DivAngle):
        t0 = None
        if cls.getDebug():
            print(" Source(BeamLineElement).calculatet0:", \
                  "power, strhlRATIO, r0, Thickness, DivAngle:", \
                   power, strhlRATIO, r0, Thickness, DivAngle)

        qi      = 1.                                #!! Fix h/wired number
        PR      = 8.71E9                            #!! Fix h/wired number
        theta_e = DivAngle * mth.pi / 180.
        intnsty = power*strhlRATIO / (mth.pi * r0**2) * 1E-4
        if cls.getDebug():
            print("     ----> qi, PR, theta_e, intnsty:", \
                  qi, PR, theta_e, intnsty)
        
        #.. Power conversion efficiency, eta:
        eta = 1.2E-15 * intnsty**(3./4.)            #!! Fix h/wired number
        if eta > 0.5: eta = 0.5                     #!! Fix h/wired number
        if cls.getDebug():
            print("     ----> eta:", eta)

        Kinfnty = 2.*electronMASSSI * speed_of_light**2 * \
            mth.sqrt(eta*power*strhlRATIO/PR)
        if cls.getDebug():
            print("     ----> Kinfnty:", Kinfnty)

        vmax = mth.sqrt(2.*Kinfnty/protonMASSSI)
        if cls.getDebug():
            print("     ----> vmax:", vmax)

        t0 = (r0 + Thickness*mth.tan(theta_e))/vmax

        if cls.getDebug():
            print(" <---- t0:", t0)

        return t0, Kinfnty

    def paramsFROMlegacyPARAMS(self, legacyPARAMS):
        r0 = mth.sqrt(legacyPARAMS[0]**2 + legacyPARAMS[1]**2)
        params = [ legacyPARAMS[8],  \
                   legacyPARAMS[6],  \
                   3./5.,            \
                   r0,               \
                   legacyPARAMS[4],  \
                   10.,              \
                   legacyPARAMS[3],  \
                   legacyPARAMS[4],  \
                   legacyPARAMS[13], \
                   legacyPARAMS[14], \
                   legacyPARAMS[15] ]

        print(" Source(BeamLineElement).paramsFROMlegacyPARAMS:")
        print("     ---->", \
              "Attempt to convert parameters to present structure.")
        print("          Legacy params:", legacyPARAMS)
        print(" <---- Converted params:", params)

        return params
        
    def getParticle(self):
        if self.getDebug():
            print(" BeamLineElement(Source).getParticle: start")
            print("     ----> Mode, parameters:", \
                  self.getMode(), self.getParameters())

        #.. Initialise o/p to None:
        X        = None
        Y        = None
        KE       = None
        cosTheta = None
        Phi      = None
        xp       = None
        yp       = None

        #-------- Laser driven:
        if self.getMode() == 0:
            #.. Proton kinetic energy:
            if self.getDebug():
                print("     ----> Calling getLaserDrivenProtonEnergy:")
                
            KE     = self.getLaserDrivenProtonEnergy()  # [MeV]
            if self.getDebug():
                print("     <---- KE:", KE)

            #.. position at production:
            X      = rnd.gauss(0., self.getParameters()[3])
            Y      = rnd.gauss(0., self.getParameters()[3])
            
            #.. x' and y' at production:
            upmax  = mth.sin(np.radians(self.g_theta(KE)))
            if self.getDebug():
                print("     ----> upmax:", upmax)

            Accept = False
            iCnt   = 0
            while not Accept:
                iCnt += 1
                if iCnt > 1E6:
                    raise KillInfiniteLoop(" iCnt: " + str(iCnt))
                
                xp    = rnd.uniform(-upmax, upmax)
                yp    = rnd.uniform(-upmax, upmax)
                if self.getDebug():
                    print("     ----> xp, yp:", xp, yp)
                
                grp   = self.getgofrp(upmax, xp, yp)
                if self.getDebug():
                    print("     ----> grp:", grp)

                Accept = False
                if rnd.random() < grp:
                    if self.getParameters()[10] == -9999.:
                        Accept = True
                    else:
                        rp = mth.sqrt(xp**2 + yp**2)
                        if self.getDebug():
                            print("     ----> rp, rpamx:", \
                                  rp, self.getParameters()[10])
                        if rp < self.getParameters()[10]:
                            Accept = True
                    
            if xp == yp:
                wrnngs.warn(" BeamLineElement.Source.getParticle:", \
                            " eqaual xp and yp")
            
            cosTheta, Phi = None, None # Backward compatibility!

        elif self.getMode() == 1:
            X             = rnd.gauss(0., self.getParameters()[0])
            Y             = rnd.gauss(0., self.getParameters()[1])
            cosTheta, Phi = self.getFlatThetaPhi()
            KE            = rnd.gauss(self.getParameters()[3], \
                                      self.getParameters()[4])
        elif self.getMode() == 2:
            X             = rnd.gauss(0., self.getParameters()[0])
            Y             = rnd.gauss(0., self.getParameters()[1])
            cosTheta, Phi = self.getFlatThetaPhi()
            KE            = rnd.uniform(self.getParameters()[3], \
                                        self.getParameters()[4])
        elif self.getMode() == 3:
            self.getParameters()[0].readBeamDataRecord()
            x = 1.
            y = 0.
            z = x/y
            exit()
        elif self.getMode() == 4:
            KE            = rnd.gauss(self.getParameters()[0], \
                                      self.getParameters()[1])
            rd            = self.getParameters()[2] * mth.sqrt(rnd.random())
            phi           = 2. * mth.pi * rnd.random()
            X             = rd * mth.cos(phi)
            Y             = rd * mth.sin(phi)
            xp            = 0.
            yp            = 0.
        else:
            raise badSOURCEmode(self.getMode())

        if self.getDebug():
            print("     ----> X, Y, KE, cosTheta, Phi, xp, yp:", \
                  X, Y, KE, cosTheta, Phi, xp, yp)
            
        if self.getDebug():
            print(" <---- BeamLineElement(Source).getParticle, done.", \
                  '  --------  --------  --------  --------  --------')

        return X, Y, KE, cosTheta, Phi, xp, yp

    #..  Used for Modes 1 and 2:
    def getFlatThetaPhi(self):
        cosTheta = rnd.uniform(self.getParameters()[2], 1.)
        Phi      = rnd.uniform( 0., 2.*mth.pi)
        return cosTheta, Phi

    #..  Calculate cumulative probability for parabolic distribution
    #    with max extent upmax:
    #..  Used for Mode 0
    def getgofrp(self, upmax, xp, yp):
        Nrm = 1./upmax**2
        rp2 = xp**2 + yp**2 
        grp = Nrm * (upmax**2 - rp2)
        
        return grp
    
    #..  Gaussian Angular Distribution with RMS determined by linearly
    #    decreasing paramterisation presented in documentation:
    #..  Used for Mode 0
    def g_theta(self, energy):
        Kmax  = self.getParameters()[7]
        if self.getDebug():
            print(" BeamLineElement(Source).g_theta:")
            print("     ----> enegy, Kmax:", energy, Kmax)
            print("     ----> theta_S(0), theta_S(slope):", \
                  self.getParameters()[8], self.getParameters()[9])
            
        theta = self.getParameters()[8] - \
            self.getParameters()[9] * energy / Kmax
        
        if self.getDebug():
            print(" <---- theta:", theta)

        return theta
    
    # Defines the function to solve for f(x) = 0
    #..  Used for Mode 0
    @classmethod
    def DurationBYt0equation(self, X, t_laser, t_0):
        if self.getDebug():
            print(" Source(BeamLineElement).DurationBYt0equation start:")
            print("     ----> X, t_laser, t_0:", 
                  X, t_laser, t_0)
        if   X >  0.9999: X = 0.9999
        elif X < -0.9999: X = -0.9999

        if self.getDebug():
            print("     ----> X:", X)

        eps= (X * (1 + (0.5 / (1 - (X**2) ) ) ) ) + \
                 (0.25 * np.log((1 + X) / (1 - X))) - (t_laser/t_0)

        if self.getDebug():
            print(" <---- eps:", eps)

        return eps

    # Generates energy values for the distribution
    #..  Used for Mode 0
    def getLaserDrivenProtonEnergy(self):
        if not Source.LsrDrvnIni:
            if self.__Debug:
                print( \
                    " BeamLineElement(Source).getLaserDrivenProtonEnergy:", \
                    " initialise")
                
            Source.LsrDrvnIni = True

            self.getLaserCumProbParam()

        Te   = self.getParameters()[5]
        Kmin = self.getParameters()[6]
        Kmax = self.getParameters()[7]

        Gamma = self.getderivedParameters()[0]

        if self.getDebug():
            print("     ----> Get kinetic energy:")
            
        GE = rnd.random()

        if self.getDebug():
            print("         ----> Te, Kmin, Kmax, Gamma, GE:", \
                  Te, Kmin, Kmax, Gamma, GE)

        sqrtK = ( mth.sqrt(Kmin) - mth.sqrt(Te/2.) * mth.log(1.-GE/Gamma))
        K     = sqrtK**2
        
        if self.getDebug():
            print("     <---- K:", K, " MeV")

        return K

    # Returns derived parameters for the calculation of the cumulative
    # probability.
    #   derivedParamters[0] - Gamma = Normalisation constant; 1/(integral
    #                                 of probability density function.
    #..  Used for Mode 0
    def getLaserCumProbParam(self):
        if self.getDebug():
            print(" Source(BeamLineElement).getLaserCumProbParam: start.")
            
        Te   = self.getParameters()[5]
        Kmin = self.getParameters()[6]
        Kmax = self.getParameters()[7]

        if self.getDebug():
            print("     ----> Te, Kmin, Kmax:", Te, Kmin, Kmax)
            
        intgrl = mth.sqrt(2./Te)*(mth.sqrt(Kmax) - mth.sqrt(Kmin))
        gamMAX = 1. - mth.exp(-intgrl)
        
        if self.getDebug():
            print("     ----> intgrl, gamMAX:", intgrl, gamMAX)
            
        Gamma  = 1./gamMAX

        if self.getDebug():
            print(" <---- Gamma:", Gamma)
            
            
        self.getderivedParameters().append(Gamma)

    # Returns kinetic PDF for kinetic energy distribution and kinetic
    # energy for 100 points between Kmin and Kmax.
    #..  Used for Mode 0
    def getLaserDrivenProtonEnergyProbDensity(self):
        if self.getDebug():
            print( \
        " Source(BeamLineElement).getLaserDrivenProtonEnergyProbDensity:")

        Te   = self.getParameters()[5]
        Kmin = self.getParameters()[6]
        Kmax = self.getParameters()[7]

        if self.getDebug():
            print("     ----> Te, Kmin, Kmmax:", Te, Kmin, Kmax)
        
        K   = np.linspace(Kmin, Kmax, 100)
        dK2 = (Kmax - Kmin) / 100. / 2.

        if self.getDebug():
            print("     ----> K, dK2:", K, dK2)
        
        # Approximate required distribution:
        eta = self.getderivedParameters()[0]
        g_K = []
        Ke  = []
        if self.getDebug():
            print("     ----> eta, g_K, Ke:", eta, g_K, Ke)

        for iK in range(len(K)):
            Ke.append( K[iK] + dK2 )
            
            g_K_val = (eta / mth.sqrt(Ke[-1])) \
                * mth.exp(-mth.sqrt(2. * Ke[-1] / Te))

            g_K.append( g_K_val )

        g_K /= np.sum(g_K)   # Normalize the probability distribution

        if self.getDebug():
            print("     ----> g_K, sum:", g_K, np.sum(g_K))

        return Ke, g_K
    
    # Returns cumulative PDF given kinetic energy (E).  Used to plot
    # cumulative probability for checks.
    #..  Used for Mode 0
    def getLaserCumProb(self, E):
        CumProb = 1.
        if E >= self.getParameters()[7]:
            return CumProb

        T_e      = self.getParameters()[5]
        E_min    = self.getParameters()[6]
        
        gam = mth.sqrt(2./T_e)*(mth.sqrt(E) - mth.sqrt(E_min))
        Gam = 1. - mth.exp(-gam)

        Gamma = self.getderivedParameters()[0]

        CumProb = Gamma * Gam

        return CumProb
    
    # Returns trace space given position (x,y), kinetic energy (K),
    # cosine of the polar angle (cTheta) and the azimuthal angle (Phi).
    # If xprime (xp) and yp (yprime) are given, these are used instead
    # of cTheta and Phi.
    #..  Used for Mode 0
    def getTraceSpace(self, x, y, K, cTheta, Phi, xp=None, yp=None):
        if self.getDebug():
            print(" Source(BeamLineElement).getTraceSpace: start.")
            print("     ----> x, y, K, cTheta, Phi:", \
                  x, y, K, cTheta, Phi)
            
        iRefPrtcl = BL.BeamLine.getcurrentReferenceParticle()
        particleMASS = PhysCnst.PhysicalConstants().getparticleMASS(   \
                                            iRefPrtcl.getSpecies() )
        p0        = iRefPrtcl.getMomentumIn(0)
        E0        = mth.sqrt( particleMASS**2 + p0**2)
        b0        = p0/E0
        if self.getDebug():
            print("     ----> p0, E0, b0, ( K0 ):", p0, E0, b0, \
                  "(", E0-particleMASS, ")")

        E = particleMASS+K
        p = mth.sqrt(E**2 - particleMASS**2)
        if self.getDebug():
            print("     ----> K, E, p:", K, E, p)

        if cTheta != None:
            sTheta = mth.sqrt(1.-cTheta**2)
            xPrime = sTheta * mth.cos(Phi) * p / p0
            yPrime = sTheta * mth.sin(Phi) * p / p0
        else:
            sTheta = None
            
        if xp != None:
            xPrime = xp * p / p0
        if yp != None:
            yPrime = yp * p / p0
            
        if self.getDebug():
            print("     ----> sTheta, xPrime, yPrime:", 
                  sTheta, xPrime, yPrime)

        delta     = (E - E0) / p0
        
        if self.getDebug():
            print("     ----> E, delta:", E, delta)

        TrcSpc = np.array([x, xPrime, y, yPrime, 0., delta])

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Trace space:", TrcSpc)

        return TrcSpc

#--------  Utilities:
    def tabulateParameters(self, filename="LaTeX.tex"):
        LTX.TableHeader(filename, '|l|c|l|', \
                        self.ModeText[self.getMode()])

        Line = "\\textbf{Parameter} & \\textbf{Value} & \\textbf{Unit}"
        LTX.TableLine(filename, Line)
        Line = "\\hline"
        LTX.TableLine(filename, Line)

        iPrm = 0
        for Prm in self.ParamLaTeX[self.getMode()]:
            Line = Prm + "&"                                     + \
                str(self.getParameters()[iPrm]) + "&" + \
                self.ParamUnit[self.getMode()][iPrm]
            iPrm += 1
            LTX.TableLine(filename, Line)

        LTX.TableTrailer(filename)

    @classmethod
    def CheckSourceParam(cls, _Mode, _Param):
        if cls.getDebug():
            print(' BeamLineElement(Source).CheckSourceParam:', \
                  ' mode, params:', \
                  _Mode, _Param)

        ValidMode  = False
        ValidParam = False
        
        ValidMode  = cls.CheckMode(_Mode)
        if ValidMode:
            ValidParam = cls.CheckParam(_Mode, _Param)

        ValidSourceParam = ValidMode and ValidParam
        if cls.getDebug():
            print("     <---- ValidMode, ValidParam, ValidSourcParam:", \
                  ValidMode, ValidParam, ValidSourceParam)

        if cls.getDebug():
            print(' <---- BeamLineElement(Source).CheckSourceParam,', \
                  ' done.  --------', \
                  '  --------  --------  --------  --------  --------')

        return ValidSourceParam

    @classmethod
    def CheckMode(cls, _Mode):
        if cls.getDebug():
            print(' BeamLineElement(Source).CheckMode: mode:', _Mode)
            
        #.. Do this with a loop, detault fail
        ValidMode = False
        for iMode in cls.ModeList:
            if cls.getDebug():
                print("     ----> iMode:", iMode, "cf _Mode", _Mode)
            
            if _Mode == iMode:
                ValidMode = True
                break
            
        if cls.getDebug():
            print("     <---- iMode, _Mode, ValidMode:", \
                  iMode, _Mode, ValidMode)
            
        if cls.getDebug():
            print(' <---- BeamLineElement(Source).CheckMode: done.', \
            '  --------  --------  --------  --------  --------  ')

        return ValidMode
    
    @classmethod
    def CheckParam(cls, _Mode, _Param):
        if cls.getDebug():
            print(' BeamLineElement(Source).CheckParam: mode, params:', \
                  _Mode, _Param)

        #.. Detault fail
        ValidParam = False
        if cls.getDebug():
            print("     ----> len(_Param):", len(_Param))
            print("         ----> cf     :", len(cls.ParamList[_Mode]))

        iMtch = 0
        if _Mode == 3 and len(_Param) == 0:
            if cls.getDebug():
                print("     ----> Mode, len(param):", _Mode, \
                      len(_Param))
        elif len(_Param) == len(cls.ParamList[_Mode]):
            for i in range(len(_Param)):
                if cls.getDebug():
                    print("         ----> i, ParamList, _Param:", \
                          i, cls.ParamList[_Mode][i], _Param[i])
                if isinstance(_Param[i], cls.ParamList[_Mode][i]) or \
                   _Param[i] == None:
                    iMtch += 1
                else:
                    if cls.getDebug():
                        print("             ----> No match for", \
                              "i, ParamList, _Param:", \
                              i, cls.ParamList[_Mode][i], _Param[i])
                    
            if cls.getDebug():
                print("     ----> iMtch:", iMtch)
        else:
            raise badParameters( \
                                 " BeamLineElement(Source).CheckParam:", \
                                 " bad source parameters. Exit")

        if iMtch == len(_Param):
            ValidParam = True
        if cls.getDebug():
            print("     <---- ValidParam:", ValidParam)
        if cls.getDebug():
            print(" <---- BeamLineElement(Source).CheckParam, done.", \
                  "  --------", \
                  '  --------  --------  --------  --------  --------')

        return ValidParam

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" Source(BeamLineElement).visualise: start")
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        cntr = self.getrStrt() + self.getStrt2End()/2.
        xlim = axs.get_xlim()
        ylim = axs.get_ylim()

        wdth = self.getLength()
        if wdth == 0.: wdth = (xlim[1] - xlim[0]) / 100.
        hght = (ylim[1] - ylim[0]) / 10.
        angl = 0.
        abt  = 'center'
        if CoordSys == "RPLC":
            if Proj == "xs":
                xy   = [cntr[2], cntr[0]]
            elif Proj == "ys":
                xy   = [cntr[2], cntr[1]]
        elif CoordSys == "Lab":
            if Proj == "xz":
                xy   = [cntr[2], cntr[0]]
            elif Proj == "yz":
                xy   = [cntr[2], cntr[1]]

        if self.getDebug():
            print("     ----> Centre:", cntr)
            print("     ---->     xy:", xy)
            print("     ---->   wdth:", wdth)
            print("     ---->   hght:", hght)
            print("     ---->   angl:", angl)
            print("     ---->    abt:", abt)
            
        Patch = patches.Ellipse(xy, wdth, hght, \
                                facecolor='gold', \
                                zorder=2)

        axs.add_patch(Patch)
                                   
    
#--------  I/o methods:
    def getLines(self):
        Lines = []

        Fields  = self.getName().split(":")
        Mode    = self.getMode()
        Text    = self.getModeText()
        Param   = self.getParameters()
        
        Stage   = Fields[1]
        Section = Fields[2]
        Element = Fields[2]

        Type    = Text
        Param   = "SourceMode"
        Value   = Mode
        Unit    = ""
        Comment = ""
        Lines.append([Stage, Section, Element, Type, \
                      Param, Value, Unit, Comment])

        for iPrm in range(len(self.getParameters())):
            Param = self.getParameterText()[iPrm]
            Value = self.getParameters()[iPrm]
            Unit  = self.getParameterUnit()[iPrm]
            Lines.append([Stage, Section, Element, Type, \
                          Param, Value, Unit, Comment])
            

        return Lines
    
    def writeElement(self, dataFILE):
        if self.getDebug():
            print(" Source(BeamLineElement).writeElement starts.")

        derivedCLASS = "Source"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Derived class:", bversion.decode('utf-8'))

        record = strct.pack(">i", self.getMode())
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Mode:", strct.unpack(">i", record))

        record = strct.pack(">i", len(self.getParameters()))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Number of parameters:", \
                  strct.unpack(">i", record))

        iPrm = 0
        if self.getDebug():
            print("     ----> Write parameters:")
        for typePrm in self.ParamList[self.getMode()]:
            if typePrm == float:
                record = strct.pack(">d", self.getParameters()[iPrm])
                dataFILE.write(record)
                if self.getDebug():
                    print("         ----> iPrm, value:", \
                          iPrm, strct.unpack(">d", record))
            else:
                record = strct.pack(">d", float(self.getParameters()[iPrm]))
                dataFILE.write(record)
                if self.getDebug():
                    print("         ----> iPrm, value:", \
                          iPrm, strct.unpack(">d", record))
            iPrm += 1
        if self.getDebug():
            print("     <---- Done.")
            
        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print(" <---- Source(BeamLineElement).writeElement done.")
                            
    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" Source(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        brecord = dataFILE.read(4)
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True
            
        record = strct.unpack(">i", brecord)
        Mode = record[0]
        if cls.getDebug():
            print("     ----> Mode:", Mode)

        brecord = dataFILE.read(4)
        if brecord == b'':
            if cls.getDebug():
                print(" <---- end of file, return.")
            return True
            
        record = strct.unpack(">i", brecord)
        nPrm = record[0]
        if cls.getDebug():
            print("     ----> Number of parameters:", nPrm)
            print("     ----> Data file version:", \
                  dataFILEinst.getdataFILEversion())

        """
        Data file version 6 or greater:
          Have laser-drivern source specified using the set of parameters
          defined in "revision 1" of the LhARA linear optics documentation.
          In this case the paramters are:
             [0] - Wavelength - microns
             [1] - Power - W
             [2] - Strehl ratio
             [3] - Focal spot radius - m
             [4] - Laser-pulse duration - s
             [5] - Hot electron temperature - MeV
             [6] - Minimum proton kinetic energy - MeV
             [7] - Maximum proton kinetic energy - MeV
             [8] - Intercept of sigma_{theta_S} at K=0 [degrees]
             [9] - Scaled slope of sigma_{theta_S}     [degrees]
            [10] - rp max

        Data file versions 5 and below:
          Laser-driven source specified by subset of paramters:
             [0] - Sigma of x gaussian - m
             [1] - Sigma of y gaussian - m
             [2] - Minimum cos theta to generate
             [3] - E_min: MeV min energy to generate
             [4] - E_max: MeV max energy to generate;
                     overwritten when calculated in
                     getLaserDrivenParticleEnergy
             [5] - nPnts: Number of points to sample for integration of PDF
             [6] - P_L: Laser power [W]
             [7] - E_L: Laser energy [J]
             [8] - lamda: Laser wavelength [um]
             [9] - t_laser: Laser pulse duration [s]
            [10] - d: Laser thickness [m]
            [11] - I: Laser intensity [W/cm2]
            [12] - theta_degrees: Electron divergence angle [degrees]
            [13] - Intercept of sigma_{theta_S} at K=0 [degrees]
            [14] - Scaled slope of sigma_{theta_S}     [degrees]
            [15] - rp max
          So, need to check against "ParamListOLD"

        """

        ParamListCHK = cls.ParamList
        Params       = []
        if dataFILEinst.getdataFILEversion() < 6:
            ParamListCHK = [ [float, float, float, float, float, int,    \
                               float, float, float, float, float, float, \
                               float, float, float, float],              \
                              [float, float, float, float, float],       \
                              [float, float, float, float, float],       \
                              [],                                        \
                              [float, float, float] ]
            
        for iPrm in range(nPrm):
            brecord = dataFILE.read((1*8))
            if brecord == b'':
                return True
        
            record  = strct.unpack(">d", brecord)
            var     = float(record[0])
            if cls.getDebug():
                print(" iPrm, var:", iPrm, var)
            if ParamListCHK[Mode][iPrm] == int:
                var = int(var)
            Params.append(var)
            
        if cls.getDebug():
            print("     ----> Parameters:", Params)

        return EoF, Mode, Params


"""
To do:
------
 - put z/s position of components of doublet
 - put length of the doublet

Derived class QuadDoublet:
==========================

  QuadDoublet class derived from BeamLineElement to contain parameters
  for an quad doublet.


  Class attributes:
  -----------------
    instances : List of instances of QuadDoublet class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
  Calling arguments:
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.

  _TrnsMtrx : Transfer matrix:


  Instance attributes to define quadrupole:
  -----------------------------------------
  _xxx  : 
  _yyy  : 

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising quad
                parameterrs.

  Set methods:
     setDebug: Set debug flag
          Input: Debug (bool) 

  Get methods:
      getDebug

  Utilities:

"""
class QuadDoublet(BeamLineElement):
    instances = []
    __Debug   = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None,
                 _FDorDF=None, _Q1par=None, _d=None, _Q2par=None):
        
        if self.getDebug():
            print(' QuadDoublet.__init__: ', \
                  'creating the QuadDoublet object: ')

        QuadDoublet.instances.append(self)

        self.setAll2None()

        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, _Name, _rStrt, _vStrt, _drStrt, \
                                 _dvStrt)

        if _FDorDF != "FD" and _FDorDF != "DF":
            raise badBeamLineElement("QuadDoublet: bad specification", \
                                     " for FDorDF")

        if isinstance(_Q1par,list):
            if len(_Q1par) != 3:
                raise badBeamLineElement("QuadDoublet: bad specification", \
                                         " for Q1par")
        else:
            raise badBeamLineElement("QuadDoublet: bad specification", \
                                     " for Q1par")

        if not(isinstance(_d, float)):
            raise badBeamLineElement("QuadDoublet: bad specification", \
                                     " for Q1par")

        if isinstance(_Q2par,list):
            if len(_Q2par) != 3:
                raise badBeamLineElement("QuadDoublet: bad specification", \
                                         " for Q2par")
        else:
            raise badBeamLineElement("QuadDoublet: bad specification", \
                                     " for Q2par")
        
        self.setFDorDF(_FDorDF)
        self.setQ1par(_Q1par)
        self.setSeparation(_d)
        self.setQ2par(_Q2par)

        if self.getFDorDF() == "FD":
            iQ1 = FocusQuadrupole(_Name+":FQ1", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ1par()[0], self.getQ1par()[1], self.getQ1par()[2])
        else:
            iQ1 = DefocusQuadrupole(_Name+":DQ1", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ1par()[0], self.getQ1par()[1], self.getQ1par()[2])
        iD = Drift(_Name+":sep", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                            self.getSeparation())
        if self.getFDorDF() == "FD":
            iQ2 = DefocusQuadrupole(_Name+":DQ2", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ2par()[0], self.getQ2par()[1], self.getQ2par()[2])
        else:
            iQ2 = FocusQuadrupole(_Name+":FQ2", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ2par()[0], self.getQ2par()[1], self.getQ2par()[2])

        if self.getDebug():
            print("     ----> Dump componenets:")
            print(iQ1)
            print(iD)
            print(iQ2)

        BeamLineElement.removeInstance(iQ1)
        BeamLineElement.removeInstance(iD)
        BeamLineElement.removeInstance(iQ2)
        self.setQ1(iQ1)
        self.setD(iD)
        self.setQ2(iQ2)
                   
        if self.getDebug():
            print("     ----> New QuadDoublet instance: \n", self)
            print(" <---- Done.")

    def __repr__(self):
        return "QuadDoublet()"

    def __str__(self):
        print(" QuadDoublet:")
        print(" ------------")
        print("     ---->     Debug flag:", QuadDoublet.getDebug())
        print("     ---->         FDorDF:", self.getFDorDF())
        print("     ---->          Q1par:", self.getQ1par())
        print("     ---->     Separation:", self.getSeparation())
        print("     ---->          Q2par:", self.getQ2par())
        print("     ---->     Components:", self.getQ1().getName(), \
              self.getD().getName(), self.getQ2().getName())
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- QuadDoublet parameter dump complete."

    def SummaryStr(self):
        Str  = "QuadDoublet  : " + BeamLineElement.SummaryStr(self) + \
            "; FDorDF = ", self.getFDorDF() + \
            "; Q1par = ", str(self.getQ1par()) + \
            "; separation = ", str(self.getSeparation()) + \
            "; Q2par = ", str(self.getQ2par())
        return Str

    
# -------- "Set methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug
        
    def setAll2None(self):
        self._FDorDF     = None
        self._Q1par      = None
        self._Separation = None
        self._Q2par      = None
        self._TrnsMtrx   = None
        

    def setFDorDF(self, _FDorDF):
        if _FDorDF != "FD" and _FDorDF != "DF":
            raise badParameter( \
                    "BeamLineElement.QuadDoublet.setFDorDF:", \
                                " bad FDorDF:", _FDorDF)
        self._FDorDF = _FDorDF
        
    def setSeparation(self, _d):
        if not(isinstance(_d, float)):
            raise badParameter( \
                    "BeamLineElement.QuadDoublet.setSeparation:", \
                    " bad separation:", _d)
               
        self._Separation = _d
        
    def setQ1par(self, _Q1par):
        if isinstance(_Q1par,list):
            if len(_Q1par) != 3:
                raise badBeamLineElement(
                    "BeamLineElement.QuadDoublet.setQ1par:", \
                    " for Q1par")
        else:
            raise badBeamLineElement(
                "BeamLineElement.QuadDoublet.setQ1par:", \
                " for Q1par")
        
        self._Q1par = _Q1par
        
    def setQ2par(self, _Q2par):
        if isinstance(_Q2par,list):
            if len(_Q2par) != 3:
                raise badBeamLineElement(
                    "BeamLineElement.QuadDoublet.setQ2par:", \
                    " for Q2par")
        else:
            raise badBeamLineElement(
                "BeamLineElement.QuadDoublet.setQ2par:", \
                " for Q2par")
        
        self._Q2par = _Q2par

    def setQ1(self, iQ1):
        if not isinstance(iQ1, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadDoublet.setQ1:", \
                " not a beamline element")
        self._iQ1 = iQ1
            
    def setD(self, iD):
        if not isinstance(iD, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadDoublet.setD:", \
                " not a beamline element")
        self._iD = iD
            
    def setQ2(self, iQ2):
        if not isinstance(iQ2, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadDoublet.setQ2:", \
                " not a beamline element")
        self._iQ2 = iQ2
            
    def setTransferMatrix(self, _R):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        
        if self.getDebug():
            print(" QuadDoublet(BeamLineElement).setTransferMatrix:")

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Trace space:", _R)

        self.getQ1().setTransferMatrix(_R)
        TrnsfrQ1 = self.getQ1().getTransferMatrix()

        self.getD().setTransferMatrix()
        TrnsfrD  = self.getD().getTransferMatrix()

        self.getQ2().setTransferMatrix(_R) 
        TrnsfrQ2 = self.getQ2().getTransferMatrix()

        if self.getDebug():
            print("     ----> Transfer matrix for Q1:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrQ1)
            print("     ----> Transfer matrix for D:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrD)
            print("     ----> Transfer matrix for Q2:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrQ2)

        
        TrnsMtrx = TrnsfrD.dot(TrnsfrQ1)
        TrnsMtrx = TrnsfrQ2.dot(TrnsMtrx)

        if self.getDebug():
            print("     ----> Transfer matrix for QuadDoublet:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

        
# -------- "Get methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug

    def getFDorDF(self):
        return self._FDorDF

    def getQ1par(self):
        return self._Q1par

    def getSeparation(self):
        return self._Separation

    def getQ2par(self):
        return self._Q2par

    def getQ1(self):
        return self._iQ1

    def getD(self):
        return self._iD

    def getQ2(self):
        return self._iQ2

    
# -------- Utilities:
    
    
"""
Derived class QuadTriplet:
==========================

  QuadTriplet class derived from BeamLineElement to contain parameters
  for an quad doublet.


  Class attributes:
  -----------------
    instances : List of instances of QuadTriplet class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
  Calling arguments:
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of principal axis of element.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.

  _TrnsMtrx : Transfer matrix:


  Instance attributes to define quadrupole:
  -----------------------------------------
  _xxx  : 
  _yyy  : 

    
  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising quad
                parameterrs.

  Set methods:
     setDebug: Set debug flag
          Input: Debug (bool) 

  Get methods:
      getDebug

  Utilities:

"""
class QuadTriplet(BeamLineElement):
    instances = []
    __Debug   = False

    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _FDForDFD=None, _Q1par=None, _d1=None, _Q2par=None,     \
                                              _d2=None, _Q3par=None     ):
        if self.getDebug():
            print(' QuadTriplet.__init__: ', \
                  'creating the QuadTriplet object: ')

        QuadTriplet.instances.append(self)

        self.setAll2None()

        # BeamLineElement class initialization:
        BeamLineElement.__init__(self, _Name, _rStrt, _vStrt, _drStrt, \
                                 _dvStrt)

        if _FDForDFD != "FDF" and _FDForDFD != "DFD":
            raise badBeamLineElement("QuadTriplet: bad specification", \
                                     " for FDForDFD")

        if isinstance(_Q1par,list):
            if len(_Q1par) != 3:
                raise badBeamLineElement("QuadTriplet: bad specification", \
                                         " for Q1par")
        else:
            raise badBeamLineElement("QuadTriplet: bad specification", \
                                     " for Q1par")

        if not(isinstance(_d1, float)):
            raise badBeamLineElement("QuadTriplet: bad specification", \
                                     " for d1")

        if isinstance(_Q2par,list):
            if len(_Q2par) != 3:
                raise badBeamLineElement("QuadTriplet: bad specification", \
                                         " for Q2par")
        else:
            raise badBeamLineElement("QuadTriplet: bad specification", \
                                     " for Q2par")
        
        if not(isinstance(_d2, float)):
            raise badBeamLineElement("QuadTriplet: bad specification", \
                                     " for d2")

        if isinstance(_Q2par,list):
            if len(_Q3par) != 3:
                raise badBeamLineElement("QuadTriplet: bad specification", \
                                         " for Q3par")
        else:
            raise badBeamLineElement("QuadTriplet: bad specification", \
                                     " for Q3par")
        
        self.setFDForDFD(_FDForDFD)
        self.setQ1par(_Q1par)
        self.setSeparation1(_d1)
        self.setQ2par(_Q2par)
        self.setSeparation2(_d2)
        self.setQ3par(_Q3par)

        lenTRPLT = 0.
        if self.getFDForDFD() == "FDF":
            iQ1 = FocusQuadrupole(_Name+":FQ1", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ1par()[0], self.getQ1par()[1], self.getQ1par()[2])
        else:
            iQ1 = DefocusQuadrupole(_Name+":DQ1", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ1par()[0], self.getQ1par()[1], self.getQ1par()[2])
        iD1 = Drift(_Name+":sep1", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                            self.getSeparation1())
        lenTRPLT += iQ1.getLength() + iD1.getLength()
        if self.getFDForDFD() == "FDF":
            iQ2 = DefocusQuadrupole(_Name+":DQ2", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ2par()[0], self.getQ2par()[1], self.getQ2par()[2])
        else:
            iQ2 = FocusQuadrupole(_Name+":FQ2", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ2par()[0], self.getQ2par()[1], self.getQ2par()[2])
        iD2 = Drift(_Name+":sep2", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                            self.getSeparation2())
        lenTRPLT += iQ2.getLength() + iD2.getLength()
        if self.getFDForDFD() == "FDF":
            iQ3 = FocusQuadrupole(_Name+":FQ3", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ3par()[0], self.getQ3par()[1], self.getQ3par()[2])
        else:
            iQ3 = DefocusQuadrupole(_Name+":DQ3", \
                            _rStrt, _vStrt, _drStrt, _dvStrt, \
                self.getQ3par()[0], self.getQ3par()[1], self.getQ3par()[2])
        lenTRPLT += iQ3.getLength()
        
        if self.getDebug():
            print("     ----> Dump componenets:")
            print(iQ1)
            print(iD1)
            print(iQ2)
            print(iD2)
            print(iQ3)
            print("     <---- Total length:", lenTRPLT)

        self.setLength(lenTRPLT)
        self.setStrt2End(np.array([0., 0., self.getLength()]))
        self.setRot2LbEnd(self.getRot2LbStrt())

        BeamLineElement.removeInstance(iQ1)
        BeamLineElement.removeInstance(iD1)
        BeamLineElement.removeInstance(iQ2)
        BeamLineElement.removeInstance(iD2)
        BeamLineElement.removeInstance(iQ3)
        self.setQ1(iQ1)
        self.setD1(iD1)
        self.setQ2(iQ2)
        self.setD2(iD2)
        self.setQ3(iQ3)
                   
        if self.getDebug():
            print("     ----> New QuadTriplet instance: \n", self)
            print(" <---- Done.")

    def __repr__(self):
        return "QuadTriplet()"

    def __str__(self):
        print(" QuadTriplet:")
        print(" ----------------")
        print("     ---->     Debug flag:", QuadTriplet.getDebug())
        print("     ---->         FDorDF:", self.getFDForDFD())
        print("     ---->          Q1par:", self.getQ1par())
        print("     ---->    Separation1:", self.getSeparation1())
        print("     ---->          Q2par:", self.getQ2par())
        print("     ---->    Separation2:", self.getSeparation2())
        print("     ---->          Q3par:", self.getQ3par())
        print("     ---->     Components:", self.getQ1().getName(), \
              self.getD1().getName(), self.getQ2().getName(), \
              self.getD2().getName(), self.getQ3().getName(), \
              )
        with np.printoptions(linewidth=500,precision=7,suppress=True):
            print("     ----> Transfer matrix: \n", self.getTransferMatrix())
        BeamLineElement.__str__(self)
        return " <---- QuadTriplet parameter dump complete."

    def SummaryStr(self):
        Str  = "QuadTriplet  : " + BeamLineElement.SummaryStr(self) + \
            "; FDForDFD = ", self.getFDForDFD() + \
            "; Q1par = ", str(self.getQ1par()) + \
            "; separation1 = ", str(self.getSeparation1()) + \
            "; Q2par = ", str(self.getQ2par()) + \
            "; separation2 = ", str(self.getSeparation2()) + \
            "; Q3par = ", str(self.getQ3par())
        return Str

    
# -------- "Set methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug
        
    def setAll2None(self):
        self._FDForDFF    = None
        self._Q1par       = None
        self._Separation1 = None
        self._Q2par       = None
        self._Separation2 = None
        self._Q3par       = None
        self._TrnsMtrx    = None

    def setFDForDFD(self, _FDForDFD):
        if _FDForDFD != "FDF" and _FDForDFD != "DFD":
            raise badParameter( \
                    "BeamLineElement.QuadDoublet.setFDForDFD:", \
                                " bad FDForDFD:", _FDForDFD)
        self._FDForDFD = _FDForDFD
        
    def setQ1par(self, _Q1par):
        if isinstance(_Q1par,list):
            if len(_Q1par) != 3:
                raise badBeamLineElement(
                    "BeamLineElement.QuadTriplet.setQ1par:", \
                    " for Q1par")
        else:
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setQ1par:", \
                " for Q1par")
        
        self._Q1par = _Q1par
        
    def setSeparation1(self, _d1):
        if not(isinstance(_d1, float)):
            raise badParameter( \
                    "BeamLineElement.QuadDoublet.setSeparation:", \
                    " bad separation 1:", _d1)
               
        self._Separation1 = _d1
        
    def setQ2par(self, _Q2par):
        if isinstance(_Q2par,list):
            if len(_Q2par) != 3:
                raise badBeamLineElement(
                    "BeamLineElement.QuadTriplet.setQ2par:", \
                    " for Q2par")
        else:
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setQ2par:", \
                " for Q2par")
        
        self._Q2par = _Q2par

    def setSeparation2(self, _d2):
        if not(isinstance(_d2, float)):
            raise badParameter( \
                    "BeamLineElement.QuadDoublet.setSeparation:", \
                    " bad separation 2:", _d2)
               
        self._Separation2 = _d2
        
    def setQ3par(self, _Q3par):
        if isinstance(_Q3par,list):
            if len(_Q3par) != 3:
                raise badBeamLineElement(
                    "BeamLineElement.QuadTriplet.setQ3par:", \
                    " for Q3par")
        else:
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setQ3par:", \
                " for Q3par")
        
        self._Q3par = _Q3par

    def setQ1(self, iQ1):
        if not isinstance(iQ1, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setQ1:", \
                " not a beamline element")
        self._iQ1 = iQ1
            
    def setD1(self, iD1):
        if not isinstance(iD1, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setD1:", \
                " not a beamline element")
        self._iD1 = iD1
            
    def setQ2(self, iQ2):
        if not isinstance(iQ2, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setQ2:", \
                " not a beamline element")
        self._iQ2 = iQ2
            
    def setD2(self, iD2):
        if not isinstance(iD2, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setD2:", \
                " not a beamline element")
        self._iD2 = iD2
            
    def setQ3(self, iQ3):
        if not isinstance(iQ3, BeamLineElement):
            raise badBeamLineElement(
                "BeamLineElement.QuadTriplet.setQ3:", \
                " not a beamline element")
        self._iQ3 = iQ3
            
    def setTransferMatrix(self, _R):
        if BeamLineElement.getNOdispersion():
            _R = np.array([0., 0., 0., 0., 0., 0.])
        
        if self.getDebug():
            print(" QuadTriplet(BeamLineElement).setTransferMatrix:")

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Trace space:", _R)

        self.getQ1().setTransferMatrix(_R)
        TrnsfrQ1 = self.getQ1().getTransferMatrix()

        self.getD1().setTransferMatrix()
        TrnsfrD1 = self.getD1().getTransferMatrix()

        self.getQ2().setTransferMatrix(_R) 
        TrnsfrQ2 = self.getQ2().getTransferMatrix()

        self.getD2().setTransferMatrix()
        TrnsfrD2 = self.getD2().getTransferMatrix()

        self.getQ3().setTransferMatrix(_R) 
        TrnsfrQ3 = self.getQ3().getTransferMatrix()

        if self.getDebug():
            print("     ----> Transfer matrix for Q1:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrQ1)
            print("     ----> Transfer matrix for D1:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrD1)
            print("     ----> Transfer matrix for Q2:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrQ2)
            print("     ----> Transfer matrix for D2:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrD2)
            print("     ----> Transfer matrix for Q3:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsfrQ3)

        TrnsMtrx = np.matmul(TrnsfrD2, TrnsfrQ3)
        if self.getDebug():
            print("     ----> Transfer matrix for D2*Q3:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)
        TrnsMtrx = np.matmul(TrnsfrQ2, TrnsMtrx)
        if self.getDebug():
            print("     ----> Transfer matrix for Q2*D2*Q3:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)
        TrnsMtrx = np.matmul(TrnsfrD1, TrnsMtrx)
        if self.getDebug():
            print("     ----> Transfer matrix for D1*Q2*D2*Q3:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)
        TrnsMtrx = np.matmul(TrnsfrQ1, TrnsMtrx)
        if self.getDebug():
            print("     ----> Transfer matrix for Q1*D1*Q2*D2*Q3:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)
                
        if self.getDebug():
            print("     ----> Transfer matrix for QuadTriplet:")
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(TrnsMtrx)

        self._TrnsMtrx = TrnsMtrx

        
# -------- "Get methods"
# Methods believed to be self-documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug

    def getFDForDFD(self):
        return self._FDForDFD
        
    def getQ1par(self):
        return self._Q1par
        
    def getSeparation1(self):
        return self._Separation1
        
    def getQ2par(self):
        return self._Q2par

    def getSeparation2(self):
        return self._Separation2
        
    def getQ3par(self):
        return self._Q3par

    def getQ1(self):
        return self._iQ1
            
    def getD1(self):
        return self._iD1
            
    def getQ2(self):
        return self._iQ2
            
    def getD2(self):
        return self._iD2
            
    def getQ3(self):
        return self._iQ3

    
# -------- Utilities:
#
#..  Thin lens approximation functions:
    @classmethod
    def getTLAr(cls, _d=None, _f=None):
        if cls.getDebug():
            print(" BeamLineElement.QuadTriplet.TLA_r:", \
                  " inter-quad separation, focal length:", _d, _f)

        r = None

        if _d == None or _f == None:
            raise badParameter( \
                        "BeamLineElement.QuadTriplet.TLA_r:, d, f :" + \
                        str(_d) + ", " + str(_f))
        s1 = _f
        s2 = _d

        r1 = s2**2 * (2.*s1**2 + 3.*s1*s2 + 2.*s2**2)**2
        r2 = 8.*s1**2*s2**3*(s1+s2)

        if r1 < r2:
            r = None
        else:
            r = mth.sqrt(r1 - r2)

        return r
    
    @classmethod
    def getTLAr(cls, _d=None, _f=None):
        if cls.getDebug():
            print(" BeamLineElement.QuadTriplet.TLA_r:", \
                  " inter-quad separation, focal length:", _d, _f)

        r = None

        if _d == None or _f == None:
            raise badParameter( \
                        "BeamLineElement.QuadTriplet.TLA_r:, d, f :" + \
                        str(_d) + ", " + str(_f))
        s1 = _f
        s2 = _d

        r1 = s2**2 * (2.*s1**2 + 3.*s1*s2 + 2.*s2**2)**2
        r2 = 8.*s1**2*s2**3*(s1+s2)

        if r1 < r2:
            r = None
        else:
            r = mth.sqrt(r1 - r2)

        return r

    @classmethod
    def getTLAf1(cls, _d=None, _f=None):

        f1 = None

        if _d == None or _f == None:
            raise badParameter( \
                        "BeamLineElement.QuadTriplet.TLA_r:, d, f :" + \
                        str(_d) + ", " + str(_f))
        s1 = _f
        s2 = _d

        r = cls.getTLAr(s2, s1)

        t = -r + s2*(2.*s1**2 + 3.*s1*s2 +2.*s2**2)
        b = 2.*(s1+s2)

        f1 = mth.sqrt(t/b)

        return f1

    @classmethod
    def getTLAf1nf2(cls, _d=None, _f=None, FDF=True):

        f1 = None

        if _d == None or _f == None:
            raise badParameter( \
                        "BeamLineElement.QuadTriplet.TLA_r:, d, f :" + \
                        str(_d) + ", " + str(_f))
        s1 = _f
        s2 = _d

        r  = cls.getTLAr(s2, s1)
        f1 = cls.getTLAf1(s2, s1)

        sgn = -1.
        if not FDF: sgn = 1

        t = r - s2*(2.*s1**2 + 5.*s1*s2 +2.*s2**2)
        b = 4.*(s1+s2) * sgn*f1

        f2 = t/b

        return f1, f2

    
"""
Derived class RPLCswitch:
=========================

  RPLCswitch class derived from BeamLineElement to manage a change in the
  relation of RPLC coordinates to laboratiry coordinates.  __init__
  sets parameters of the RPLCswitch.  The transfer matrix is the
  effective transformation from the previous RPLC->laboratory
  relationship to the new one.


  Class attributes:
  -----------------
    instances : List of instances of BeamLineElement class
  __Debug     : Debug flag


  Parent class instance attributes:
  ---------------------------------
  Calling arguments:
   _Name : Name
   _rStrt : numpy array; x, y, z position (in m) of start of element.
   _vStrt : numpy array; theta, phi of y and z axes of new RPLC
            coordinate system.
  _drStrt : "error", displacement of start from nominal position.
  _dvStrt : "error", deviation in theta and phy from nominal axis.

  _Length   : 0
  _TrnsMtrx : Transfer matrix.


  Methods:
  --------
  Built-in methods __init__, __repr__ and __str__.
      __init__ : Creates instance of beam-line element class.
      __repr__: One liner with call.
      __str__ : Dump of constants

    SummaryStr: No arguments, returns one-line string summarising RPLCswitch
                parameterrs.

  Set methods:
 setTransferMatrix : "Calculate" amd set transfer matrix.

  Get methods:
     getDebug: get debug flag, bool
    getLength: Returns length of RPLCswitch (presently 0)

  Utilities:


"""
class RPLCswitch(BeamLineElement):
    instances  = []
    __Debug    = False

#--------  "Built-in methods":    
    def __init__(self, _Name=None, \
                 _rStrt=None, _vStrt=None, _drStrt=None, _dvStrt=None, \
                 _3Drotation=False):
        if self.getDebug():
            print(' RPLCswitch.__init__: ', \
                  'creating the RPLCswitch object')
            print("     ----> vStrt:", _vStrt)

        RPLCswitch.instances.append(self)
        
        #.. BeamLineElement class initialisation:
        BeamLineElement.__init__(self, _Name, _rStrt, _vStrt, _drStrt, _dvStrt)

        self.setStrt2End(np.array([0., 0., 0.]))
        self.setRot2LbEnd(self.getRot2LbStrt())

        self.set3Drotation(_3Drotation)
        
        self.setTransferMatrix()
                
        if self.getDebug():
            print("     ----> New RPLCswitch instance: \n", \
                  self)
            
    def __repr__(self):
        return "RPLCswitch()"
    
    def __str__(self):
        print(" RPLCswitch:")
        print(" -----------")
        print("     ----> Debug flag:", RPLCswitch.getDebug())
        print("     ----> 3Drotation:", self.get3Drotation())
        BeamLineElement.__str__(self)
        return " <---- RPLCswitch parameter dump complete."

    def SummaryStr(self):
        Str  = "RPLCswitch       : " + BeamLineElement.SummaryStr(self)
        return Str

    
#--------  "Set methods".
#.. Methods believed to be self documenting(!)
    @classmethod
    def setDebug(cls, Debug=False):
        cls.__Debug = Debug

    def set3Drotation(self, _3Drotation):
        if not isinstance(_3Drotation, bool):
            raise badParameter()

        self._3Drotation = _3Drotation
        
    def setTransferMatrix(self):
        if self.getDebug():
            print(" RPLC(BeamLineElement).setTransferMatrix; start:")

        iLst  = BeamLineElement.getinstances()[ \
                                len(BeamLineElement.getinstances())-2 \
                                               ]
        if self.get3Drotation():
            invRE = np.linalg.inv(iLst.getRot2LbEnd())

            effctvRot = np.matmul(invRE, self.getRot2LbStrt())

            offDiag   = np.zeros((3,3))

            TrnsMtrx = np.block([ \
                                  [effctvRot, offDiag],  \
                                  [  offDiag, effctvRot] \
                                 ])
        else:
            if self.getDebug():
                with np.printoptions(linewidth=500,precision=7,suppress=True):
                    print("     ----> vStrt(iLst), vStrt:", \
                          iLst.getvStrt(), self.getvStrt())
                    
            PhiLst   = iLst.getvStrt()[0][1]
            Phi      = self.getvStrt()[0][1]
            if self.getDebug():
                print("     ----> PhiLst:", PhiLst, \
                      "\n              Phi:", Phi)
                
            dPhi     = Phi - PhiLst
            cdPhi    = mth.cos(dPhi)
            sdPhi    = mth.sin(dPhi)

            TrnsMtrx = np.array([ \
                                  [cdPhi,   0., -sdPhi,     0., 0., 0.], \
                                  [  0., cdPhi,     0., -sdPhi, 0., 0.], \
                                  [sdPhi,   0.,  cdPhi,     0., 0., 0.], \
                                  [  0., sdPhi,   0.,    cdPhi, 0., 0.], \
                                  [  0.,    0.,   0.,       0., 1., 0.], \
                                  [  0.,    0.,   0.,       0., 0., 1.]  \
                                 ])
                             
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" <---- RPLCswitch(BeamLineElement): ", \
                      "returns Transfer matrix: \n", \
                      TrnsMtrx)
                
        self._TrnsMtrx = TrnsMtrx

    def calcRot2LbEnd(self):
        if self.getDebug():
            print(" RPLC(BeamLineElement).calcRot2LbEnd; start:")

        iLst  = BeamLineElement.getinstances()[ \
                                len(BeamLineElement.getinstances())-2 \
                                               ]
            
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> vStrt(iLst), vStrt:", \
                          iLst.getvStrt(), self.getvStrt())
                
        PhiLst   = iLst.getvStrt()[0][1]
        Phi      = self.getvStrt()[0][1]
        if self.getDebug():
            print("     ----> PhiLst:", PhiLst, \
                  "              Phi:", Phi)
                
        dPhi     = Phi - PhiLst
        cdPhi    = mth.cos(dPhi)
        sdPhi    = mth.sin(dPhi)

        RotMtrx = np.array([ \
                             [ cdPhi,   sdPhi, 0.], \
                             [-sdPhi,   cdPhi, 0.], \
                             [  0.,        0., 1.] \
                            ])
        
        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("     ----> Rotation matrix: \n", RotMtrx)
                
        return RotMtrx

#--------  "get methods"
#.. Methods believed to be self documenting(!)
    @classmethod
    def getDebug(cls):
        return cls.__Debug
    
    def getLength(self):
        return 0.

    def get3Drotation(self):
        return self._3Drotation

    def visualise(self, axs, CoordSys, Proj):
        if self.getDebug():
            print(" RPLCswitch(BeamLineElement).visualise: start")
            print("     ----> self.getrStrt():", self.getrStrt())
            print("     ----> self.getStrt2End():", self.getStrt2End())

        if CoordSys == "Lab":
            return


        cntr = self.getrStrt()
        ylim = axs.get_ylim()

        if self.getDebug():
            print("     ----> ylim:", ylim)

        sz = cntr[2]
        xyup  = ylim[1]
        xydn  = ylim[0]
        

        if self.getDebug():
            print("     ----> Centre:", cntr)
            print("     ---->     sz:", sz)
            print("     ---->   xyup:", xyup)
            print("     ---->   xydn:", xydn)

        axs.axvline(sz, xyup, 1., \
                    color="seagreen", \
                    linewidth=0.2, \
                    zorder=2)
        axs.axvline(sz, xydn, 0., \
                    color="seagreen", \
                    linewidth=0.2, \
                    zorder=2)
                                   
        if self.getDebug():
            print(" RPLCswitch(BeamLineElement).visualise: start")

    
#--------  I/o methods:               <--------  Here
    def writeElement(self, dataFILE):
        if self.getDebug():
            print(" RPLCswitch(BeamLineElement).writeElement starts.")

        derivedCLASS = "RPLCswitch"
        bversion = bytes(derivedCLASS, 'utf-8')
        record   = strct.pack(">i", len(derivedCLASS))
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Length of derived class record:", \
                  strct.unpack(">i", record))
        record   = bversion
        dataFILE.write(record)
        if self.getDebug():
            print("     ----> Derived class:", bversion.decode('utf-8'))

        if self.getDebug():
            print("     <---- Done.")
            
        BeamLineElement.writeElement(self, dataFILE)
        
        if self.getDebug():
            print(" <---- RPLCswitch(BeamLineElement).writeElement done.")
        
    @classmethod
    def readElement(cls, dataFILEinst):
        if cls.getDebug():
            print(" RPLCswitch(BeamLineElement).readElement starts.")

        dataFILE = dataFILEinst.getdataFILE()

        EoF = False

        if cls.getDebug():
            print(" <---- end of file, return.")
            
        return EoF


#--------  Utilities:
    def Transport(self, _R=None):
        if not isinstance(_R, np.ndarray) or np.size(_R) != 6:
            raise badParameter( \
            " RPLCswitch(BeamLineElement).Transport:" + \
                                "bad input vector:", \
                                _R)

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" RPLCswitch(BeamLineElement).Transport: \n", \
                      "     ----> input trace space:", _R)
            print("          ----> Outside:", self.OutsideBeamPipe(_R))
            print("          ----> 3D-rotatation:", self.get3Drotation())
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print("         ----> Transfer matrix: \n", \
                      self.getTransferMatrix())

            
        if self.OutsideBeamPipe(_R) or \
           self.ExpansionParameterFail(_R) or \
           abs(_R[4]) > 2.5:
            if abs(_R[4]) < 2.5:
                print(" RPLCswitch.Transport:       outside beam pipe:", \
                                         self.OutsideBeamPipe(_R), "\n", \
                      "                       expansion paramter fail:", \
                                  self.ExpansionParameterFail(_R), "\n", \
                      "                                 TraceSpace[4]:", \
                                                                _R[4]    )
            _Rprime = None
        else:
            if self.get3Drotation():
                phsSpc      = \
                    Prtcl.Particle.RPLCTraceSpace2PhaseSpace(_R).reshape(6)
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7, \
                                         suppress=True):
                        print("     ----> PhaseSpace      :", phsSpc) 
                phsSpcprime = self.getTransferMatrix().dot(phsSpc)
                if self.getDebug():
                    with np.printoptions(linewidth=500,precision=7, \
                                         suppress=True):
                        print("     ----> PhaseSpace prime:", phsSpcprime)
                _Rprime     = \
                    Prtcl.Particle.RPLCPhaseSpace2TraceSpace(phsSpcprime)
            else:
                detTrnsfrMtrx = np.linalg.det(self.getTransferMatrix())
                error         = abs(1. - abs(detTrnsfrMtrx))
                if error > 1.E-6:
                    print(" RPLCswitch(BeamLineElement).Transport:(4):" \
                          " detTrnsfrMtrx:", detTrnsfrMtrx)
                _Rprime = self.getTransferMatrix().dot(_R)

        if self.getDebug():
            with np.printoptions(linewidth=500,precision=7,suppress=True):
                print(" <---- Rprime:", _Rprime)

        if not isinstance(_Rprime, np.ndarray):
            pass

        return _Rprime

    
#--------  Exceptions:
class badBeamLineElement(Exception):
    pass

class badParameter(Exception):
    pass

class badSourceSpecification(Exception):
    pass

class secondFacility(Exception):
    pass

class badParameters(Exception):
    pass

class ReferenceParticleNotSpecified(Exception):
    pass

class FailToCreateTraceSpaceAtSource(Exception):
    pass

class KillInfiniteLoop(Exception):
    pass

class cantFINDfile(Exception):
    pass

class need2convertdvStrt(Exception):
    pass

class badSOURCEmode(Exception):
    pass

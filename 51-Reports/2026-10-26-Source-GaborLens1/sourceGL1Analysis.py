#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sourceGL1Analysis.py:
=====================

  Simulation of the LhARA Stage 1 beam line (laser-driven source and Gabor
  lenses, 11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv) and analysis
  of the laser-driven proton source and of the first Gabor lens.

  Usage, from the top directory of the repository after "source startup.bash":
    python <this script> simulate [-n nEvts] [-s seed] [-d datafile]
    python <this script> analyse  [-d datafile] [-o outputdir]

  simulate : track nEvts events through the whole beam line with a fixed
             random seed, so that the run can be reproduced.  The BeamIO
             data file is written to datafile.
  analyse  : read datafile; write tables to <outputdir>/results and figures
             to <outputdir>/figures (default outputdir: this directory).

  Inside drifts and Gabor lenses the beam is propagated with the transfer
  matrices of the beam-line elements themselves (Gabor lens matrices are
  evaluated by the code for each particle energy), so that the envelope can
  be plotted inside the elements; the result at the exit of each element is
  checked against the trace space recorded by the simulation.

Created on Mon 26Oct26: Version history:
----------------------------------------
 1.0: 26Oct26: First implementation

"""

import os
import sys
import csv
import json
import math
import random
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

import BeamIO            as bmIO
import BeamLine          as BL
import BeamLineElement   as BLE
import Particle          as Prtcl
import PhysicalConstants as PhysCnst

BEAMLINE = "11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv"
DATAFILE = "99-Scratch/LhARA-LsrDrvn-Gabor.dat"
LASTNAME = "LhARA:1:Capture:Drift:4"     #.. End of the capture section
GL1NAME  = "LhARA:1:Capture:Gabor lens:1"
BANDS    = [5., 10., 15., 20., 30.]      #.. Energy bands (+/- 1 MeV), MeV

#--------  Simulation:
def simulate(nEvts, seed, datafile):
    import Simulation as Simu
    beamlinefile = os.path.join(os.getenv('LhARAOpticsPATH'), BEAMLINE)
    Smltn = Simu.Simulation(nEvts, beamlinefile, None, \
                            os.path.abspath(datafile))
    #.. Source sampling uses the "random" module; seed it after the
    #   Simulation instance (which seeds with the time) has been created:
    random.seed(seed)
    Smltn.RunSim()

#--------  Reading:
def readData(datafile):
    ibmIOr = bmIO.BeamIO(None, datafile)
    Names, sOf, nReach = None, {}, None
    K, TS, TSend, Kend, TSlast, iLast = [], [], [], [], [], []
    while True:
        EndOfFile = ibmIOr.readBeamDataRecord()
        if EndOfFile:
            break
        iPrtcl = Prtcl.Particle.getinstances()[-1]
        if isinstance(iPrtcl, Prtcl.ReferenceParticle):
            continue
        if Names is None:
            Elmnts  = BLE.BeamLineElement.getinstances()
            AllNames = [e.getName() for e in Elmnts[1:]]
            Names   = AllNames[:AllNames.index(LASTNAME)+1]
            nReach  = np.zeros(len(AllNames))
            iRef    = BL.BeamLine.getcurrentReferenceParticle()
            p0      = np.linalg.norm(iRef.getPrIn()[0][:3])
            E0      = iRef.getPrIn()[0][3]
            mass    = PhysCnst.PhysicalConstants().getparticleMASS( \
                                                         iRef.getSpecies())
        Locs = iPrtcl.getLocation()
        Trc  = iPrtcl.getTraceSpace()
        sLst = iPrtcl.gets()
        nReach[:len(Locs)] += 1
        for iL, Name in enumerate(Locs):
            if Name not in sOf:
                sOf[Name] = sLst[iL]
        Row = np.full((len(Names), 6), np.nan)
        for iL, Name in enumerate(Locs[:len(Names)]):
            Row[Names.index(Name)] = Trc[iL]
        TS.append(Row)
        K.append(E0 + Row[0][5]*p0 - mass)
        if len(Locs) == len(AllNames):
            TSend.append(np.array(Trc[-1]))
            Kend.append(K[-1])
        TSlast.append(np.array(Trc[-1]))
        iLast.append(len(Locs) - 1)
        Prtcl.Particle.cleanParticles()
    return {"Names": Names, "AllNames": AllNames, "sOf": sOf, \
            "nReach": nReach, "K": np.array(K), "TS": np.array(TS), \
            "TSend": np.array(TSend), "Kend": np.array(Kend), \
            "TSlast": np.array(TSlast), "iLast": np.array(iLast), \
            "p0": p0, "E0": E0, "mass": mass}

def classifyLosses(D, b0):
    #.. Reason why each lost particle was not transported through the next
    #   element, following the checks in BeamLineElement.Transport:
    VCMVr  = BLE.Facility.getinstances().getVCMVr()
    nLoc   = len(D["AllNames"])
    Counts = {}
    Lost   = D["iLast"] < nLoc - 1
    for TSl, iL, Kp in zip(D["TSlast"][Lost], D["iLast"][Lost], D["K"][Lost]):
        Next = D["AllNames"][iL+1]
        r    = math.sqrt(TSl[0]**2 + TSl[2]**2)
        D2   = 1. + 2.*TSl[5]/b0 + TSl[5]**2
        eps  = (TSl[1]**2 + TSl[3]**2) / (2.*D2) if D2 > 0. else 2.
        if "Aperture" in Next:
            Why = "aperture"
        elif r >= VCMVr:
            Why = "vacuum chamber"
        elif eps > 1.:
            Why = "expansion parameter (angle)"
        elif abs(TSl[4]) > 5.:
            Why = "longitudinal cut |z| > 5 m"
        else:
            Why = "other"
        C = Counts.setdefault(Why, {"N": 0, "Kmin": 1.E9, "Kmax": 0., \
                                    "smin": 1.E9, "smax": 0.})
        C["N"]   += 1
        C["Kmin"] = min(C["Kmin"], Kp)
        C["Kmax"] = max(C["Kmax"], Kp)
        s         = D["sOf"].get(D["AllNames"][iL], 0.)
        C["smin"] = min(C["smin"], s)
        C["smax"] = max(C["smax"], s)
    return Counts

#--------  Utilities:
def setLength(iBLE, Length):
    if hasattr(iBLE, "setLength"):
        iBLE.setLength(float(Length))
    else:
        iBLE._Length = float(Length)

def lensMatrices(iBLE, Lengths, deltas):
    #.. 4x4 transverse transfer matrices of a Gabor lens, evaluated by the
    #   code, for each length and each delta = (E - E0)/p0:
    L0  = iBLE.getLength()
    Out = np.empty((len(Lengths), len(deltas), 4, 4))
    for iL, Length in enumerate(Lengths):
        setLength(iBLE, Length)
        for iD, delta in enumerate(deltas):
            iBLE.setTransferMatrix(np.array([0., 0., 0., 0., 0., delta]))
            Out[iL, iD] = iBLE.getTransferMatrix()[0:4, 0:4]
    setLength(iBLE, L0)
    return Out

def applyMatrices(Mgrid, deltas, TSin):
    #.. Interpolate the matrices in delta for each particle and apply them:
    d   = TSin[:, 5]
    M   = np.empty((len(d), 4, 4))
    for i in range(4):
        for j in range(4):
            M[:, i, j] = np.interp(d, deltas, Mgrid[:, i, j])
    Out = TSin.copy()
    Out[:, 0:4] = np.einsum('nij,nj->ni', M, TSin[:, 0:4])
    return Out

def beamParameters(TS):
    #.. RMS sizes, divergences, emittances and Twiss parameters:
    Res = {}
    for u, i in (("x", 0), ("y", 2)):
        C   = np.cov(TS[:, i], TS[:, i+1], bias=True)
        eps = math.sqrt(max(np.linalg.det(C), 0.))
        Res["sig" + u]    = math.sqrt(C[0, 0])
        Res["sig" + u + "p"] = math.sqrt(C[1, 1])
        Res["eps" + u]    = eps
        Res["alpha" + u]  = -C[0, 1]/eps if eps > 0. else float('nan')
        Res["beta" + u]   =  C[0, 0]/eps if eps > 0. else float('nan')
    return Res

def thinElement(Name):
    return ("Aperture" in Name) or ("Source" in Name)

#--------  Analysis:
def analyse(datafile, outdir):
    figdir = os.path.join(outdir, "figures")
    resdir = os.path.join(outdir, "results")
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(resdir, exist_ok=True)

    D     = readData(datafile)
    Names = D["Names"]
    TS    = D["TS"]
    K     = D["K"]
    p0, E0, mass = D["p0"], D["E0"], D["mass"]
    N0    = len(K)
    sOf   = D["sOf"]
    Summary = {"nEvents": int(N0), "datafile": datafile}

    Elmnts = {e.getName(): e for e in BLE.BeamLineElement.getinstances()}
    iSrc   = [e for e in Elmnts.values() if isinstance(e, BLE.Source)][0]
    iGL1   = Elmnts[GL1NAME]
    iRef   = BL.BeamLine.getcurrentReferenceParticle()
    pAll   = [np.linalg.norm(P[:3]) for P in iRef.getPrOut()]
    Summary["p0_MeV"]          = p0
    Summary["p0_min_max_MeV"]  = [min(pAll), max(pAll)]
    Summary["NOdispersion"]    = bool(BLE.BeamLineElement.getNOdispersion())

    #.. ----> Source parameters (as used by the code):
    P = iSrc.getParameters()
    Thickness, DivAngle = 5.E-6, 25.      #.. Defaults (Source.setDEFAULTparams)
    t0, Kinf = BLE.Source.calculatet0(P[1], P[2], P[3], Thickness, DivAngle)
    Kinf_MeV = Kinf / 1.602176634E-13
    Intensity = P[1]*P[2] / (math.pi * P[3]**2) * 1.E-4
    a0 = (P[0] * 1.602176634E-19) / \
         (2.*math.pi * 9.1093837E-31 * 299792458.**2 * P[3]) * \
         math.sqrt(2.*P[1]*P[2] / (math.pi * 8.8541878128E-12 * 299792458.))
    X  = math.sqrt(P[7] / Kinf_MeV)
    SrcRows = [ \
        ["Wavelength", P[0]*1.E6, "um", "default"], \
        ["Power", P[1], "W", "default (10 J / 40 fs)"], \
        ["Strehl ratio", P[2], "", "default"], \
        ["Focal spot radius r0", P[3]*1.E6, "um", "default"], \
        ["Pulse duration", P[4]*1.E15, "fs", "default"], \
        ["Target thickness", Thickness*1.E6, "um", "default"], \
        ["Electron divergence half angle", DivAngle, "deg", "default"], \
        ["Intensity", Intensity, "W/cm2", "derived"], \
        ["a0", a0, "", "derived"], \
        ["Hot electron temperature Te", P[5], "MeV", "derived"], \
        ["K_infinity", Kinf_MeV, "MeV", "derived"], \
        ["Ballistic time t0", t0*1.E15, "fs", "derived"], \
        ["t_laser/t0", P[4]/t0, "", "derived"], \
        ["X = sqrt(Kmax/K_infinity)", X, "", "derived"], \
        ["Kmin", P[6], "MeV", "default"], \
        ["Kmax", P[7], "MeV", "derived"], \
        ["SigmaThetaS0", P[8], "deg", "default"], \
        ["SlopeThetaS", P[9], "deg", "default"], \
        ["rpmax", P[10], "", "beam-line file"]]
    with open(os.path.join(resdir, "source_parameters.csv"), "w", \
              newline="") as File:
        W = csv.writer(File)
        W.writerow(["Parameter", "Value", "Unit", "Origin"])
        W.writerows(SrcRows)
    Summary["source"] = {r[0]: r[1] for r in SrcRows}

    #.. ----> Source distributions:
    TS0 = TS[:, 0, :]
    pK  = np.sqrt((K + mass)**2 - mass**2)
    xp  = TS0[:, 1] * p0 / pK              #.. Direction slopes at source
    yp  = TS0[:, 3] * p0 / pK
    rp  = np.sqrt(xp**2 + yp**2)
    Summary["source_K_mean_MeV"]   = float(np.mean(K))
    Summary["source_K_median_MeV"] = float(np.median(K))
    Summary["source_K_min_max"]    = [float(K.min()), float(K.max())]
    Summary["source_frac_K_14_16"] = float(np.mean((K > 14.) & (K < 16.)))
    Summary["source_rp_max_mrad"]  = float(rp.max()*1.E3)
    Summary["source_sig_xp_mrad"]  = float(np.std(xp)*1.E3)
    Summary["source_beam"]         = beamParameters(TS0)

    #.. TNSA model (same functional form as the sampler):
    Kgrid = np.linspace(P[6], P[7], 400)
    Model = np.exp(-np.sqrt(2.*Kgrid/P[5])) / np.sqrt(Kgrid)
    Model /= np.trapezoid(Model, Kgrid)
    thetaMax = P[8] - P[9]*Kgrid/P[7]                 #.. Source.g_theta
    upmax    = np.sin(np.radians(thetaMax))
    a        = P[10] / upmax
    AccFrac  = 2.*a**2 - a**4                         #.. Fraction r' < rpmax
    Summary["accepted_fraction_at_1_15_Kmax"] = [ \
        float(np.interp(v, Kgrid, AccFrac)) for v in (P[6], 15., P[7])]
    Summary["accepted_fraction_spectrum_weighted"] = \
        float(np.trapezoid(AccFrac*Model, Kgrid))

    fig, axs = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    Bins = np.linspace(P[6], P[7], 50)
    axs[0].hist(K, bins=Bins, density=True, histtype="stepfilled", \
                alpha=0.5, label="Simulated (%d protons)" % N0)
    axs[0].plot(Kgrid, Model, "k-", lw=1.5, \
                label=r"Model $\propto e^{-\sqrt{2K/T_e}}/\sqrt{K}$")
    axs[0].axvline(15., color="r", ls="--", lw=1, label="Reference 15 MeV")
    axs[0].set_yscale("log")
    axs[0].set_xlabel("Kinetic energy K (MeV)")
    axs[0].set_ylabel("Probability density (1/MeV)")
    axs[0].set_title(r"Source energy spectrum ($T_e$=%.2f MeV, " \
                     r"$K_{max}$=%.2f MeV)" % (P[5], P[7]))
    axs[0].legend(fontsize=8)
    WModel = AccFrac*Model / np.trapezoid(AccFrac*Model, Kgrid)
    axs[1].plot(Kgrid, Model, "k-", lw=1.5, \
                label="Sampled spectrum (all angles)")
    axs[1].plot(Kgrid, WModel, "b--", lw=1.5, \
                label=r"Spectrum within $r'<$%.0f mrad" % (P[10]*1.E3))
    axs[1].set_yscale("log")
    axs[1].set_xlabel("Kinetic energy K (MeV)")
    axs[1].set_ylabel("Probability density (1/MeV)")
    axs[1].set_title("Effect of the angular cut on the spectrum (model)")
    axs[1].legend(fontsize=8)
    fig.savefig(os.path.join(figdir, "fig1_source_spectrum.png"), dpi=150)
    plt.close(fig)

    Sub = np.arange(min(N0, 5000))
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), layout="constrained")
    h, e, _ = axs[0].hist(rp*1.E3, bins=40, range=(0., P[10]*1.E3), \
                          histtype="stepfilled", alpha=0.5, label="Simulated")
    c = 0.5*(e[1:] + e[:-1])
    axs[0].plot(c, c * N0*(e[1]-e[0]) * 2./(P[10]*1.E3)**2, "k-", \
                label=r"Uniform disc ($\propto r'$)")
    axs[0].set_xlabel(r"$r' = \sqrt{x_p^2+y_p^2}$ at source (mrad)")
    axs[0].set_ylabel("Protons per bin")
    axs[0].set_title("Emission angle at the source")
    axs[0].legend(fontsize=8)
    axs[1].scatter(xp[Sub]*1.E3, yp[Sub]*1.E3, s=1, alpha=0.5)
    axs[1].set_aspect("equal")
    axs[1].set_xlabel(r"$x_p$ (mrad)")
    axs[1].set_ylabel(r"$y_p$ (mrad)")
    axs[1].set_title("Emission directions (5000 protons)")
    axs[2].plot(Kgrid, thetaMax, "k-", label=r"$\theta_{max}(K)$ of model")
    axs[2].axhline(np.degrees(P[10]), color="r", ls="--", \
                   label=r"Cut $r'<$%.0f mrad (%.2f$^\circ$)" % \
                   (P[10]*1.E3, np.degrees(P[10])))
    axs[2].set_yscale("log")
    axs[2].set_xlabel("Kinetic energy K (MeV)")
    axs[2].set_ylabel("Angle (degrees)")
    axs[2].set_title("Model emission cone versus angular cut")
    axs[2].legend(fontsize=8)
    fig.savefig(os.path.join(figdir, "fig2_source_angles.png"), dpi=150)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), layout="constrained")
    axs[0].hist(TS0[:, 0]*1.E6, bins=60, density=True, \
                histtype="stepfilled", alpha=0.5, label="Simulated x")
    xg = np.linspace(-6., 6., 200)
    axs[0].plot(xg, np.exp(-0.5*(xg/(P[3]*1.E6))**2) / \
                (math.sqrt(2.*math.pi)*P[3]*1.E6), "k-", \
                label=r"Gaussian, $\sigma$ = $r_0$ = %.1f $\mu$m" % (P[3]*1.E6))
    axs[0].set_xlabel(r"x at source ($\mu$m)")
    axs[0].set_ylabel(r"Probability density (1/$\mu$m)")
    axs[0].set_title("Transverse position at the source")
    axs[0].legend(fontsize=8)
    sc = axs[1].scatter(TS0[Sub, 0]*1.E6, TS0[Sub, 1]*1.E3, c=K[Sub], \
                        s=1, cmap="viridis")
    axs[1].set_xlabel(r"x ($\mu$m)")
    axs[1].set_ylabel(r"x' = $p_x/p_0$ (mrad)")
    axs[1].set_title("Trace space (x, x') at the source")
    fig.colorbar(sc, ax=axs[1], label="K (MeV)")
    axs[2].scatter(K[Sub], np.abs(TS0[Sub, 1])*1.E3, s=1, alpha=0.4)
    axs[2].plot(Kgrid, P[10]*1.E3*np.sqrt((Kgrid+mass)**2-mass**2)/p0, \
                "r-", label=r"$r'_{max}\,p/p_0$")
    axs[2].set_xlabel("Kinetic energy K (MeV)")
    axs[2].set_ylabel(r"|x'| at source (mrad)")
    axs[2].set_title(r"x' = $x_p\,p/p_0$ scales with momentum")
    axs[2].legend(fontsize=8)
    fig.savefig(os.path.join(figdir, "fig3_source_transverse.png"), dpi=150)
    plt.close(fig)

    #.. ----> Beam parameters at the exit of each element of the section:
    Rows = []
    for iN, Name in enumerate(Names):
        Ok  = ~np.isnan(TS[:, iN, 0])
        BP  = beamParameters(TS[Ok, iN, :])
        #.. Transverse positions rescaled by p0/p (see report):
        TSc = TS[Ok, iN, :].copy()
        TSc[:, 0] *= p0/pK[Ok]
        TSc[:, 2] *= p0/pK[Ok]
        BPc = beamParameters(TSc)
        Rows.append([Name, sOf[Name], int(Ok.sum()), Ok.sum()/N0, \
                     BP["sigx"]*1.E3, BP["sigy"]*1.E3, \
                     BP["sigxp"]*1.E3, BP["sigyp"]*1.E3, \
                     BP["epsx"]*1.E6, BP["epsy"]*1.E6, \
                     BP["alphax"], BP["betax"], BP["alphay"], BP["betay"], \
                     BPc["sigx"]*1.E3, BPc["sigy"]*1.E3, BPc["epsx"]*1.E6])
    Header = ["Location", "s_m", "N", "Transmission", "sigx_mm", "sigy_mm", \
              "sigxp_mrad", "sigyp_mrad", "epsx_mm_mrad", "epsy_mm_mrad", \
              "alphax", "betax_m", "alphay", "betay_m", \
              "sigx_corr_mm", "sigy_corr_mm", "epsx_corr_mm_mrad"]
    with open(os.path.join(resdir, "beam_parameters_capture.csv"), "w", \
              newline="") as File:
        W = csv.writer(File)
        W.writerow(Header)
        W.writerows(Rows)
    Summary["capture_table"] = [dict(zip(Header, r)) for r in Rows]

    #.. ----> Envelope inside the elements, using the element matrices:
    Kd     = np.linspace(P[6]*0.99, P[7]*1.01, 600)
    deltas = (Kd + mass - E0) / p0
    nSub   = 20
    sCurve, TSCurve = [0.], [TS[:, 0, :]]
    MaxDev = {}
    for iN in range(1, len(Names)):
        Name  = Names[iN]
        iBLE  = Elmnts[Name]
        Lngth = iBLE.getLength() if hasattr(iBLE, "getLength") else 0.
        TSin  = TS[:, iN-1, :]
        if thinElement(Name) or Lngth <= 0.:
            continue
        fr = np.linspace(0., 1., nSub+1)[1:]
        if isinstance(iBLE, BLE.GaborLens):
            Mg = lensMatrices(iBLE, fr*Lngth, deltas)
        for iF, f in enumerate(fr):
            if isinstance(iBLE, BLE.GaborLens):
                TSf = applyMatrices(Mg[iF], deltas, np.nan_to_num(TSin))
            else:
                TSf = TSin.copy()
                TSf[:, 0] += f*Lngth*TSin[:, 1]
                TSf[:, 2] += f*Lngth*TSin[:, 3]
            sCurve.append(sOf[Names[iN-1]] + f*Lngth)
            TSCurve.append(TSf)
        Ok = ~np.isnan(TS[:, iN, 0])
        Dev = np.nanmax(np.abs(TSCurve[-1][Ok, 0:4] - TS[Ok, iN, 0:4]), \
                        axis=0) / np.nanstd(TS[Ok, iN, 0:4], axis=0)
        MaxDev[Name] = float(Dev.max())
    Summary["envelope_check_max_relative_deviation"] = MaxDev
    sCurve = np.array(sCurve)

    def rms(TSs, col, Sel):
        return np.array([np.std(T[Sel, col]) for T in TSs])
    AllSel = np.ones(N0, dtype=bool)
    fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True, \
                            layout="constrained")
    for ax, col, lab in ((axs[0], 0, r"$\sigma_x$ (mm)"), \
                         (axs[1], 1, r"$\sigma_{x'}$ (mrad)")):
        ax.plot(sCurve, rms(TSCurve, col, AllSel)*1.E3, "k-", lw=2.2, \
                label="All protons")
        if col == 0:
            ax.plot(sCurve, np.array([np.std(T[:, 0]*p0/pK) \
                    for T in TSCurve])*1.E3, "k--", lw=1.5, \
                    label=r"All protons, x rescaled by $p_0/p$")
        for Kb in BANDS:
            Sel = np.abs(K - Kb) < 1.
            if Sel.sum() > 50:
                ax.plot(sCurve, rms(TSCurve, col, Sel)*1.E3, lw=1.1, \
                        label=r"K = %.0f $\pm$ 1 MeV" % Kb)
        for Name in Names:
            iBLE = Elmnts[Name]
            if isinstance(iBLE, BLE.GaborLens):
                s1 = sOf[Name]
                ax.axvspan(s1 - iBLE.getLength(), s1, color="orange", \
                           alpha=0.15)
            if "Aperture" in Name:
                ax.axvline(sOf[Name], color="grey", ls=":", lw=1)
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    axs[0].set_title("RMS envelope from the source to the end of the " \
                     "capture section (shaded: Gabor lenses; dotted: " \
                     "apertures)", fontsize=10)
    axs[0].legend(fontsize=8, ncol=2)
    axs[1].set_xlabel("s (m)")
    fig.savefig(os.path.join(figdir, "fig4_envelope_capture.png"), dpi=150)
    plt.close(fig)

    #.. ----> Trace space at entrance and exit of Gabor lens 1:
    iEnt = Names.index(GL1NAME) - 1
    iExt = Names.index(GL1NAME)
    Spots = [(iEnt, "Entrance of Gabor lens 1"), \
             (iExt, "Exit of Gabor lens 1"), \
             (iExt+1, "0.15 m after Gabor lens 1")]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3), layout="constrained")
    for ax, (iN, Ttl) in zip(axs, Spots):
        sc = ax.scatter(TS[Sub, iN, 0]*1.E3, TS[Sub, iN, 1]*1.E3, c=K[Sub], \
                        s=1, cmap="viridis", vmin=P[6], vmax=P[7])
        ax.set_title("%s (s = %.3f m)" % (Ttl, sOf[Names[iN]]), fontsize=10)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("x' (mrad)")
        ax.grid(alpha=0.3)
    fig.colorbar(sc, ax=axs, label="K (MeV)")
    fig.savefig(os.path.join(figdir, "fig5_trace_space_GL1.png"), dpi=150)
    plt.close(fig)

    #.. Emittance in narrow energy bands (chromatic emittance growth):
    BandEps = []
    for Kb in BANDS:
        Sel = np.abs(K - Kb) < 0.25
        if Sel.sum() > 20:
            BandEps.append([Kb, int(Sel.sum()), \
                beamParameters(TS[Sel, iEnt, :])["epsx"]*1.E6, \
                beamParameters(TS[Sel, iExt, :])["epsx"]*1.E6])
    Summary["band_emittance_GL1_entrance_exit_mm_mrad"] = BandEps

    #.. ----> Chromatic focusing of Gabor lens 1:
    LGL1 = iGL1.getLength()
    Mfull = lensMatrices(iGL1, [LGL1], deltas)[0]
    d     = sOf[Names[iEnt]]                  #.. Source to lens entrance
    M11, M12 = Mfull[:, 0, 0], Mfull[:, 0, 1]
    M21, M22 = Mfull[:, 1, 0], Mfull[:, 1, 1]
    Ratio = (M21*d + M22) / (M11*d + M12)      #.. x'/x at exit, point source
    sqrtk = np.arccos(np.clip(M11, -1., 1.)) / LGL1
    focal = -1. / M21
    i15   = np.argmin(np.abs(Kd - 15.))
    def zeros(F):
        return [float(Kd[z] - F[z]*(Kd[z+1]-Kd[z])/(F[z+1]-F[z])) \
                for z in np.where(np.diff(np.sign(F)) != 0)[0]]
    Kpar  = zeros(M21*d + M22)      #.. x' = 0 at exit: parallel beam
    Kimg  = zeros(M11*d + M12)      #.. x  = 0 at exit: imaged onto axis
    Summary["GL1"] = { \
        "Length_m": LGL1, "Strength_internal": iGL1.getStrength(), \
        "ElectronDensity_m-3": iGL1.getElectronDensity(), \
        "d_source_to_entrance_m": d, \
        "k_15MeV_m-2": float(sqrtk[i15]**2), \
        "phase_advance_15MeV_rad": float(sqrtk[i15]*LGL1), \
        "focal_length_15MeV_m": float(focal[i15]), \
        "xp_over_x_exit_15MeV_per_m": float(Ratio[i15]), \
        "K_parallel_after_GL1_MeV": Kpar, \
        "K_imaged_on_axis_at_exit_MeV": Kimg, \
        "M_15MeV": Mfull[i15, 0:2, 0:2].tolist(), \
        "max_xy_coupling": float(np.abs(Mfull[:, 0:2, 2:4]).max())}
    with open(os.path.join(resdir, "GL1_chromatic.csv"), "w", \
              newline="") as File:
        W = csv.writer(File)
        W.writerow(["K_MeV", "delta", "M11", "M12_m", "M21_per_m", "M22", \
                    "sqrt_k_per_m", "phase_rad", "focal_length_m", \
                    "xp_over_x_exit_per_m"])
        for i in range(0, len(Kd), 10):
            W.writerow([Kd[i], deltas[i], M11[i], M12[i], M21[i], M22[i], \
                        sqrtk[i], sqrtk[i]*LGL1, focal[i], Ratio[i]])

    xE, xpE = TS[:, iExt, 0], TS[:, iExt, 1]
    Big = np.abs(xE) > 0.05*np.nanstd(xE)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3), layout="constrained")
    axs[0].scatter(K[Big], xpE[Big]/xE[Big], s=1, alpha=0.3, \
                   label="Simulated protons")
    axs[0].plot(Kd, Ratio, "r-", lw=1.5, label="Gabor-lens transfer matrix")
    axs[0].axhline(0., color="k", lw=0.8)
    axs[0].axvline(15., color="grey", ls="--", lw=1)
    axs[0].set_ylim(-6., 2.)
    axs[0].set_xlabel("Kinetic energy K (MeV)")
    axs[0].set_ylabel(r"x'/x at exit of Gabor lens 1 (m$^{-1}$)")
    axs[0].set_title(r"> 0 diverging, = 0 parallel, < 0 converging", \
                     fontsize=10)
    axs[0].legend(fontsize=8)
    axs[1].plot(Kd, focal, "b-", label="Focal length f = -1/M$_{21}$ (m)")
    axs[1].plot(Kd, sqrtk*LGL1, "g-", \
                label=r"Phase advance $\sqrt{k}L$ (rad)")
    axs[1].axvline(15., color="grey", ls="--", lw=1)
    axs[1].set_xlabel("Kinetic energy K (MeV)")
    axs[1].set_title("Gabor lens 1 focusing versus energy (f diverges " \
                     r"where $\sqrt{k}L = \pi$)", fontsize=10)
    axs[1].set_ylim(-1., 4.)
    axs[1].legend(fontsize=8)
    axs[1].grid(alpha=0.3)
    fig.savefig(os.path.join(figdir, "fig6_chromatic_GL1.png"), dpi=150)
    plt.close(fig)

    #.. ----> Whole beam line: transmission and spectrum at the end:
    AllNames = D["AllNames"]
    nReach   = D["nReach"]
    with open(os.path.join(resdir, "transmission_full_line.csv"), "w", \
              newline="") as File:
        W = csv.writer(File)
        W.writerow(["Index", "Location", "s_m", "N", "Fraction"])
        for i, Name in enumerate(AllNames):
            W.writerow([i, Name, sOf.get(Name, ""), int(nReach[i]), \
                        nReach[i]/N0])
    Summary["N_end"]          = int(nReach[-1])
    Summary["fraction_end"]   = float(nReach[-1]/N0)
    Summary["first_big_loss"] = []
    for i in range(1, len(AllNames)):
        Lost = nReach[i-1] - nReach[i]
        if Lost > 0.05*N0:
            Summary["first_big_loss"].append([AllNames[i], int(Lost)])
    Summary["loss_classification"] = classifyLosses(D, p0/E0)
    with open(os.path.join(resdir, "loss_classification.csv"), "w", \
              newline="") as File:
        W = csv.writer(File)
        W.writerow(["Reason", "N", "Fraction_of_generated", "Kmin_MeV", \
                    "Kmax_MeV", "s_first_m", "s_last_m"])
        for Why, C in Summary["loss_classification"].items():
            W.writerow([Why, C["N"], C["N"]/N0, C["Kmin"], C["Kmax"], \
                        C["smin"], C["smax"]])
    Kend = D["Kend"]
    if len(Kend) > 1:
        Summary["end_K_mean_rms_MeV"] = [float(np.mean(Kend)), \
                                         float(np.std(Kend))]
        Summary["end_beam"] = beamParameters(D["TSend"])

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.3), layout="constrained")
    sPlot = [sOf.get(Name, np.nan) for Name in AllNames]
    axs[0].step(sPlot, nReach/N0, where="post", color="k")
    axs[0].set_yscale("log")
    axs[0].set_xlabel("s (m)")
    axs[0].set_ylabel("Fraction of generated protons")
    axs[0].set_title("Transmission along the whole Stage 1 beam line", \
                     fontsize=10)
    axs[0].grid(alpha=0.3, which="both")
    axs[1].hist(K, bins=Bins, histtype="step", color="k", lw=1.5, \
                label="Source (%d)" % N0)
    if len(Kend) > 0:
        axs[1].hist(Kend, bins=Bins, histtype="stepfilled", alpha=0.6, \
                    label="End of beam line (%d)" % len(Kend))
    axs[1].set_yscale("log")
    axs[1].set_xlabel("Kinetic energy K (MeV)")
    axs[1].set_ylabel("Protons per bin")
    axs[1].set_title("Energy spectrum: source and end of beam line", \
                     fontsize=10)
    axs[1].legend(fontsize=8)
    fig.savefig(os.path.join(figdir, "fig7_full_line.png"), dpi=150)
    plt.close(fig)

    Ok = ~np.isnan(TS[:, iExt, 0])
    electronAnalysis(P, K, Ok, TS[:, iExt, 0], Kend, Summary, figdir, resdir)

    with open(os.path.join(resdir, "summary.json"), "w") as File:
        json.dump(Summary, File, indent=1, default=float)
    print(json.dumps(Summary, indent=1, default=float))

#--------  Hot electrons at the source:
QE   = 1.602176634E-19
ME   = 9.1093837E-31
MP   = 1.67262192E-27
CL   = 299792458.
EPS0 = 8.8541878128E-12
MEC2 = 0.51099895                       #.. MeV

def KmaxModel(P, Thickness, DivAngle):
    #.. Cut-off energy as in BeamLine.parseSource (Source.calculatet0 and
    #   Source.DurationBYt0equation):
    t0, Kinf = BLE.Source.calculatet0(P[1], P[2], P[3], Thickness, DivAngle)
    X = float(fsolve(BLE.Source.DurationBYt0equation, 0.5, \
                     args=(P[4], t0))[0])
    return Kinf*X**2 / (QE*1.E6), t0

def tnsaPDF(Kv, Te, Kmin, Kmax):
    Kg  = np.linspace(Kmin, Kmax, 4000)
    Nrm = np.trapezoid(np.exp(-np.sqrt(2.*Kg/Te))/np.sqrt(Kg), Kg)
    Kc  = np.clip(Kv, Kmin, None)
    f   = np.exp(-np.sqrt(2.*Kc/Te))/np.sqrt(Kc) / Nrm
    return np.where((Kv >= Kmin) & (Kv <= Kmax), f, 0.)

def electronAnalysis(P, K, Ok, xGL1, Kend, Summary, figdir, resdir):
    Thickness, DivAngle = 5.E-6, 25.         #.. Defaults
    #.. Hot-electron temperature (as Source.calculateTe):
    Intensity = P[1]*P[2] / (math.pi*P[3]**2) * 1.E-4      #.. W/cm2
    a0   = (P[0]*QE) / (2.*math.pi*ME*CL**2*P[3]) * \
           math.sqrt(2.*P[1]*P[2] / (math.pi*EPS0*CL))
    y0   = a0*P[0] / (2.*math.pi)
    rLe  = P[3] / (2.*math.log(2.))
    Corr = math.sqrt(1. - (1. - rLe/y0)**2) if rLe < y0 else 1.
    TePond = MEC2*(math.sqrt(1. + a0**2) - 1.)
    TeCorr = MEC2*(math.sqrt(1. + (a0*Corr)**2) - 1.)
    #.. Electron number, density and sheath (Fuchs et al. model, TN §5):
    etaRaw = 1.2E-15 * Intensity**0.75
    eta    = min(etaRaw, 0.5)
    Elaser = P[1]*P[4]*P[2]                  #.. J, within the Strehl spot
    TeJ    = P[5]*QE*1.E6
    Ne     = eta*Elaser/TeJ
    Bsh    = P[3] + Thickness*math.tan(math.radians(DivAngle))
    Ssh    = math.pi*Bsh**2
    ne0    = Ne / (CL*P[4]*Ssh)
    ncrit  = EPS0*ME*(2.*math.pi*CL/P[0])**2 / QE**2
    E0     = math.sqrt(ne0*TeJ/EPS0)
    lamD   = math.sqrt(EPS0*TeJ/(ne0*QE**2))
    cs     = math.sqrt(TeJ/MP)
    wpi    = math.sqrt(ne0*QE**2/(MP*EPS0))
    def Mora(t):
        tau = wpi*t / math.sqrt(2.*math.e)
        return 2.*P[5]*math.log(tau + math.sqrt(tau**2 + 1.))**2
    KmaxCode, t0 = KmaxModel(P, Thickness, DivAngle)
    Rows = [ \
        ["Intensity (Strehl weighted)", Intensity, "W/cm2", ""], \
        ["Normalised vector potential a0", a0, "", ""], \
        ["Transverse acceleration length y0", y0*1.E6, "um", ""], \
        ["Effective spot radius r0/(2 ln 2)", rLe*1.E6, "um", ""], \
        ["Small-spot correction pt/pt0", Corr, "", "Dover et al. 2020"], \
        ["Te, ponderomotive (no correction)", TePond, "MeV", ""], \
        ["Te, with small-spot correction", TeCorr, "MeV", \
         "code: %.4f" % P[5]], \
        ["Conversion efficiency, raw scaling", etaRaw, "", \
         "1.2e-15 I^0.75"], \
        ["Conversion efficiency used", eta, "", "capped at 0.5"], \
        ["Laser energy in focal spot", Elaser, "J", "E_laser x Strehl"], \
        ["Number of hot electrons Ne", Ne, "", "eta E / Te"], \
        ["Sheath radius B = r0 + d tan(theta_e)", Bsh*1.E6, "um", ""], \
        ["Sheath area", Ssh, "m2", ""], \
        ["Hot-electron density ne0", ne0, "m-3", ""], \
        ["ne0 / critical density", ne0/ncrit, "", \
         "n_c = %.3e m-3" % ncrit], \
        ["Sheath field scale sqrt(ne0 Te/eps0)", E0, "V/m", ""], \
        ["Debye length", lamD*1.E6, "um", ""], \
        ["Ion-acoustic speed c_s", cs/CL, "c", ""], \
        ["Ion plasma frequency x t_laser", wpi*P[4], "", ""], \
        ["Ballistic time t0", t0*1.E15, "fs", ""], \
        ["Kmax, code model", KmaxCode, "MeV", "code: %.4f" % P[7]], \
        ["Kmax, Mora (t_acc = t_laser)", Mora(P[4]), "MeV", "PRL 90 185002"], \
        ["Kmax, Mora (t_acc = 1.3 t_laser)", Mora(1.3*P[4]), "MeV", ""]]
    with open(os.path.join(resdir, "electron_source_model.csv"), "w", \
              newline="") as File:
        W = csv.writer(File)
        W.writerow(["Quantity", "Value", "Unit", "Note"])
        W.writerows(Rows)
    Summary["electrons"] = {r[0]: r[1] for r in Rows}

    #.. Cut-off energy versus electron divergence and target thickness:
    Angles  = np.arange(5., 50.01, 2.5)
    KmaxTab = {dum: np.array([KmaxModel(P, dum*1.E-6, a)[0] \
                              for a in Angles]) for dum in (2., 5., 10.)}
    i25     = int(np.argmin(np.abs(Angles - DivAngle)))
    dKdth   = (KmaxTab[5.][i25+1] - KmaxTab[5.][i25-1]) / \
              (Angles[i25+1] - Angles[i25-1])
    Summary["dKmax_dthetae_MeV_per_deg"] = float(dKdth)

    #.. Reweighting the simulated protons to other spectra:
    w0  = tnsaPDF(K, P[5], P[6], P[7])
    w0e = tnsaPDF(Kend, P[5], P[6], P[7]) if len(Kend) > 0 else None
    def analytic15(Te, Kmax):
        Kg = np.linspace(14., 16., 200)
        return float(np.trapezoid(tnsaPDF(Kg, Te, P[6], Kmax), Kg))
    def reweigh(Te, Kmax):
        w    = tnsaPDF(K, Te, P[6], Kmax) / w0
        F15  = np.sum(w*((K > 14.) & (K < 16.))) / np.sum(w)
        wx   = w[Ok]
        xm   = np.sum(wx*xGL1[Ok]) / np.sum(wx)
        sigx = math.sqrt(np.sum(wx*(xGL1[Ok]-xm)**2) / np.sum(wx))
        Kmn  = np.sum(w*K) / np.sum(w)
        Fend = float('nan')
        if w0e is not None:
            Fend = np.sum(tnsaPDF(Kend, Te, P[6], Kmax)/w0e) / np.sum(w)
        return [float(F15), sigx, float(Fend), float(Kmn)]
    TeRows = []
    for f in (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3):
        TeRows.append([f, f*P[5]] + reweigh(f*P[5], P[7]))
    with open(os.path.join(resdir, "electron_sensitivity_Te.csv"), "w", \
              newline="") as File:
        W = csv.writer(File)
        W.writerow(["Te_over_Te0", "Te_MeV", "Fraction_14_16MeV", \
                    "sigx_GL1_exit_m", "Fraction_end", "Kmean_MeV"])
        W.writerows(TeRows)
    AngRows = []
    for a, Km in zip(Angles, KmaxTab[5.]):
        Row = [a, Km, analytic15(P[5], Km)]
        Row += reweigh(P[5], Km) if Km <= P[7] else \
               [float('nan')]*4
        AngRows.append(Row)
    with open(os.path.join(resdir, "electron_sensitivity_angle.csv"), \
              "w", newline="") as File:
        W = csv.writer(File)
        W.writerow(["theta_e_deg", "Kmax_MeV", "Fraction_14_16MeV_model", \
                    "Fraction_14_16MeV_reweighted", "sigx_GL1_exit_m", \
                    "Fraction_end", "Kmean_MeV"])
        W.writerows(AngRows)
    Summary["electron_sensitivity_Te"]    = TeRows
    Summary["electron_sensitivity_angle"] = AngRows
    Summary["Kmax_vs_angle_d2_d10um"]     = {"angles": Angles.tolist(), \
        "d2": KmaxTab[2.].tolist(), "d10": KmaxTab[10.].tolist()}

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.4), layout="constrained")
    for dum, st in ((2., "--"), (5., "-"), (10., ":")):
        axs[0].plot(Angles, KmaxTab[dum], "k" + st, \
                    label=r"target thickness %.0f $\mu$m" % dum)
    axs[0].plot([DivAngle], [P[7]], "ro", label="Default (25$^\circ$, " \
                r"5 $\mu$m)")
    axs[0].axhline(15., color="grey", ls="-.", lw=1, label="15 MeV")
    axs[0].set_xlabel(r"Hot-electron divergence half angle $\theta_e$ (deg)")
    axs[0].set_ylabel(r"Proton cut-off energy $K_{max}$ (MeV)")
    axs[0].set_title("Cut-off energy versus electron sheath size", \
                     fontsize=10)
    axs[0].legend(fontsize=8)
    axs[0].grid(alpha=0.3)
    Kg = np.linspace(P[6], P[7], 400)
    for f, st in ((0.8, "--"), (1.0, "-"), (1.2, ":")):
        axs[1].plot(Kg, tnsaPDF(Kg, f*P[5], P[6], P[7]), "k" + st, \
                    label=r"$T_e$ = %.1f MeV" % (f*P[5]))
    axs[1].set_yscale("log")
    axs[1].set_xlabel("Proton kinetic energy K (MeV)")
    axs[1].set_ylabel("Probability density (1/MeV)")
    axs[1].set_title("Proton spectrum versus hot-electron temperature", \
                     fontsize=10)
    axs[1].legend(fontsize=8)
    Tr = np.array(TeRows)
    axs[2].plot(Tr[:, 0], Tr[:, 2]*100., "bo-", \
                label="Protons with 14 < K < 16 MeV (%)")
    axs[2].set_xlabel(r"$T_e / T_{e0}$")
    axs[2].set_ylabel("Fraction (%)", color="b")
    ax2 = axs[2].twinx()
    ax2.plot(Tr[:, 0], Tr[:, 3]*1.E3, "rs-", \
             label=r"$\sigma_x$ at exit of Gabor lens 1")
    ax2.set_ylabel(r"$\sigma_x$ at exit of Gabor lens 1 (mm)", color="r")
    axs[2].set_title("Sensitivity to the electron temperature " \
                     "(reweighted simulation)", fontsize=10)
    axs[2].grid(alpha=0.3)
    fig.savefig(os.path.join(figdir, "fig8_electrons.png"), dpi=150)
    plt.close(fig)

#--------  Main:
def main(argv):
    prsr = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    prsr.add_argument("step", choices=["simulate", "analyse"])
    prsr.add_argument("-n", "--nEvts", type=int, default=50000)
    prsr.add_argument("-s", "--seed", type=int, default=20261026)
    prsr.add_argument("-d", "--datafile", default=DATAFILE)
    prsr.add_argument("-o", "--outputdir", \
                      default=os.path.dirname(os.path.abspath(__file__)))
    args = prsr.parse_args(argv)
    if os.getenv('LhARAOpticsPATH') is None:
        print(" sourceGL1Analysis: source startup.bash first.")
        return 1
    if args.step == "simulate":
        simulate(args.nEvts, args.seed, args.datafile)
    else:
        analyse(args.datafile, args.outputdir)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

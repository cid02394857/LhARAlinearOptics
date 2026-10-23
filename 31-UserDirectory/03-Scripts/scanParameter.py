#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scanParameter.py:
=================

  Example of a parameter scan.  One parameter in a beam-line specification
  (csv) file is set to each of a list of values in turn.  For each value
  the beam is simulated (03-Scripts/runBeamSim.py), the beam envelope is
  evaluated (03-Scripts/plotBeam.py) and the number of particles and the
  RMS beam size at the end of the beam line are tabulated.

  Usage:
    scanParameter.py -b <beamlinefile> -e <Element> -p <Parameter>
                     -v <value1,value2,...> [-k <occurrence>] [-n <nEvts>]
                     [-i <inputfile> [-z]]

    -e, -p : "Element" and "Parameter" columns of the line to change; if
             several lines match, -k selects which (1 = first, default).
    -i, -z : optional input particle file; -z if it is in BDSIM format.

  Example: scan the strength of the final-focus F quadrupole of the LhARA
  Stage 1 beam delivery, reading particles from a BDSIM-format file:
    scanParameter.py -b 11-Parameters/LhARA-Stage1-BeamDelivery.csv \\
                     -e Fquad -p kq -v 18,19,20,21,22,23,24 \\
                     -i 12-Data4Tests/BDSIM-Data4Tests.dat -z

  Files are looked for in the working directory, then $HOMEPATH, then
  $LhARAOpticsPATH.  Output is written to $REPORTPATH (99-Scratch):
    scan-<Element>-<Parameter>.csv : value, particles at end, sigma_x, sigma_y
    scan-<Element>-<Parameter>.pdf : the same, plotted versus value

  Requires the environment set by startup.bash (or user_startup.bash).

Created on Fri 23Oct26: Version history:
----------------------------------------
 1.0: 23Oct26: First implementation

"""

import os
import sys
import ast
import csv
import argparse
import subprocess

import matplotlib.pyplot as plt

def findFile(fileName):
    #.. Working directory, then $HOMEPATH, then the repository:
    for Dir in ["", os.getenv('HOMEPATH', ""), os.getenv('LhARAOpticsPATH', "")]:
        Path = os.path.join(Dir, fileName)
        if os.path.isfile(Path):
            return os.path.abspath(Path)
    return None

def writeBeamLine(srcFile, dstFile, Element, Parameter, Occurrence, Value):
    with open(srcFile, newline='', encoding='utf-8-sig') as File:
        Rows = list(csv.reader(File))
    iElm = Rows[0].index("Element")
    iPrm = Rows[0].index("Parameter")
    iVal = Rows[0].index("Value")

    nFound = 0
    for Row in Rows[1:]:
        if len(Row) > iVal and Row[iElm] == Element and Row[iPrm] == Parameter:
            nFound += 1
            if nFound == Occurrence:
                Row[iVal] = str(Value)
                break
    if nFound < Occurrence:
        return False

    with open(dstFile, 'w', newline='', encoding='utf-8-sig') as File:
        csv.writer(File).writerows(Rows)
    return True

def runScript(Script, Args, logFile):
    Cmd = [sys.executable, \
           os.path.join(os.getenv('LhARAOpticsPATH'), '03-Scripts', Script)] \
          + Args
    with open(logFile, 'w') as Log:
        rc = subprocess.run(Cmd, stdout=Log, stderr=subprocess.STDOUT, \
                            stdin=subprocess.DEVNULL).returncode
    return rc == 0

def readEndOfBeamLine(beamCSV):
    #.. Last line of the summary written by plotBeam.py:
    with open(beamCSV, newline='') as File:
        Last = list(csv.DictReader(File))[-1]
    nPrtcls = float(Last["nPrtcls"]) if Last["nPrtcls"] not in ("", "None") \
              else 0.
    sigmaxy = [None, None]
    if Last["sigmaxy"] not in ("", "None"):
        sigmaxy = ast.literal_eval(Last["sigmaxy"])
    return Last["Name"], nPrtcls, sigmaxy[0], sigmaxy[1]

def main(argv):
    prsr = argparse.ArgumentParser( \
                  description="Scan one parameter of a beam-line file.")
    prsr.add_argument("-b", "--beamlinefile", required=True, \
                      help="beam-line specification file (.csv)")
    prsr.add_argument("-e", "--element", required=True, \
                      help='"Element" column, e.g. Fquad')
    prsr.add_argument("-p", "--parameter", required=True, \
                      help='"Parameter" column, e.g. kq')
    prsr.add_argument("-v", "--values", required=True, \
                      help="comma-separated list of values")
    prsr.add_argument("-k", "--occurrence", type=int, default=1, \
                      help="which of the matching lines to change (1=first)")
    prsr.add_argument("-n", "--nEvts", type=int, default=1000, \
                      help="number of events for each value")
    prsr.add_argument("-i", "--inputfile", default=None, \
                      help="input particle file (.dat)")
    prsr.add_argument("-z", "--bdsimfile", action='store_true', \
                      help="input file is in BDSIM format")
    args = prsr.parse_args(argv)

    if os.getenv('LhARAOpticsPATH') is None or os.getenv('REPORTPATH') is None:
        print(" scanParameter: environment not set; source startup.bash first.")
        return 1

    beamlinefile = findFile(args.beamlinefile)
    if beamlinefile is None:
        print(" scanParameter: beam-line file not found:", args.beamlinefile)
        return 1
    inputfile = None
    if args.inputfile is not None:
        inputfile = findFile(args.inputfile)
        if inputfile is None:
            print(" scanParameter: input file not found:", args.inputfile)
            return 1
    Values = [float(Value) for Value in args.values.split(",")]

    Scratch = os.getenv('REPORTPATH')
    os.makedirs(Scratch, exist_ok=True)
    Tag     = "scan-" + args.element.replace(" ", "") + "-" + args.parameter
    tmpCSV  = os.path.join(Scratch, Tag + "-tmp.csv")
    datFile = os.path.join(Scratch, Tag + "-tmp.dat")
    bmCSV   = os.path.join(Scratch, Tag + "-tmp-beam.csv")
    logFile = os.path.join(Scratch, Tag + ".log")

    print(" scanParameter: start;", args.element, args.parameter, \
          "in", beamlinefile)
    Results = []
    for Value in Values:
        if not writeBeamLine(beamlinefile, tmpCSV, args.element, \
                             args.parameter, args.occurrence, Value):
            print("     ----> No line", args.occurrence, "with Element =", \
                  args.element, "and Parameter =", args.parameter)
            return 1

        simArgs = ["-b", tmpCSV, "-o", datFile, "-n", str(args.nEvts)]
        if inputfile is not None:
            simArgs += ["-i", inputfile]
            if args.bdsimfile:
                simArgs += ["-z", "True"]
        if not runScript("runBeamSim.py", simArgs, logFile) or \
           not runScript("plotBeam.py", \
                         ["-i", datFile, "-n", str(args.nEvts), "-o", bmCSV], \
                         logFile):
            print("     ----> Failed for value", Value, "; see", logFile)
            return 1

        Name, nPrtcls, sigx, sigy = readEndOfBeamLine(bmCSV)
        Results.append([Value, nPrtcls, sigx, sigy])
        print("     ----> Value:", Value, " particles at", Name, ":", \
              nPrtcls, " sigma_x, sigma_y (m):", sigx, sigy)

    outCSV = os.path.join(Scratch, Tag + ".csv")
    with open(outCSV, 'w', newline='') as File:
        Writer = csv.writer(File)
        Writer.writerow(["Value", "nPrtcls_end", "sigma_x_end", "sigma_y_end"])
        Writer.writerows(Results)

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, \
                            layout="constrained")
    axs[0].plot([R[0] for R in Results], \
                [1000.*R[2] if R[2] is not None else float('nan') \
                 for R in Results], 'o-', label=r'$\sigma_x$')
    axs[0].plot([R[0] for R in Results], \
                [1000.*R[3] if R[3] is not None else float('nan') \
                 for R in Results], 's-', label=r'$\sigma_y$')
    axs[0].set_ylabel("RMS beam size at end (mm)")
    axs[0].legend()
    axs[1].plot([R[0] for R in Results], [R[1] for R in Results], 'o-')
    axs[1].set_ylabel("Particles at end")
    axs[1].set_xlabel(args.element + " " + args.parameter)
    fig.savefig(os.path.join(Scratch, Tag + ".pdf"))
    plt.close(fig)

    print(" <---- scanParameter: done; results in", outCSV)
    return 0

"""
   Execute main"
"""
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

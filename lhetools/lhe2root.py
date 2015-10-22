#! /usr/bin/env python

# Jim Henderson January 2013
# James.Henderson@cern.ch
#
# Usage:
# lhe2root.py <input_file.lhe> <OPTIONAL: output_file_name.root>
#
# PLEASE NOTE: This conversion was generated to convert les houches 1.0, it may not work on other versions
#              Please check the les houches version # at the top of the .lhe file

import os, sys
import ROOT as r
from ROOT import TTree, TFile, AddressOf, gROOT
import numpy as np
import gzip

# Get the input lhe file
if len(sys.argv) < 3:
    print "Run as python lhetools.py input.lhe output.root"
    sys.exit(1)

fn = sys.argv[1]
ofn = sys.argv[2]

try:
    input_file = None
    if fn.endswith("gz"):
        input_file = gzip.open( sys.argv[1], 'r')
    elif fn.endswith("lhe"):
        input_file = file( sys.argv[1], 'r')
    else:
        raise Exception("Unknown file format.")
except Exception as e:
    print "\nThe entered file cannot be opened, please enter a valid .lhe (.lhe.gz) file. Exiting. \n"
    print e
    sys.exit(1)

try:
    output_file = TFile(ofn, "RECREATE")
except:
    print "Cannot open output file named: " + ofn + "\nPlease enter a valid output file name as the 2nd arguement. Exiting"
    sys.exit(1)
    pass

output_tree = TTree("events", "events")
print "Setup complete \nOpened file " + str(fn) + "  \nConverting to .root format and outputing to " + ofn

NMAX = 100
# Setup output branches
PIDX_v = np.zeros(NMAX, "i")  
STATUS_v = np.zeros(NMAX, "i")  
PID_v = np.zeros(NMAX, "i")  
MID1_v = np.zeros(NMAX, "i")  
MID2_v = np.zeros(NMAX, "i")  
P_X_v =  np.zeros(NMAX, "d")
P_Y_v = np.zeros(NMAX, "d")
P_Z_v = np.zeros(NMAX, "d")
E_v = np.zeros(NMAX, "d")
M_v = np.zeros(NMAX, "d")
n_particles = np.zeros(1, "i")

# Create a struct which acts as the TBranch for non-vectors
output_tree.Branch('n_particles', n_particles, 'n_particles/I')
output_tree.Branch("PIDX", PID_v, "PIDX[n_particles]/I")
output_tree.Branch("STATUS", STATUS_v, "STATUS[n_particles]/I")
output_tree.Branch("PID", PID_v, "PID[n_particles]/I")
output_tree.Branch("MID1", MID1_v, "MID1[n_particles]/I")
output_tree.Branch("MID2", MID2_v, "MID2[n_particles]/I")
output_tree.Branch("P_X", P_X_v, "P_X[n_particles]/D")
output_tree.Branch("P_Y", P_Y_v, "P_Y[n_particles]/D")
output_tree.Branch("P_Z", P_Z_v, "P_Z[n_particles]/D")
output_tree.Branch("E", E_v, "E[n_particles]/D")
output_tree.Branch("M", M_v, "M[n_particles]/D")

skippedLines = []
in_ev = 0 # To know when to look for particles we must know when we are inside an event
in_ev_1 = 0 # The first line after <event> is information so we must skip that as well

for line in input_file:
    
    if line[:1] == "#":
        continue
    
    if in_ev_1 == 1:
        in_ev_1 = 0
        #s.weight = float( line.split()[2] )
        in_ev = 1
        continue
    
    if line.strip().startswith("<event>"):
        in_ev_1 = 1
        continue
    
    # Some versions of les houches have a pdf line that we don't care about here
    if line.startswith("#pdf"):
        continue
    
    if in_ev == 1 and (line.strip().startswith("</event>") or line.strip().startswith("<rwgt>")):
        output_tree.Fill()
        # Reset variables
        n_particles[0] = 0
        PIDX_v[:] = 0
        STATUS_v[:] = 0
        PID_v[:] = 0
        MID1_v[:] = 0
        MID2_v[:] = 0
        P_X_v[:] = 0.0
        P_Y_v[:] = 0.0
        P_Z_v[:] = 0.0
        E_v[:] = 0.0
        M_v[:] = 0.0
        in_ev = 0
        continue
    
    if in_ev == 1:
        # Check the status of this particle
        try:
            spl = line.split()
            if int(spl[1]) in [-1,1,2]:
                # We have a final state particle on this line
                PIDX_v[n_particles[0]] = n_particles[0]
                STATUS_v[n_particles[0]] = int(spl[1])
                MID1_v[n_particles[0]] = int(spl[2])
                MID2_v[n_particles[0]] = int(spl[3])
                PID_v[n_particles[0]] =  int(spl[0])
                P_X_v[n_particles[0]] =  float(spl[6])
                P_Y_v[n_particles[0]] =   float(spl[7])
                P_Z_v[n_particles[0]] =   float(spl[8])
                E_v[n_particles[0]] =  float(spl[9])
                M_v[n_particles[0]] =  float(spl[10])
                n_particles[0] += 1
                pass
        except:
            if line not in skippedLines:
                print "Problem with line: ", line.replace("\n","")
                print "Skipping..."
                skippedLines.append( line )
                pass

output_tree.Write()
output_file.Close()

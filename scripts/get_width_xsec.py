#!/usr/bin/env python3
"""
Extract width and xsec from gridpack log and madevents
Save the width/xsec into a csv file
Christina Wang 10/22/2025

"""


import numpy as np
import sys
import re
import csv
import os
from utilities import *

def get_decay_value(filepath):
    """
    Extract the decay width value from a line starting with 'DECAY 9900012'
    Example line: DECAY 9900012 1.525792e-15
    """
    pattern = re.compile(r"^DECAY\s+9900012\s+([0-9eE\+\-\.]+)")
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            match = pattern.match(line.strip())
            if match:
                return float(match.group(1))
    return None
def get_xsec(file_path):
    """
    return cross section and uncertainty in units of pb
    """

    
    # Regex pattern to match the section and extract the cross-section
    pattern = r"=== Results Summary for run: pilotrun tag: WJets ===\s+Cross-section\s*:\s*([0-9.eE+-]+)\s*\+-\s*([0-9.eE+-]+)\s*pb"
    
    with open(file_path, "r") as file:
        content = file.read()
    
    match = re.search(pattern, content, re.MULTILINE)
    if match:
        cross_section = match.group(1)
        uncertainty = match.group(2)
        return cross_section, uncertainty
    else:
        return None,None




def ctauFromWidth(width):
    return hbar*c/width

def widthFromCtau(ctau):
    return hbar*c/ctau

hbar = 6.58211915e-25 # hbar in GeV s
c = 299792458000.0 # speed of light in mm/s
# Output CSV file
output_file = "../data/decay_results_coupling0p01.csv"

BASE_PATH=f"{get_git_root()}/genproductions_scripts/bin/MadGraph5_aMCatNLO"

# Write header
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["name", "width [GeV]", "ctau [mm]", "xsec[pb]","uncertainty on xsec[pb]"])

    for decay in ["e", "mu", "tau"]:
        for m in np.arange(1., 11, 1):
            m_str = str(m).replace('.', 'p')
            for coupling in [0.01]: 
            #for ctau in [10,100,1000,10000]:
                #name = f"HNL_{decay}_mN_{m_str}_ctau_{ctau}_13p6TeV"
                name = f"HNL_{decay}_mN_{m_str}_coupling_0p01_13p6TeV"
                filepath = f"{BASE_PATH}/{name}/{name}_gridpack/work/gridpack/process/madevent/Cards/param_card.dat"
                width = get_decay_value(filepath)
                gridpack_log = f"{BASE_PATH}/{name}/{name}_gridpack/work/gridpack/gridpack_generation.log"
                xsec,unc = get_xsec(gridpack_log)



                if width and xsec and unc:
                    ctau = ctauFromWidth(width)
                    writer.writerow([name, width, ctau, xsec, unc])
                else:
                    print(f"{name}: ERROR!!")
                    print(width, xsec,unc)
                    writer.writerow([name, "N/A", "N/A","N/A","N/A","N/A"])

print(f"\n Results saved to {os.path.abspath(output_file)}")



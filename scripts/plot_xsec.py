import os
import matplotlib.pyplot as plt
import numpy as np
from utilities import *

# Load data from CSV
# CSV should have columns: mass, ctau, cross_section

# Get unique masses
os.makedirs(f"{get_git_root()}/plots/", exist_ok=True)

masses = [1,2,3,4,5,6,8,10]
ctau = [10,100,1000, 10000]
for flavor in ['e','mu','tau']:
    plt.figure(figsize=(8,6))
    xsec = {}
    for mass in masses:
        xsec[mass] = []
        for ct in ctau:
            xsec[mass].append(get_xsec(flavor,mass, coupling = None, ctau=ct))
        plt.plot(ctau, xsec[mass],  marker='o', label=f'Mass = {mass} GeV')
    
    plt.xlabel('ctau [mm]')
    plt.ylabel('Cross-section [pb]')
    plt.title('Cross-section vs ctau')
    plt.xscale('log')  # Optional: if ctau spans many orders of magnitude
    plt.yscale('log')  # Optional: if cross-section spans many orders of magnitude
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'{get_git_root()}/plots/cross_section_{flavor}.png')  # Save figure
    

import csv

import subprocess

def get_git_root():
    return subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"],
        universal_newlines=True
    ).strip()



def get_ctau(flavor, mass, coupling, csv_file="/data/decay_results_coupling0p01.csv"):
    """
    Read a CSV file and return ctau for a given HNL flavor, mass, and coupling.

    Args:
        csv_file (str): Path to the CSV file.
        flavor (str): Flavor of HNL ('e', 'mu', or 'tau').
        mass (float): HNL mass (e.g. 5.0).
        coupling (float): Coupling value (e.g. 0.01).

    Returns:
        float or None: The ctau value in mm, or None if not found.
    """
    print(csv_file)
    if type(mass) == int: m_str = str(mass)+'p0'
    else: m_str = str(mass).replace('.', 'p')
    coupling_str = str(coupling).replace('.', 'p')
    default_coupling = "0p01"
    target_name = f"HNL_{flavor}_mN_{m_str}_coupling_{default_coupling}_13p6TeV"
    ctau_default = None
    with open(f"{get_git_root()}/{csv_file}", newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["name"] == target_name: ctau_default = float(row["ctau [mm]"])
    if ctau_default is None: print("mass and flavor combination not found in csv")
    return ctau_default * float(default_coupling.replace('p','.'))**2/coupling**2
def get_coupling(flavor, mass, ctau,  csv_file=f"/scripts/decay_results_coupling0p01.csv"):
    """
    Read a CSV file and return coupling for a given HNL flavor, mass, and ctau

    Args:
        csv_file (str): Path to the CSV file.
        flavor (str): Flavor of HNL ('e', 'mu', or 'tau').
        mass (float): HNL mass (e.g. 5.0).
        ctau (float): ctau in mm

    Returns:
        float or None: The coupling
    """
    if type(mass) == int: m_str = str(mass)+'p0'
    else: m_str = str(mass).replace('.', 'p')
    default_coupling = "0p01"
    target_name = f"HNL_{flavor}_mN_{m_str}_coupling_{default_coupling}_13p6TeV"
    ctau_default = None
    with open(f"{get_git_root()}/{csv_file}", newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["name"] == target_name: ctau_default = float(row["ctau [mm]"])
    if ctau_default is None: print("mass and flavor combination not found in csv")
    return (ctau_default/ctau* float(default_coupling.replace('p','.'))**2)**0.5




def get_xsec(flavor, mass, coupling = None, ctau = None, csv_file="/scripts/decay_results_coupling0p01.csv"):
    """
    Read a CSV file and return xsec for a given HNL flavor, mass, providing coupling or ctau.

    Args:
        csv_file (str): Path to the CSV file.
        flavor (str): Flavor of HNL ('e', 'mu', or 'tau').
        mass (float): HNL mass (e.g. 5.0).
        coupling (float): Coupling value (e.g. 0.01)
        ctau (float): ctau in mm

    Returns:
        float or None: The xsec value in pb, or None if not found.
    """
    if coupling is None and ctau is None: raise ValueError("Need to provide coupling or ctau")
    if coupling and ctau:raise ValueError("Need to provide coupling or ctau")


    if type(mass) == int: m_str = str(mass)+'p0'
    else: m_str = str(mass).replace('.', 'p')
    coupling_str = str(coupling).replace('.', 'p')
    default_coupling = "0p01"
    target_name = f"HNL_{flavor}_mN_{m_str}_coupling_{default_coupling}_13TeV"
    xsec_default = None
    ctau_default = None
    with open(f"{get_git_root()}/{csv_file}", newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["name"] == target_name: 
                xsec_default = float(row["xsec[pb]"])
                ctau_default = float(row["ctau [mm]"])
    if xsec_default is None or ctau_default is None: 
        raise ValueError(f"{flavor} {mass} {target_name} mass and flavor combination not found in csv")


    if coupling: return xsec_default / float(default_coupling.replace('p','.'))**2 * coupling**2
    if ctau: return xsec_default * ctau_default / ctau
    return None

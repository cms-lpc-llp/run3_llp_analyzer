#!/usr/bin/env python3
import os
import json
from pathlib import Path

if __name__ == '__main__':
    '''
    File for setting up static jupyter kernel with necessary python packages. We iterate through the python packages on CMSSW, find the required path for PYTHONPATH, and append to an existing kernel.json
    file used to configure the jupyter kernel. I found this was needed to use jupyter notebooks for VSCode with both PyRoot and local python packages - setting the PYTHONPATH to point to where pyroot is seems
    to unset the path to the other python packages, so this is my quick and dirty solution. Ignore the errors, and the occasional pip install for pkgs not covered may be required.
    '''
    #pkg_strings = ["numpy", "pandas", "matplotlib", "uproot", "coffea", "mplhep", "json", "importlib"]
    path_str = "/cvmfs/cms.cern.ch/el8_amd64_gcc12/external"
    external_path = Path(path_str)
    python_dirs = list(external_path.glob("*py3*"))
    pkg_strings = [str(pkg)[str(pkg).rfind("py3-")+4:] for pkg in python_dirs]
    print(pkg_strings)
    #exit()
    paths=["/cvmfs/cms.cern.ch/el8_amd64_gcc12/external/python3/3.9.14-c10287ae9cadff55334e60003302c349/lib/python3.9"] #start with generic path, has json, etc
    for pkg in pkg_strings:
        pkg_underscore = pkg.replace("-","_")
        output = os.popen(f"python3 -c \"import {pkg_underscore}; print({pkg_underscore}.__file__)\"")
        output_str = output.read()
        output_str = output_str[:output_str.rfind(pkg)-1]
        print(output_str)
        if output_str==None:print("No site packages found");continue
        if output_str not in paths: paths.append(output_str)
        
    print("outputlist ", paths)
    
    large_string = ""
    for single_str in paths: 
        if len(large_string)>0 and large_string[0]!=":":
            large_string=single_str+":"+large_string
        else:
            large_string=single_str+large_string
    #print(large_string)
    original_kernel_path = "/uscms/homes/a/amalbert/.local/share/jupyter/kernels/pyroot_local/kernel.json" #user specific
    with open(original_kernel_path, "r") as f:
        original_dict = json.load(f)
        new_PYTHONPATH = large_string+":"+original_dict["env"]["PYTHONPATH"]
    #print(new_PYTHONPATH)

    new_dict = original_dict.copy()
    new_dict["env"]["PYTHONPATH"] = new_PYTHONPATH
    #print("newdict:", new_dict)
    with open(original_kernel_path, "w") as f:
        json.dump(new_dict, f)
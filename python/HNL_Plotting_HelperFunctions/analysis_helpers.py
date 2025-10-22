#!/usr/bin/env python3

import sys
import os

import coffea
import awkward as ak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot
import numpy as np
import hist
import yaml
import pandas as pd
import dask_awkward as dak
import dask
import dask.dataframe as dd
from dask.distributed import Client
import matplotlib.pyplot as plt
import ROOT


def get_ROOT_1D_fromhistpkg(hist):
    '''
    Function to take hist.Hist and convert to a root hist
    '''
    edges = hist.axes[0].edges
    matplotlib_hist = hist.view()

    root_hist = ROOT.TH1F("num_hist", "numerator hist", np.size(edges)-1, edges)

    for idx in range(len(matplotlib_hist)):
        root_hist.SetBinContent(idx, matplotlib_hist[idx])

    return root_hist

def get_ROOT_hist_from_ak(array, xlow, xhigh, nbins):
    '''
    Function for returning ROOT float histogram from awkward array
    '''
    
    root_hist = ROOT.TH1F("hist", "hist", nbins, xlow, xhigh)
    root_hist.Sumw2()
    for cluster in list(array):
        root_hist.Fill(cluster)
    return root_hist
        

def hist_to_fraction(counts, bins):
    '''
    Function to take matplotlib histogram and re-express so each bin has the fraction of events
    '''
    total = np.sum(counts)
    frac = counts / total
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bin_widths = np.diff(bins)
    errors = np.sqrt(counts)/total

    plt.clf() 
    plot = plt.step(
        bins[:-1],
        frac,
    )
    errors = plt.errorbar(bin_centers, frac, yerr=errors)

    return plot, errors


def get_ROOT_2D_from_matplotlib(hist1, hist2):
    hist1_root = get_ROOT_1D_from_matplotlib(hist1)
    hist2_root = get_ROOT_1D_from_matplotlib(hist2)



def convert_density_to_fraction():
    pass
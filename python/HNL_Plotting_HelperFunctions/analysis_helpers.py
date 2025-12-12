#!/usr/bin/env python3

import sys
import os

import coffea
import awkward as ak
from coffea import processor

#from coffea.nanoevents.methods import vector
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

def auto_legend_position(h):
    max_bin = h.GetMaximumBin()
    x = h.GetXaxis().GetBinCenter(max_bin)
    # If max is on left, put legend on right, etc.
    if x < 0.5*(h.GetXaxis().GetXmax() + h.GetXaxis().GetXmin()):
        return True
    else:
        return False

def get_ROOT_hist_from_ak(array, xlow, xhigh, nbins):
    '''
    Function for returning ROOT float histogram from awkward array
    '''
    
    root_hist = ROOT.TH1F("hist", "hist", nbins, xlow, xhigh)
    root_hist.Sumw2()
    for cluster in list(array):
        root_hist.Fill(cluster)
    return root_hist


def get_ROOT_1D_fromhistpkg(hist, fraction=True):
    '''
    Function to take hist.Hist and convert to a root hist
    '''
    edges = hist.axes[0].edges
    #print(edges)
    hist_counts = hist.view()
    if fraction:
        total = np.sum(hist_counts)
    else:
        total = 1.0
    frac = hist_counts / total
    error = np.sqrt(hist_counts) / total


    root_hist = ROOT.TH1F("num_hist", "numerator hist", np.size(edges)-1, edges)

    for idx in range(len(hist_counts)):
        root_hist.SetBinContent(idx+1, frac[idx])
        root_hist.SetBinError(idx+1, error[idx])

    return root_hist
        

def hist_dict_to_fraction_dict(hist_dict):
    '''
    Function to convert a dictionary of histograms to a dictionary of fraction histograms
    '''
    frac_dict = {}
    
    for key, hist in hist_dict.items():
        #print(hist.metadata['title'])
        frac_dict[key] = (get_ROOT_1D_fromhistpkg(hist), hist.metadata['title'], hist.metadata['x_label'])

    return frac_dict

def plot_one_root_hists(root_hist_1, title, xlabel, label1, logscale=True, fraction=True, outdir='one_hist_plots/'):
    '''
    Function to take 3 root histograms and plot them on same canvas
    '''
    c=ROOT.TCanvas()
    ROOT.gStyle.SetOptStat(0)
    if fraction:
        fraction_label = "Fraction of Events"
    else:
        fraction_label = "Counts"
    if logscale:
        c.SetLogy()
    if auto_legend_position(root_hist_1):
        legend = ROOT.TLegend(0.6,0.7,0.9,0.8)
    else:
        legend = ROOT.TLegend(0.1,0.7,0.4,0.8)
    root_hist_1.SetLineColor(ROOT.kRed)
    root_hist_1.SetLineWidth(2)
    root_hist_1.Draw("HIST E")
    root_hist_1.SetTitle(f"{title};{xlabel};{fraction_label}")
    root_hist_1.GetYaxis().SetRangeUser(1e-4, 10)
    legend.AddEntry(root_hist_1, label1, "lp")
    legend.Draw()
    c.SaveAs(f"{outdir}/{title}.png")


def plot_hists_one_dicts(hist_dict_1, label1, logscale=True, fraction=True, outdir='one_hist_plots/'):
    '''
    Function to take root histograms in 2 dictionaries and plot corresponding histograms on same canvas
    e.g. use for e, and electron-triggered data
    '''
    os.makedirs(outdir, exist_ok=True)
    for key in hist_dict_1.keys():
        root_hist_1, title, x_label = hist_dict_1[key]
        plot_one_root_hists(root_hist_1, title, x_label, label1, logscale=logscale, fraction=fraction, outdir=outdir)
    return

def plot_two_root_hists(root_hist_1, root_hist_2, title, xlabel, label1, label2, logscale=True, fraction=True, outdir='three_hist_plots/'):
    '''
    Function to take 3 root histograms and plot them on same canvas
    '''
    c=ROOT.TCanvas()
    ROOT.gStyle.SetOptStat(0)
    if fraction:
        fraction_label = "Fraction of Events"
    else:
        fraction_label = "Counts"
    if logscale:
        c.SetLogy()
    if auto_legend_position(root_hist_1):
        legend = ROOT.TLegend(0.6,0.7,0.9,0.85)
    else:
        legend = ROOT.TLegend(0.1,0.7,0.4,0.85)
    root_hist_1.SetLineColor(ROOT.kRed)
    root_hist_2.SetLineColor(ROOT.kBlack)
    root_hist_1.SetLineWidth(2)
    root_hist_2.SetLineWidth(2)
    root_hist_1.Draw("HIST E")
    root_hist_2.Draw("HIST E SAME")
    root_hist_1.SetTitle(f"{title};{xlabel};{fraction_label}")
    root_hist_1.GetYaxis().SetRangeUser(1e-4, 10)
    legend.AddEntry(root_hist_1, label1, "lp")
    legend.AddEntry(root_hist_2, label2, "lp")
    legend.Draw()
    c.SaveAs(f"{outdir}/{title}.png")


def plot_hists_2_dicts(hist_dict_1, hist_dict_2, label1, label2, logscale=True, fraction=True, outdir='two_hist_plots/'):
    '''
    Function to take root histograms in 2 dictionaries and plot corresponding histograms on same canvas
    e.g. use for e, and electron-triggered data
    '''
    os.makedirs(outdir, exist_ok=True)
    for key in hist_dict_1.keys():
        root_hist_1, title, x_label = hist_dict_1[key]
        root_hist_2, _, _ = hist_dict_2[key]
        plot_two_root_hists(root_hist_1, root_hist_2, title, x_label, label1, label2, logscale=logscale, fraction=fraction, outdir=outdir)
    return



def plot_three_root_hists(root_hist_1, root_hist_2, root_hist_3, title, xlabel, label1, label2, label3, logscale=True, fraction=True, outdir='three_hist_plots/'):
    '''
    Function to take 3 root histograms and plot them on same canvas
    '''
    c=ROOT.TCanvas()
    ROOT.gStyle.SetOptStat(0)
    if fraction:
        fraction_label = "Fraction of Events"
    else:
        fraction_label = "Counts"
    if logscale:
        c.SetLogy()
    if auto_legend_position(root_hist_1):
        legend = ROOT.TLegend(0.6,0.7,0.9,0.9)
    else:
        legend = ROOT.TLegend(0.1,0.7,0.4,0.9)
    root_hist_1.SetLineColor(ROOT.kRed)
    root_hist_2.SetLineColor(ROOT.kBlack)
    root_hist_3.SetLineColor(ROOT.kViolet)
    root_hist_1.SetLineWidth(2)
    root_hist_2.SetLineWidth(2)
    root_hist_3.SetLineWidth(2)
    root_hist_1.Draw("HIST E")
    root_hist_2.Draw("HIST E SAME")
    root_hist_3.Draw("HIST E SAME")
    root_hist_1.SetTitle(f"{title};{xlabel};{fraction_label}")
    root_hist_1.GetYaxis().SetRangeUser(1e-4, 10)
    legend.AddEntry(root_hist_1, label1, "lp")
    legend.AddEntry(root_hist_2, label2, "lp")
    legend.AddEntry(root_hist_3, label3, "lp")
    legend.Draw()
    c.SaveAs(f"{outdir}/{title}.png")


def plot_hists_3_dicts(hist_dict_1, hist_dict_2, hist_dict_3, label1, label2, label3, logscale=True, fraction=True, outdir='three_hist_plots/'):
    '''
    Function to take root histograms in 3 dictionaries and plot corresponding histograms on same canvas
    e.g. use for e-type HNL, tau decay to electron, and electron-triggered data
    '''
    os.makedirs(outdir, exist_ok=True)
    for key in hist_dict_1.keys():
        root_hist_1, title, x_label = hist_dict_1[key]
        root_hist_2, _, _ = hist_dict_2[key]
        root_hist_3, _, _ = hist_dict_3[key]
        plot_three_root_hists(root_hist_1, root_hist_2, root_hist_3, title, x_label, label1, label2, label3, logscale=logscale, fraction=fraction, outdir=outdir)
    return

def plot_four_root_hists(root_hist_1, root_hist_2, root_hist_3, root_hist_4, title, xlabel, label1, label2, label3, label4, logscale=True, fraction=True, outdir='four_hist_plots/'):
    '''
    Function to take 4 root histograms and plot them on same canvas
    '''
    c=ROOT.TCanvas()
    ROOT.gStyle.SetOptStat(0)
    if fraction:
        fraction_label = "Fraction of Events"
    else:
        fraction_label = "Counts"
    if logscale:
        c.SetLogy()
    if auto_legend_position(root_hist_1):
        legend = ROOT.TLegend(0.6,0.7,0.9,0.95)
    else:
        legend = ROOT.TLegend(0.1,0.7,0.4,0.95)
    root_hist_1.SetLineColor(ROOT.kRed)
    root_hist_2.SetLineColor(ROOT.kBlack)
    root_hist_3.SetLineColor(ROOT.kViolet)
    root_hist_4.SetLineColor(ROOT.kBlue)
    root_hist_1.SetLineWidth(2)
    root_hist_2.SetLineWidth(2)
    root_hist_3.SetLineWidth(2)
    root_hist_4.SetLineWidth(2)
    root_hist_1.Draw("HIST E")
    root_hist_2.Draw("HIST E SAME")
    root_hist_3.Draw("HIST E SAME")
    root_hist_4.Draw("HIST E SAME")
    root_hist_1.SetTitle(f"{title};{xlabel};{fraction_label}")
    root_hist_1.GetYaxis().SetRangeUser(1e-4, 10)
    legend.AddEntry(root_hist_1, label1, "lp")
    legend.AddEntry(root_hist_2, label2, "lp")
    legend.AddEntry(root_hist_3, label3, "lp")
    legend.AddEntry(root_hist_4, label4, "lp")
    legend.Draw()
    c.SaveAs(f"{outdir}/{title}.png")


def plot_hists_4_dicts(hist_dict_1, hist_dict_2, hist_dict_3, hist_dict_4, label1, label2, label3, label4, logscale=True, fraction=True, outdir='four_hist_plots/'):
    '''
    Function to take root histograms in 4 dictionaries and plot corresponding histograms on same canvas
    e.g. use for e-type HNL, tau decay to electron, and electron-triggered data
    '''
    os.makedirs(outdir, exist_ok=True)
    for key in hist_dict_1.keys():
        root_hist_1, title, x_label = hist_dict_1[key]
        root_hist_2, _, _ = hist_dict_2[key]
        root_hist_3, _, _ = hist_dict_3[key]
        root_hist_4, _, _ = hist_dict_4[key]
        plot_four_root_hists(root_hist_1, root_hist_2, root_hist_3, root_hist_4, title, x_label, label1, label2, label3, label4, logscale=logscale, fraction=fraction, outdir=outdir)
    return

def make_datacard_2tag(outDataCardsDir,modelName,  signal_rate, norm, bkg_rate, observation, sig_unc, sig_unc_name,signal_region, prefix, mass, ctau):
    a,b,c,d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    nSig = len(signal_rate.keys())
    os.makedirs(outDataCardsDir+modelName, exist_ok=True)
    text_file = open(outDataCardsDir+modelName+"/mN{0}_ctau{1}_norm{2:.5f}.txt".format(mass, ctau, norm), "w")
    text_file.write('# signal norm {0} \n'.format(norm))

    text_file.write('imax {0} \n'.format(4))
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')



    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \t chD \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n'.format(observation[0],observation[1],observation[2],observation[3]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig) +'\t chD '*(1+nSig) +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * 4 + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * 4 + '\n')
    rate_string = 'rate'
    for i in range(4):# 4 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')

    text_file.write(prefix+'A   rateParam       chA     bkg     (@0*@2/@1)     '+prefix+'B,'+prefix+'C,'+prefix+'D \n')
    #text_file.write(prefix+'A   rateParam       chA     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(a, a*7))
    if b == 0: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, c*7))
    else: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7))
    text_file.write(prefix+'C   rateParam       chC     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(c, c*7))
    if d == 0:text_file.write(prefix+'D   rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(d, c*7))
    else: text_file.write(prefix+'D   rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(d, d*7))


    for k,v in signal_rate.items():
        text_file.write('norm rateParam * {0} {1:.2f}  \n'.format(k,norm))


  #### uncertainties ####
    for k,v in sig_unc.items():assert(len(sig_unc_name)==len(v))
    for i in range(len(sig_unc_name)):
        if 'mc_stats' in sig_unc_name[i]:
            for j, bin in enumerate(['A', 'B', 'C', 'D']):#bin
                    for l, k in enumerate(sig_unc.keys()): #channels
                        before = (len(sig_unc.keys())+1)*j+l
                        after = (len(sig_unc.keys())+1)*4-before-1
                        if sig_unc[k][i][j] > 0.0: text_file.write(sig_unc_name[i]+'_'+k+'_'+bin+' \t gmN ' +str(int(sig_unc[k][i][j]))+ '  '+'\t -  '*before + str(signal_rate[k][j]/int(sig_unc[k][i][j])) + '\t - '*after +'\n')

        else:

            unc_text = sig_unc_name[i]+' \t lnN'
            if len(sig_unc[list(sig_unc.keys())[0]][i])==4:#symmetric uncertainties
                for j in range(4):#bin
                    for k,v in sig_unc.items():
                        if v[i][j] == 0.0:unc_text += ' \t -'
                        else: unc_text += ' \t '+str(v[i][j]+1)
                    unc_text += '\t - '
            else:#asymmetric
                for j in range(4):#bin A, B, C, D
                    for k,v in sig_unc.items():
                        if  v[i][j] == 0.0 and v[i][j+4] == 0.0: unc_text += ' \t -'
                        else:unc_text += ' \t {0}/{1}'.format(1-v[i][j],1+v[i][j+4])
                    unc_text += '\t -'
            text_file.write(unc_text + ' \n')
    # for i in range(len(bkg_unc_name)):
    #     bkg_unc_text = bkg_unc_name[i] + ' \t lnN ' + '\t - '*(4*nSig+3) + '\t ' + str(1+bkg_unc[i]) + ' \n'
    #     text_file.write(bkg_unc_text)

    text_file.close()

def convert_density_to_fraction():
    pass
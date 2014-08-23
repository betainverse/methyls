#!/usr/bin/env python
"""
Read in a set of sparky peaklist files with columns for Data Height and S/N
for methyl resonances, with files for the 'yes' and 'no' condition at each
of several delay lengths. From these values, compute and chart order parameters
for each methyl group, with error bars, using nonlinear fitting and monte carlo
methods.
"""
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from math import pi
from scipy.stats import norm

# User-editable constant
tauC=10.595e-9 #s, the global molecular tumbling time, (if you have not measured your protein'ss tauC, use the following website to calculate the rough value, http://nickanthis.com/tools/tau)

# Sample Name for titles
sample_name = '20140624 HCN apo CzrA'

# Input file information

# Location of files
FileDirectory = '/Users/edmonds/git/matlab/orderparameter/'
testfile = '20140624_apoCzra_DQ_y_40c_ph6_2ms.list'
# Filenames for Yes condition, keyed by delay length (seconds)
Yfiles = {0.0020:'20140624_apoCzra_DQ_y_40c_ph6_2ms.list', 
          0.0050:'20140624_apoCzra_DQ_y_40c_ph6_5ms.list', 
          0.0080:'20140624_apoCzra_DQ_y_40c_ph6_8ms.list', 
          0.0120:'20140624_apoCzra_DQ_y_40c_ph6_12ms.list',
          0.0170:'20140624_apoCzra_DQ_y_40c_ph6_17ms.list',
          0.0220:'20140624_apoCzra_DQ_y_40c_ph6_22ms.list',
          0.0270:'20140624_apoCzra_DQ_y_40c_ph6_27ms.list'}

# Filenames for No condition, keyed by delay length (seconds)
Nfiles = {0.0020:'20140624_apoCzra_DQ_n_40c_ph6_2ms.list', 
          0.0050:'20140624_apoCzra_DQ_n_40c_ph6_5ms.list', 
          0.0080:'20140624_apoCzra_DQ_n_40c_ph6_8ms.list', 
          0.0120:'20140624_apoCzra_DQ_n_40c_ph6_12ms.list',
          0.0170:'20140624_apoCzra_DQ_n_40c_ph6_17ms.list',
          0.0220:'20140624_apoCzra_DQ_n_40c_ph6_22ms.list',
          0.0270:'20140624_apoCzra_DQ_n_40c_ph6_27ms.list'}


def parsepeaklist(filepath):
    # Given the path to a file, return a pandas DataFrame object indexed by
    # the peak assignment, with columns for Data Height and S/N.
    # sep='\s\s+' ensures that 'Data Height' is treated as one cell, while
    # multiple whitespace characters are combined into single delimiters.
    # For some reason, engine='python' is required for regex separators like \s+.
    # skiprows=[1] removes the empty line after the header found in typical
    # sparky peaklist files.
    # discard the w1 and w2 columns. 
    return pd.read_table(filepath, sep='\s\s+',index_col='Assignment',engine='python',skiprows=[1])[['Data Height','S/N']]

def parselists(Yfiles, Nfiles):
    # Given a pair of dictionaries of filenames for each delay length in the
    # Yes condition and the No condition, obtain DataFrames of data from each file
    delays = Yfiles.keys()
    delays.sort()
    Ydataframes={}
    Ndataframes={}
    for d in delays:
        Ydataframes[d] = parsepeaklist(FileDirectory+Yfiles[d])
        Ndataframes[d] = parsepeaklist(FileDirectory+Nfiles[d])
    # Display input summary. Write to a file?
#    for d in delays:
#        print "\nInput for delay = %0.3f, yes condition:"%d
#        print Ydataframes[d]
#        print "\nInput for delay = %0.3f, no condition:"%d
#        print Ndataframes[d]
    return Ydataframes,Ndataframes

def computeratiossigmas(Yframes,Nframes):
    # Recombine DataFrames to obtain new DataFrame objects containing
    # peak height ratios and error bars for each peak at each delay.
    delays = Yframes.keys()
    delays.sort()
    ratios = DataFrame()
    ratios.columns.names = ['Delay']
    #ratios.name = 'ratios'
    sigmas = DataFrame()
    sigmas.columns.names = ['Delay']
    for d in delays:
        ratios[d] = -1*Yframes[d]['Data Height']/Nframes[d]['Data Height']
        sigmas[d] = ratios[d]*np.sqrt(1/Yframes[d]['S/N']**2+1/Nframes[d]['S/N']**2)
    print "ratios:"
    print ratios
    print "sigmas:"
    print sigmas
    return ratios,sigmas

def fitFunc(t, eta, delta):
    return 0.5*eta*np.tanh(t*np.sqrt(eta**2+delta**2))/(np.sqrt(eta**2+delta**2)-delta*np.tanh(t*np.sqrt(eta**2+delta**2)))

def eta2S2axis(eta):
    # Constants
    mu0=1.2566e-6 #T*m/A, ideal vacuum apedimity constant
    gammaH=2.675e8 #s-1*T-1, proton gyromagnetic ratio
    rHH=1.813e-10 #m, the distance between pairs of methyl protons
    h=6.626E-34 #J*s, Planck constant
    S2axis=(10.0/9)*((4*pi/mu0)**2)*4*(rHH**6)*eta/(tauC*(h/(2*pi))**2*gammaH**4)
    return S2axis

def plotcurve(assignment,allratios,allsigmas):
    delays = allratios.columns.values
    ratios = allratios.ix[assignment]
    sigmas = allsigmas.ix[assignment]

    fitParams, fitCovariances = curve_fit(fitFunc, delays, ratios)
    eta = fitParams[0]
    S2axis = eta2S2axis(eta)
    
    print 'S2axis: %.2f'%S2axis
    print ' fit coefficients [eta,delta]:\n', fitParams
    print ' Covariance matrix:\n', fitCovariances
    plt.ylabel('Peak height ratio '+r'$\frac{I_a}{I_b}$')
    plt.xlabel('delay (s)')
    sigS2axis=0
    S2expression = r'$S_{axis}^2 = %.2f\pm%.2f$'%(S2axis,sigS2axis)
    plt.title(assignment)
    plt.errorbar(delays, ratios, fmt = 'bo', yerr = sigmas)
    plt.plot(delays, fitFunc(delays, fitParams[0], fitParams[1]))
    plt.annotate(S2expression,xy=(10,-10),xycoords='axes points',
                 horizontalalignment='left',verticalalignment='top')
    # save plot to a file
    plt.savefig(assignment+'.pdf', bbox_inches=0, dpi=600)

def S2error(assignment,allratios,allsigmas):
    delays = allratios.columns.values
    ratios = allratios.ix[assignment]
    sigmas = allsigmas.ix[assignment]
    S2s = []
    for k in range(5000):
        generatedratios = np.random.normal(ratios,sigmas)
        fitParams, fitCovariances = curve_fit(fitFunc, delays, generatedratios)
        S2s.append(eta2S2axis(fitParams[0]))
    mu,std = norm.fit(S2s)
    #x = np.linspace(0,1,100)
    #pdf_fitted = norm.pdf(x,loc=mu,scale=std)
    #plt.figure()
    #plt.hist(S2s,30,normed=True)
    #plt.plot(x,pdf_fitted, 'k', linewidth=2)
    #plt.show()
    return std    

def plot1curve(assignment,allratios,allsigmas,ax):
    delays = allratios.columns.values
    ratios = allratios.ix[assignment]
    sigmas = allsigmas.ix[assignment]
    fitParams, fitCovariances = curve_fit(fitFunc, delays, ratios)
    eta = fitParams[0]
    S2axis = eta2S2axis(eta)
    sigS2axis=S2error(assignment,allratios,allsigmas)
    S2expression = r'$S_{axis}^2 = %.2f\pm%.2f$'%(S2axis,sigS2axis)
    ax.errorbar(delays, ratios, fmt = 'bo', yerr = sigmas)
    ax.plot(delays, fitFunc(delays, fitParams[0], fitParams[1]))
    plt.setp(ax.get_xticklabels(),rotation='vertical')
    ax.annotate(assignment+'\n'+S2expression,xy=(10,-10),xycoords='axes points',
                horizontalalignment='left',verticalalignment='top')
    #ax.title(assignment)
    
    
def plot3curves(allratios,allsigmas):
    delays=allratios.columns.values
    assignments = allratios[delays[0]].keys()
    #assignments = ['M32CE-HE','M78CE-HE','M90CE-HE','M32CE-HE','M78CE-HE','M90CE-HE','M32CE-HE','M78CE-HE','M90CE-HE','M32CE-HE','M78CE-HE','M90CE-HE']
    assignments = ['M32CE-HE','M78CE-HE','M90CE-HE','M32CE-HE','M78CE-HE','M90CE-HE']
    rows = 2
    cols = 3
    f, axes = plt.subplots(rows,cols,sharex=True,sharey=True)
    f.subplots_adjust(wspace=0.05,hspace=0.05)
    #axes = (ax1,ax2,ax3)
    i=0
    j=0
    for ass in assignments:
#        print ass
        plot1curve(ass,allratios,allsigmas,axes[i,j])
        j=j+1
        if j>=cols:
            j=0
            i=i+1
#    for ax in axes:
#        plot1curve(assignments[i],allratios,allsigmas,ax)
#        i=i+1
    big_ax=f.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_color('none')
    big_ax.spines['bottom'].set_color('none')
    big_ax.spines['left'].set_color('none')
    big_ax.spines['right'].set_color('none')
    plt.ylabel('Peak height ratio '+r'$\frac{I_a}{I_b}$')
    plt.xlabel('delay (s)',labelpad=20)
    #f.set_tight_layout(True)
    #plt.setp([a.get_yticklabels() for a in f.axes[:-1]], visible=False)
    plt.savefig('stuff.pdf')
    plt.show()

    

def main():
    filepath = FileDirectory+testfile
    Ydataframes,Ndataframes = parselists(Yfiles,Nfiles)
    ratios,sigmas=computeratiossigmas(Ydataframes,Ndataframes)
    plot3curves(ratios,sigmas)
    #print parsepeaklist(filepath)
    #S2error('M32CE-HE',ratios,sigmas)

main()

# To do:
#    Plot all peaks with subplot
#    Compute sigma
#    Plot sigma

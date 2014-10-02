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
from math import pi,ceil
from scipy.stats import norm
import re

# User-editable constant
tauC=10.595e-9 #s, the global molecular tumbling time, (if you have not measured your protein'ss tauC, use the following website to calculate the rough value, http://nickanthis.com/tools/tau)

# Sample Name for titles
sample_name = '20140928 ILVAM apo CzrA'
methylsperpage = 9
monte_carlo_iterations = 5 #Use 5000 for real, 50 to test

# Input file information

# Location of files
FileDirectory = '/Users/edmonds/git/methyls/orderparameter/apo_CzrA_peaklists/'
#testfile = '20140624_apoCzra_DQ_y_40c_ph6_2ms.list'
# Filenames for Yes condition, keyed by delay length (seconds)
Yfiles = {0.0030:'20140928_apo_WT_CzrA_ILVAM_DQy_40c_pd6_3ms.list', 
          0.0050:'20140928_apo_WT_CzrA_ILVAM_DQy_40c_pd6_5ms.list', 
          0.0080:'20140928_apo_WT_CzrA_ILVAM_DQy_40c_pd6_8ms.list', 
          0.0120:'20140928_apo_WT_CzrA_ILVAM_DQy_40c_pd6_12ms.list',
          0.0170:'20140928_apo_WT_CzrA_ILVAM_DQy_40c_pd6_17ms.list',
          0.0220:'20140928_apo_WT_CzrA_ILVAM_DQy_40c_pd6_22ms.list',
          0.0270:'20140928_apo_WT_CzrA_ILVAM_DQy_40c_pd6_27ms.list'}

# Filenames for No condition, keyed by delay length (seconds)
Nfiles = {0.0030:'20140928_apo_WT_CzrA_ILVAM_DQn_40c_pd6_3ms.list', 
          0.0050:'20140928_apo_WT_CzrA_ILVAM_DQn_40c_pd6_5ms.list', 
          0.0080:'20140928_apo_WT_CzrA_ILVAM_DQn_40c_pd6_8ms.list', 
          0.0120:'20140928_apo_WT_CzrA_ILVAM_DQn_40c_pd6_12ms.list',
          0.0170:'20140928_apo_WT_CzrA_ILVAM_DQn_40c_pd6_17ms.list',
          0.0220:'20140928_apo_WT_CzrA_ILVAM_DQn_40c_pd6_22ms.list',
          0.0270:'20140928_apo_WT_CzrA_ILVAM_DQn_40c_pd6_27ms.list'}

# Noise level for both yes and no conditions. If you don't have the same noise
# in both spectra, you need to re-record the specra.
# In sparky st dialog, enter 10000 for the number of points and hit [Recompute]
Noise = {0.0030:1420,
         0.0050:4300,
         0.0080:1400,
         0.0120:3460,
         0.0170:2000,
         0.0220:2440,
         0.0270:1610}

def parsepeaklist(filepath):
    # Given the path to a file, return a pandas DataFrame object indexed by
    # the peak assignment, with columns for Data Height and S/N.
    # sep='\s\s+' ensures that 'Data Height' is treated as one cell, while
    # multiple whitespace characters are combined into single delimiters.
    # For some reason, engine='python' is required for regex separators like \s+.
    # skiprows=[1] removes the empty line after the header found in typical
    # sparky peaklist files.
    # discard the w1 and w2 columns. 
    return pd.read_table(filepath, sep='\s\s+',index_col='Assignment',engine='python',skiprows=[1])[['Data Height']]

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

## def writeratiossigmas(ratios,sigmas):
##     assignments = ratios.transpose().keys()
##     pages = int(ceil(float(len(assignments))/methylsperpage))
##     concatenated = pd.concat([ratios,sigmas],keys=['ratio','sigma']).swaplevel(0,1).sortlevel(0).transpose()
##     concatenated.to_excel(sample_name+'.xls')
##     for i in range(pages):
##         first = i*methylsperpage
##         last = (i+1)*methylsperpage
##         outfile = '%s_%d.csv'%(sample_name,i)
##         concatenated[assignments[first:min(last,len(assignments))]].to_csv(outfile)

def formatAssignment(ass):
    return re.sub(r'(\d+)C',r'\1-C',ass.split('-')[0])
def greekFormatAssignment(ass):
    return r'$%s$'%ass.replace('-CB',' \\beta').replace('-CE',' \epsilon').replace('-CG',' \gamma').replace('-CD',' \delta')


def computeratiossigmas(Yframes,Nframes,Noise):
    # Recombine DataFrames to obtain new DataFrame objects containing
    # peak height ratios and error bars for each peak at each delay.
    delays = Yframes.keys()
    delays.sort()
    ratios = DataFrame()
    ratios.columns.names = ['Delay']
    sigmas = DataFrame()
    sigmas.columns.names = ['Delay']
    for d in delays:
        ratios[d] = -1*Yframes[d]['Data Height']/Nframes[d]['Data Height']
        sigmas[d] = abs(ratios[d])*np.sqrt((Noise[d]/Yframes[d]['Data Height'])**2+(Noise[d]/Nframes[d]['Data Height'])**2)
        # if there are any negative numbers in ratios, set them to zero:
        ratios.loc[ratios[d]<0,d] = 0
    methyls = ratios.index
    ratios.index = [formatAssignment(x) for x in methyls]
    sigmas.index = [formatAssignment(x) for x in methyls]
    #sigmas.index = [x.split('-')[0] for x in methyls]
#    writeratiossigmas(ratios,sigmas)
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

def S2error(assignment,allratios,allsigmas):
    delays = allratios.columns.values
    ratios = allratios.ix[assignment]
    sigmas = allsigmas.ix[assignment]
    S2s = []
    for k in range(monte_carlo_iterations):
        generatedratios = np.random.normal(ratios,sigmas)
        fitParams, fitCovariances = curve_fit(fitFunc, delays, generatedratios)
        S2s.append(eta2S2axis(fitParams[0]))
    mu,std = norm.fit(S2s)
    return std    


def S2barplot(S2errorDF):
    S2values = S2errorDF['S2'].values
    assignments = S2errorDF.index
    #assignments = [greekFormatAssignment(x) for x in S2errorDF.index]
    S2errors = S2errorDF['S2error'].values
    fix,ax = plt.subplots(figsize=(20,5))
    h = plt.bar(xrange(len(assignments)),
                  S2values,
                  color='r',
                  label=assignments,
                  yerr=S2errors)
    plt.subplots_adjust(bottom=0.3)
    xticks_pos = [0.5*patch.get_width() + patch.get_xy()[0] for patch in h]
    plt.xticks(xticks_pos, assignments, ha='right', rotation=45)
    ax.set_ylabel('S2axis')
    #ax.set_xticks(ind+width)
    #ax.set_xticklabels(assignments)
    ax.set_title(sample_name)
    plt.savefig(sample_name+'_bar.pdf')
    plt.show()

        
def plotfakecurve(assignment,allratios,allsigmas,ax):
    delays = allratios.columns.values
    ratios = allratios.ix[assignment]
    sigmas = allsigmas.ix[assignment]
    ax.errorbar(delays, ratios, fmt = 'w.', yerr = sigmas)
    plt.setp(ax.get_xticklabels(),rotation='vertical')

def plot1curve(assignment,allratios,allsigmas,ax):
    delays = allratios.columns.values
    ratios = allratios.ix[assignment]
    sigmas = allsigmas.ix[assignment]
    fitParams, fitCovariances = curve_fit(fitFunc, delays, ratios)
    eta = fitParams[0]
    S2axis = eta2S2axis(eta)
    sigS2axis=S2error(assignment,allratios,allsigmas)
    S2expression = r'$S_{axis}^2 = %.2f\pm%.2f$'%(S2axis,sigS2axis)
    ax.errorbar(delays, ratios, fmt = 'b.', yerr = sigmas)
    ax.plot(delays, fitFunc(delays, fitParams[0], fitParams[1]))
    plt.setp(ax.get_xticklabels(),rotation='vertical')
    ax.annotate(greekFormatAssignment(assignment)+'\n'+S2expression,xy=(10,-10),xycoords='axes points',
                horizontalalignment='left',verticalalignment='top')
    return S2axis, sigS2axis

def plot3curves(allratios,allsigmas):
    delays=allratios.columns.values
    assignments = allratios[delays[0]].keys()    
    print assignments
    rows = 5
    cols = 4
    pages = ceil(float(len(assignments))/(rows*cols))
    f, axes = plt.subplots(rows,cols,sharex=True,sharey=True)
    f.set_size_inches(8,10.5)
    f.subplots_adjust(wspace=0.05,hspace=0.05)
    row=0
    col=0
    page=0
    S2values = []
    S2errors = []
    print "Computing curves:"
    for ass in assignments:
        print ass
        S2,S2err = plot1curve(ass,allratios,allsigmas,axes[row,col])
        S2values.append(S2)
        S2errors.append(S2err)
        col=col+1
        if col>=cols:
            col=0
            row=row+1
        if row >=rows:
            big_ax=f.add_subplot(111)
            big_ax.set_axis_bgcolor('none')
            big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
            big_ax.spines['top'].set_color('none')
            big_ax.spines['bottom'].set_color('none')
            big_ax.spines['left'].set_color('none')
            big_ax.spines['right'].set_color('none')
            big_ax.set_title(sample_name)
            plt.ylabel('Peak height ratio '+r'$\frac{I_a}{I_b}$')
            plt.xlabel('delay (s)',labelpad=20)
            plt.savefig('%s_curves_%d.pdf'%(sample_name,page))
            #plt.show()            
            row=0
            page=page+1
            f, axes = plt.subplots(rows,cols,sharex=True,sharey=True)
            f.set_size_inches(8,10.5)
            f.subplots_adjust(wspace=0.05,hspace=0.05)
    maxcharts = rows*cols*pages
    fakecharts = int(maxcharts-len(assignments))
    for i in range(fakecharts):
        plotfakecurve(ass,allratios,allsigmas,axes[row,col])
        col=col+1
        if col>=cols:
            col=0
            row=row+1
    big_ax=f.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_color('none')
    big_ax.spines['bottom'].set_color('none')
    big_ax.spines['left'].set_color('none')
    big_ax.spines['right'].set_color('none')
    big_ax.set_title(sample_name)
    plt.ylabel('Peak height ratio '+r'$\frac{I_a}{I_b}$')
    plt.xlabel('delay (s)',labelpad=20)
    plt.savefig('%s_curves_%d.pdf'%(sample_name,page))

    S2errorDF = DataFrame({'S2':S2values,'S2error':S2errors}, index=assignments)
    S2errorDF.to_excel(sample_name+'_S2.xls')
    return S2errorDF
    #plt.show()            
    
    

def main():
    Ydataframes,Ndataframes = parselists(Yfiles,Nfiles)
    ratios,sigmas=computeratiossigmas(Ydataframes,Ndataframes,Noise)
    S2errorDF = plot3curves(ratios,sigmas)
    S2barplot(S2errorDF)

main()

# To do:
#    Plot all peaks with subplot
#    Compute sigma
#    Plot sigma
#    Save files with data

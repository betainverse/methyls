#!/usr/bin/env python
"""
Read in a table 
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
import matplotlib.lines as mpllines
import matplotlib.ticker as mplticker

# User-editable constant
# Measure noise level in the st dialog for each spectrum, and take an average.
noise = 4000

# Sample Name for titles
sample_name = '20140905 Zn AdcR T2'

# Input file information

# Location of files
# Use 'rh' in Sparky to generate a peak height table. Use Setup to select
# which spectra to include and to indicate the delay lengths associated
# with each one. Choose either to use 'heights at the same position in each
# spectrum' or 'heights at assigned peak positions only', hit clear if
# necessary, and then invoke 'rh' again once all parameters are set.
# In particular, you must decide whether to indicate the delay lengths in
# seconds or milliseconds. So far, this script assumes delays are indicated
# in milliseconds, and converts them to seconds. 
rhTable = '/Users/edmonds/git/methyls/orderparameter/T2_peaklists/20140905_T2_DCNZnAdcR_heights.txt'
FileDirectory = '/Users/edmonds/git/methyls/orderparameter/T2_peaklists/'
testfile = '20140905_T2_DCNZnAdcR_10.list'
# Filenames for peaklists, keyed by delay length (seconds)
peaklists = {0.010:'20140905_T2_DCNZnAdcR_10.list', 
             0.030:'20140905_T2_DCNZnAdcR_30.list', 
             0.050:'20140905_T2_DCNZnAdcR_50.list', 
             0.090:'20140905_T2_DCNZnAdcR_90.list',
             0.110:'20140905_T2_DCNZnAdcR_110.list',
             0.150:'20140905_T2_DCNZnAdcR_150.list',
             0.190:'20140905_T2_DCNZnAdcR_190.list',
             0.210:'20140905_T2_DCNZnAdcR_210.list'}

# Given a column header containing the name of a spectrum and the
# associated delay, as in '20140905_T2_DCNZnAdcR_10/10', return the delay
# in a string and the original units, in this case, '10'
def header2delay(headerstring):
    return headerstring.split('/')[1]

def parse_rh_table(filepath):
    # Given the path to a file, return a pandas DataFrame object indexed by
    # the peak assignment, ignoring columns with fitted T-decay and SD,
    # collecting the columns that have data heights for each spectrum.
    # The rh tables have column headers separated by a single space, but
    # the data are often separated by more spaces, so sep='\s+' is the only
    # possibility.
    # Convert ms delays to s.
    DF = pd.read_table(filepath, sep='\s+',index_col='Assignment',engine='python',skiprows=[1])
    del DF['T-decay']
    del DF['SD']
    DF.rename(columns=lambda x: 0.001*float(x.split('/')[1]), inplace=True)
    DF.columns.name='delay (s)'
    return DF

def fitFunc(t,R2,I0):
    return I0*np.exp(-1*R2*t)

def interpolate(delays):
    return np.array([float(x)/1000 for x in range(int(max(delays)*1000))])

def plotcurve(assignment,rhDF):
    # Assume delays are already floats expressed in seconds
    delays = rhDF.columns.values
    heights = rhDF.ix[assignment]
    fitParams, fitCovariances = curve_fit(fitFunc,delays,heights)
    T2 = 1/fitParams[0]
    print 'T2: %.2f s'%T2
    print 'I0: %.2f'%fitParams[1]
    print 'Covariance matrix:\n', fitCovariances
    plt.ylabel('Peak height')
    plt.xlabel('delay (s)')
    plt.title(assignment)
    plt.errorbar(delays,heights,fmt='bo',yerr = noise*20)
    xvalues = interpolate(delays)
    plt.plot(xvalues, fitFunc(xvalues, fitParams[0], fitParams[1]))
    T2err = T2error(assignment,rhDF)
    T2expression=r'$T2 = %.1f\pm%.2f$ ms'%(1000*T2,T2err*1000)
    plt.annotate(T2expression,xy=(10,-10),xycoords='axes points',
                 horizontalalignment='left',verticalalignment='top')
    plt.show()
    #plt.savefig(assignment+'.pdf', bbox_inches=0, dpi=600)

def T2error(assignment,rhDF):
    delays = rhDF.columns.values
    heights = rhDF.ix[assignment]
    T2s = []
    for k in range(500):
        generatedheights = np.random.normal(heights,noise)
        fitParams, fitCovariances = curve_fit(fitFunc, delays, generatedheights)
        T2s.append(1/fitParams[0])
    mu,std = norm.fit(T2s)
    #x = np.linspace(0,1,100)
    #pdf_fitted = norm.pdf(x,loc=mu,scale=std)
    #plt.figure()
    #plt.hist(S2s,30,normed=True)
    #plt.plot(x,pdf_fitted, 'k', linewidth=2)
    #plt.show()
    print std
    return std    

def T2barplot(rhDF):
    delays=rhDF.columns.values
    assignments = rhDF[delays[0]].keys()
    intassignments = [int(x[1:-3]) for x in rhDF[delays[0]].keys()]
    T2values = []
    T2errors = []
    for ass in assignments:
        print ass
        heights = rhDF.ix[ass]
        T2s = []
        for k in range(50):
            generatedheights = np.random.normal(heights,noise)
            fitParams, fitCovariances = curve_fit(fitFunc, delays, generatedheights)
            T2s.append(1/fitParams[0])
        mu,std = norm.fit(T2s)
        T2values.append(mu)
        T2errors.append(std)
    fix,ax = plt.subplots()
    h = plt.bar(intassignments,
                  T2values,
                  color='r',
                  label=ass,
                  yerr=T2errors)
    #plt.subplots_adjust(bottom=0.3)
    #xticks_pos = [0.5*patch.get_width() + patch.get_xy()[0] for patch in h]
    #plt.xticks(xticks_pos, intassignments, ha='right', rotation=45)
    ax.set_ylabel('T2 (s)')
    ax.set_xlabel('Residue')
    plt.axis(xmax=ceil((max(intassignments)+1)/10.0)*10)
    ax.xaxis.set_major_locator(mplticker.MultipleLocator(10))
    for line in ax.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)
    plt.gcf().subplots_adjust(bottom=0.15)

    #ax.set_xticks(ind+width)
    #ax.set_xticklabels(assignments)
    ax.set_title(sample_name)
    plt.savefig(sample_name+'_bar.pdf')
    plt.show()
    return T2values,T2errors
        

def plot1curve(assignment,rhDF,T2errors,ax):
    delays = rhDF.columns.values
    heights = rhDF.ix[assignment]
    sigmas = allsigmas.ix[assignment]
    fitParams, fitCovariances = curve_fit(fitFunc, delays, ratios)
    eta = fitParams[0]
    S2axis = eta2S2axis(eta)
    sigS2axis=S2error(assignment,rhDF,allsigmas)
    S2expression = r'$S_{axis}^2 = %.2f\pm%.2f$'%(S2axis,sigS2axis)
    ax.errorbar(delays, ratios, fmt = 'bo', yerr = sigmas)
    ax.plot(delays, fitFunc(delays, fitParams[0], fitParams[1]))
    plt.setp(ax.get_xticklabels(),rotation='vertical')
    ax.annotate(assignment+'\n'+S2expression,xy=(10,-10),xycoords='axes points',
                horizontalalignment='left',verticalalignment='top')
    #ax.title(assignment)
    
    
def plot3curves(rhDF,T2errors):
    delays=rhDF.columns.values
    assignments = rhDF[delays[0]].keys()
    rows = 30
    cols = 5
    f, axes = plt.subplots(rows,cols,sharex=True,sharey=True)
    f.subplots_adjust(wspace=0.05,hspace=0.05)
    #axes = (ax1,ax2,ax3)
    i=0
    j=0
    for ass in assignments:
#        print ass
        plot1curve(ass,rhDF,errors,axes[i,j])
        j=j+1
        if j>=cols:
            j=0
            i=i+1
#    for ax in axes:
#        plot1curve(assignments[i],rhDF,allsigmas,ax)
#        i=i+1
    big_ax=f.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_color('none')
    big_ax.spines['bottom'].set_color('none')
    big_ax.spines['left'].set_color('none')
    big_ax.spines['right'].set_color('none')
    big_ax.set_title(sample_name)
    plt.ylabel('Peak height')
    plt.xlabel('delay (s)',labelpad=20)
    #f.set_tight_layout(True)
    #plt.setp([a.get_yticklabels() for a in f.axes[:-1]], visible=False)
    plt.savefig(sample_name+'_curves.pdf')
    plt.show()

    

def main():
    #filepath = FileDirectory+testfile
    #Ydataframes,Ndataframes = parselists(Yfiles,Nfiles)
    #ratios,sigmas=computeratiossigmas(Ydataframes,Ndataframes)
    #plot3curves(ratios,sigmas)
    #S2barplot(ratios,sigmas)
    #print parsepeaklist(filepath)
    #S2error('M32CE-HE',ratios,sigmas)
    rhDF = parse_rh_table(rhTable)
    #plotcurve('F137N-H',rhDF)
    T2barplot(rhDF)
main()

# To do:
#    Plot all peaks with subplot
#    Compute sigma
#    Plot sigma
#    Save files with data

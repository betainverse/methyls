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
from math import pi,ceil,sqrt
from scipy.stats import norm
import matplotlib.lines as mpllines
import matplotlib.ticker as mplticker

# User-editable constant
# Measure noise level in the st dialog for each spectrum, and take an average.
noise = 4000
yscalefactor = 100000

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

# Filenames for peaklists, keyed by delay length (seconds)
#peaklists = {0.010:'20140905_T2_DCNZnAdcR_10.list', 
#             0.030:'20140905_T2_DCNZnAdcR_30.list', 
#             0.050:'20140905_T2_DCNZnAdcR_50.list', 
#             0.090:'20140905_T2_DCNZnAdcR_90.list',
#             0.110:'20140905_T2_DCNZnAdcR_110.list',
#             0.150:'20140905_T2_DCNZnAdcR_150.list',
#             0.190:'20140905_T2_DCNZnAdcR_190.list',
#             0.210:'20140905_T2_DCNZnAdcR_210.list'}

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
    DF = pd.read_table(filepath, sep='\s+',index_col='Assignment',engine='python')
    del DF['T-decay']
    del DF['SD']
    DF.rename(columns=lambda x: 0.001*float(x.split('/')[1]), inplace=True)
    DF.columns.name='delay (s)'
    return DF

def fitFunc(t,R2,I0):
    return I0*np.exp(-1*R2*t)

#def fitFunc(t,T2,I0):
#    return I0*np.exp(t/T2)

def interpolate(delays):
    return np.array([float(x)/1000 for x in range(int(max(delays)*1000)) if float(x)/1000 >= min(delays)])



## def T2error(assignment,rhDF):
##     delays = rhDF.columns.values
##     heights = rhDF.ix[assignment]
##     T2s = []
##     for k in range(50):
##         generatedheights = np.random.normal(heights,noise)
##         fitParams, fitCovariances = curve_fit(fitFunc, delays, generatedheights)
##         T2s.append(1/fitParams[0])
##     mu,std = norm.fit(T2s)
##     return mu,std    

def T2error(assignment,rhDF):
    delays = rhDF.columns.values
    heights = rhDF.ix[assignment]
    fitParams, fitCovariances = curve_fit(fitFunc, delays, heights)
    R2 = fitParams[0]
    I0 = fitParams[1]
    T2 = 1/R2
    dR = sqrt(fitCovariances[0,0])
    dT = dR*T2/R2
    return T2,dT,I0


def calcT2errors(rhDF):
    delays=rhDF.columns.values
    assignments = rhDF[delays[0]].keys()
    #intassignments = [int(x[1:-3]) for x in rhDF[delays[0]].keys()]
    T2anderrs = {}
    for ass in assignments:
        print ass
        T2anderrs[ass] = T2error(ass,rhDF)
    return T2anderrs

def T2barplot(rhDF):
    delays=rhDF.columns.values
    assignments = rhDF[delays[0]].keys()
    intassignments = [int(x[1:-3]) for x in rhDF[delays[0]].keys()]
    T2anderrs = calcT2errors(rhDF)
    T2values = [T2anderrs[ass][0] for ass in assignments]
    T2errors = [T2anderrs[ass][1] for ass in assignments]
    fix,ax = plt.subplots()
    h = plt.bar(intassignments,
                  T2values,
                  color='r',
                  yerr=T2errors)
    ax.set_ylabel('T2 (s)')
    ax.set_xlabel('Residue')
    plt.axis(xmax=ceil((max(intassignments)+1)/10.0)*10)
    ax.xaxis.set_major_locator(mplticker.MultipleLocator(10))
    for line in ax.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)
    plt.gcf().subplots_adjust(bottom=0.15)
    ax.set_title(sample_name)
    plt.savefig(sample_name+'_bar.pdf')
    plt.show()
    return T2anderrs

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

def plot1curve(assignment,rhDF,ax):
    delays = rhDF.columns.values
    heights = rhDF.ix[assignment]
    yvalues = heights/yscalefactor
    T2,dT,I0 = T2error(assignment,rhDF)
    T2ms = T2*1000
    dTms = dT*1000
    T2expression = r'$%.2f\pm%.2f  ms$'%(T2ms,dTms)
    Assexpression = '%s: T2 ='%assignment[0:-3]
    xvalues = interpolate(delays)
    ax.errorbar(delays, yvalues, fmt = 'b.', yerr = noise/yscalefactor)
    ax.plot(xvalues, fitFunc(xvalues, 1/T2, I0)/yscalefactor)
    plt.setp(ax.get_xticklabels(),rotation='vertical')
    ax.annotate(Assexpression+'\n'+T2expression,xy=(15,-10),xycoords='axes points',
                horizontalalignment='left',verticalalignment='top')

    #ax.title(assignment)
    
    
def plot3curves(rhDF):
    delays=rhDF.columns.values
    assignments = rhDF[delays[0]].keys()
    rows = 3
    cols = 3
    pages = ceil(len(assignments)/(rows*cols))
    f, axes = plt.subplots(rows,cols,sharex=True,sharey=True)
    f.subplots_adjust(wspace=0.05,hspace=0.05)
    #axes = (ax1,ax2,ax3)
    i=0
    j=0
    row=0
    col=0
    page=0
    for ass in assignments:
        plot1curve(ass,rhDF,axes[row,col])
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
            plt.ylabel('Peak height / %d'%yscalefactor)
            plt.xlabel('delay (s)',labelpad=20)
            #f.set_tight_layout(True)
            #plt.setp([a.get_yticklabels() for a in f.axes[:-1]], visible=False)
            plt.savefig('%s_curves_%d.pdf'%(sample_name,page))
            #plt.show()            
            row=0
            page=page+1
            f, axes = plt.subplots(rows,cols,sharex=True,sharey=True)
            f.subplots_adjust(wspace=0.05,hspace=0.05)
        if page>=pages:
            break
        
##     for ass in assignments:
## #        print ass
##         plot1curve(ass,rhDF,axes[i,j])
##         j=j+1
##         if j>=cols:
##             j=0
##             i=i+1
##         if i >= rows:
##             break
## #    for ax in axes:
## #        plot1curve(assignments[i],rhDF,allsigmas,ax)
## #        i=i+1
##     big_ax=f.add_subplot(111)
##     big_ax.set_axis_bgcolor('none')
##     big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
##     big_ax.spines['top'].set_color('none')
##     big_ax.spines['bottom'].set_color('none')
##     big_ax.spines['left'].set_color('none')
##     big_ax.spines['right'].set_color('none')
##     big_ax.set_title(sample_name)
##     plt.ylabel('Peak height / %d'%yscalefactor)
##     plt.xlabel('delay (s)',labelpad=20)
##     #f.set_tight_layout(True)
##     #plt.setp([a.get_yticklabels() for a in f.axes[:-1]], visible=False)
##     plt.savefig(sample_name+'_curves.pdf')
##     plt.show()

    

def main():
    rhDF = parse_rh_table(rhTable)
    #plotcurve('F137N-H',rhDF)
    plot3curves(rhDF)
    T2barplot(rhDF)
main()

# To do:
#    Plot all peaks with subplot
#    Compute sigma
#    Plot sigma
#    Save files with data

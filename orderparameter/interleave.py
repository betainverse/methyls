#!/usr/bin/env python
"""
interleave.py reads a procpar file to extract parameters for invoking the
megaleave program, then arranges all the output files into several experiment 
directories. This program must be invoked from the same directory as the 
procpar and fid files.

Usage:
interleave.py output_prefix
"""
import sys
from subprocess import call
from os import system

path_to_megaleave = "/home/edmonds/bin/interleave_macro/megaleave"

def seconds2ms(secstr):
    return str(int(1000*float(secstr)))

def get_params():
    params = {}
    openfile = open('procpar','r')
    lines = openfile.readlines()
    for i in range(len(lines)):
        if lines[i][0:3] == 'ni ':
            params['ni'] = lines[i+1].split()[1]
            print 'ni = ', params['ni']
        elif lines[i][0:3] == 'np ':
            params['np'] = lines[i+1].split()[1]
            print 'np = ', params['np']
        elif lines[i][0:6] == 'array ':
            if lines[i+1].split()[1] == '\"phase,relaxT\"':
                params['flg'] = '0'
            elif lines[i+1].split()[1] == '\"relaxT,phase\"':
                params['flg'] = '1'
            elif lines[i+1].split()[1] == '\"phase,ncyc_cp\"': 
                params['flg'] = '0'
            elif lines[i+1].split()[1] == '\"ncyc_cp,phase\"': 
               params['flg'] = '1' 
            else:
                print 'array not found'
        elif lines[i][0:7] == 'relaxT ':
            delays_in_seconds = lines[i+1].split()[1:] 
            params['delays'] = [seconds2ms(s) for s in delays_in_seconds]
            params['num'] = lines[i+1].split()[0]
        elif lines[i][0:8] == 'ncyc_cp ': 
            params['delays'] = lines[i+1].split()[1:] # not really "delays"
            params['num'] = lines[i+1].split()[0] 
    return params
                
def run_megaleave(params):
    call([params['path'], 'fid',params['output'],params['num'],params['ni'],params['np'],params['flg']])

def rearrange_files(params):
    for i in range(int(params['num'])):
        orig_filename = '%s_%s'%(params['output'],i)
        new_filename = '%s_%s_%d.fid'%(params['delays'][i],params['output'],i)
        system("/bin/mkdir %s"%new_filename)
        system("/bin/mv %s %s/fid"%(orig_filename,new_filename))
        system("/bin/cp log %s/log"%new_filename)
        system("/bin/cp text %s/text"%new_filename)
        system("/bin/cp procpar %s/procpar"%new_filename)
        

def main():
    if len(sys.argv)!= 2:
        print "\nUsage:"
        print "interleave.py output_prefix "
        #print "\nWARNING: This script assumes you don't have any repeat delay values!!!\n"
        return
    outputname = sys.argv[1]    
    #numexpts = sys.argv[2]
    params = get_params()
    params['output'] = outputname
    params['path'] = path_to_megaleave
    run_megaleave(params)
    rearrange_files(params)
          

main()



from survey_tools import *
import numpy as np, re
from sys import argv

"""
Author: Trevor Dorn-Wallenstein
Date: 12/21/17
Purpose: Call python script_arctic.py region_file filt1 filt2 exp1 exp2 root1 root2 n_per_filt 
to generate the ARCTIC+APO low-level scripting commands that will tile the locations given in 
region_file, using the exposure times given in exp1 and exp2 to expose in filt1 and filt2 respectively.
File names are root1_nfield root2_nfield. n_per_filt is the number of pointings to observe before 
doing a filter change.
If called with the final parameter n_start, will start with the n_start'th field
"""

if __name__ == '__main__':
    
    region_file = argv[1]
    
    filt1 = int(argv[2])
    filt2 = int(argv[3])
    
    exp1 = float(argv[4])
    exp2 = float(argv[5])
    
    root1 = str(argv[6])
    root2 = str(argv[7])
    
    n_per_filt = int(argv[8])
    
    if len(argv) > 9:
        
        n_start = int(argv[9])
    
    else:
        
        n_start = 1
    
    #pull out the strings for the region files of the field
    with open(region_file,'r') as f:
        reg_string = f.read().split('\n')[3:]
        
    ras = []
    decs = []
    
    #Do some parsing to split out the string into ra, dec
    for pos in reg_string:
        
        split = re.split(',|\(',pos)
        
        ras.append(split[1])
        decs.append(split[2])
        
    out_str = ''
    
    n_track = n_start
    
    while n_track < len(ras)+1-n_per_filt:
        
        out_str += 'arctic set filter={}\n'.format(filt1)
    
        for i in range(n_track,n_track+n_per_filt):
            
            ra = ras[i-1]
            dec = decs[i-1]
            
            out_str += 'tcc track {0}, {1}\n'.format(ra,dec)
            out_str += 'arcticExpose object time={0} name={1}_{2}\n'.format(exp1,root1,i)
            
        out_str += 'arctic set filter={}\n'.format(filt2)
        
        for i in range(n_track,n_track+n_per_filt):
            
            ra = ras[i-1]
            dec = decs[i-1]
            
            out_str += 'tcc track {0}, {1}\n'.format(ra,dec)
            out_str += 'arcticExpose object time={0} name={1}_{2}\n'.format(exp2,root2,i)
            
        n_track += n_per_filt
        
    out_str += 'arctic set filter={}\n'.format(filt1)
    
    for i in range(n_track,len(ras)+1):
            
        ra = ras[i-1]
        dec = decs[i-1]

        out_str += 'tcc track {0}, {1}\n'.format(ra,dec)
        out_str += 'arcticExpose object time={0} name={1}_{2}\n'.format(exp1,root1,i)
        
    out_str += 'arctic set filter={}\n'.format(filt2)
    
    for i in range(n_track,len(ras)+1):
            
        ra = ras[i-1]
        dec = decs[i-1]

        out_str += 'tcc track {0}, {1}\n'.format(ra,dec)
        out_str += 'arcticExpose object time={0} name={1}_{2}\n'.format(exp2,root2,i)

    with open('arctic_script.txt','w') as f:
        f.write(out_str)

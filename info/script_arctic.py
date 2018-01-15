from survey_tools import *
import numpy as np, re
from sys import argv

"""
Author: Trevor Dorn-Wallenstein
Date: 12/21/17
Purpose: Call python script_arctic.py region_file filt1 filt2 exp1 exp2 root1 root2 to generate the
ARCTIC+APO low-level scripting commands that will tile the locations given in region_file,
using the exposure times given in exp1 and exp2 to expose in filt1 and filt2 respectively.
If called with the final parameter n_start
"""

if __name__ == '__main__':
    
    region_file = argv[1]
    
    filt1 = int(argv[2])
    filt2 = int(argv[3])
    
    exp1 = float(argv[4])
    exp2 = float(argv[5])
    
    root1 = str(argv[6])
    root2 = str(argv[7])
    
    if len(argv) > 8:
        
        n_start = int(argv[8])
    
    else:
        
        n_start = 1
    
    with open(region_file,'r') as f:
        reg_string = f.read().split('\n')[3:]
        
    ras = []
    decs = []
    
    for pos in reg_string:
        
        split = re.split(',|\(',pos)
        
        ras.append(split[1])
        decs.append(split[2])
        
    out_str = ''
    
    ns = np.arange(1,len(ras)+1)
    
    for ra,dec,n in zip(ras,decs,ns):
        
        if n >= n_start:
            
            out_str += 'MOVE {0} {1}\n'.format(ra,dec)
            out_str += 'arctic set filter={}\n'.format(filt1)
            out_str += 'arcticExpose object time={0} name={1}_{2}\n'.format(exp1,root1,n)
            out_str += 'arctic set filter={}\n'.format(filt2)
            out_str += 'arcticExpose object time={0} name={1}_{2}\n'.format(exp2,root2,n)
        
    with open('arctic_script.txt','w') as f:
        f.write(out_str)
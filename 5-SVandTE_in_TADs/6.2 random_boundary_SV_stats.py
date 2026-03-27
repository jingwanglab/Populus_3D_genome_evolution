# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:17:36 2024

@author: ztt
"""

import os
import pandas as pd
import glob

# List of species
spes = ['Ppse', 'Pwua', 'Psze', 'Plas', 'Pyun', 'Prot', 'Pqio', 'Pdav', 'Pade', 'Psim']

# Define SV types
sv_types = {
    'dup': 'dup_stats.txt',
    'inv': 'inv_stats.txt',
    'indel': 'indel_stats.txt',
    'trans': 'trans_stats.txt'
}

# Base path
base_path = '/usr_storage2/stt/tad/SV/'

for spe in spes:  
    for sv_type, output_file in sv_types.items():
        path = os.path.join(base_path, spe, 'bound_random_sv', sv_type)                
        # Get all overlap files
        files = [f for f in os.listdir(path) if f.endswith('.overlap')]       
        # Open output file
        output_path = os.path.join(base_path, spe, 'bound_random_sv', output_file)
        
        with open(output_path, 'w') as f_out:
            for fi in files:
                try:
                    # Read the overlap file
                    file_path = os.path.join(path, fi)
                    data = pd.read_csv(file_path, sep='\t', header=None)
                    data_filtered = data[data[5] != -1]                    
                    # Count the number of overlapping SVs
                    count = len(data_filtered)
                    f_out.write(f"{count}\n")

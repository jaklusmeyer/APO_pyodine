#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 15:05:14 2021

@author: pheeren
"""

import os


print('----------------------------------------------------')
print('               Welcome to pyodine!                  ')
print()

###############################################################################
## Find out the current path where the pyodine master directory sits
###############################################################################

paths = {
        'i2_dir_path': '',
        'master_dir_path': ''}

print('Setting up pathname structure...')
print()

paths['master_dir_path'] = os.getcwd()

print('Master directory path of pyodine: ', paths['master_dir_path'])
print()

paths['i2_dir_path'] = os.path.join(paths['master_dir_path'], 'iodine_atlas')

print('Iodine directory path of pyodine: ', paths['i2_dir_path'])
print()


###############################################################################
## Adapt the pathnames in:
## - 'utilities.../conf'
## - 'utilities.../pyodine_parameters'
###############################################################################

change_paths = {
        'conf.py': 'i2_dir_path',
        'pyodine_parameters.py': 'master_dir_path'}

print('Aiming to adapt the pathnames in utilities directories...')
print()

utilities_dirs = [os.path.join(paths['master_dir_path'], d) for d in os.listdir() if 'utilities' in d]

print('Utilities directories found:')
for ud in utilities_dirs:
    print('\t{}'.format(ud))

print()

for ud in utilities_dirs:
    for key, value in change_paths.items():
        try:
            with open(os.path.join(ud, key), 'r') as old_file:
                new_file_content = ''
                
                for line in old_file:
                    stripped_line = line.strip()
                    if value in stripped_line and stripped_line[stripped_line.index(value)+1] == '=':
                        stripped_line[-1] = paths[value]
                    new_file_content += stripped_line
            
            with open(os.path.join(ud, key), 'w') as new_file:
                new_file.write(new_file_content)
        
        except Exception as e:
            print(e)
    
    
    
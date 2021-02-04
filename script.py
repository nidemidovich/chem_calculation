# -*- coding: utf-8 -*-
import json
import os
import sys
from subprocess import call

from ase.calculators.gamess_us import GAMESSUS
from ase.io import read
from ase.io import write
from openbabel import pybel

# Optimize calc
opt_calc = GAMESSUS(basis=dict(gbasis='sto', ngauss=3), 
                    contrl=dict(scftyp='rhf', runtyp='optimize'),
                    statpt=dict(opttol=0.0001, nstep=20))

# Energy calc
energy_calc = GAMESSUS(basis=dict(gbasis='sto', ngauss=3),
                       contrl=dict(scftyp='rhf', runtyp='energy'))

def run_gamess(name, 
               input_file, 
               version, 
               ncpus, 
               output_name, 
               input_directory):
	try:
		os.remove(os.path.join('restart', name+'.dat'))
	except:
		pass
    
    rungms = ''
    if sys.platform == 'win32':
        rungms = 'rungms.bat'
    else:
        rungms = 'rungms'
		
	call([rungms, input_file, version, str(ncpus), output_name], 
         shell=True)
	os.rename(output_name, os.path.join(input_directory, output_name))
	print('\n\n')

# Get config
with open('config.json') as f:
    config = json.load(f)
print(config)

if sys.platform == 'win32':
    path_to_gamess = config['path']['Windows']
elif sys.platform == 'linux':
    path_to_gamess = config['path']['Linux']
else:
    path_to_gamess = config['path']['MacOS']

version = config['version']

# Get arguments from console
if len(sys.argv) <= 1:
    input_file = input('Enter path to input: ')
else:
    input_file = sys.argv[1]

try:
    ncpus = sys.argv[2] 
except:
    ncpus = 1
    
name = input_file.split('.')[0]

# smi to mol
smiles = next(pybel.readfile('smi', input_file))
smiles.make3D()
smiles.write(format='mol', filename=name+'.mol', overwrite=True)

# Create inp file for opt
mol = read(name+'.mol')
opt_calc.write_input(atoms=mol)
mol_for_opt_filename = name + '_for_opt.inp'
os.rename('gamess_us.inp', mol_for_opt_filename)
mol_for_opt_name = mol_for_opt_filename.split('.')[0]

# Geometry optimization
output_name = name+'_opt.gamout'
input_directory = os.getcwd()
os.rename(mol_for_opt_filename,
          os.path.join(path_to_gamess, mol_for_opt_filename))
os.chdir(path_to_gamess)

run_gamess(name=mol_for_opt_name,
           input_file=mol_for_opt_filename,
           version=version,
           ncpus=ncpus,
           output_name=output_name,
           input_directory=input_directory)

os.remove(mol_for_opt_filename)
os.chdir(input_directory)

# Get optimized coordinates
mol = next(pybel.readfile('gamout', output_name))
mol.write(format='mol', filename='coords_opt.mol', overwrite=True)

#Create inp file for energy calculation
mol = read('coords_opt.mol')
energy_calc.write_input(atoms=mol)
mol_for_nrg_filename = name + '_for_nrg.inp'
os.rename('gamess_us.inp', mol_for_nrg_filename)
mol_for_nrg_name = mol_for_nrg_filename.split('.')[0]

#Energy calculation
output_name = name+'_nrg.gamout'
input_directory = os.getcwd()
os.rename(mol_for_nrg_filename, 
          os.path.join(path_to_gamess, mol_for_nrg_filename))
os.chdir(path_to_gamess)

run_gamess(name=mol_for_nrg_name,
           input_file=mol_for_nrg_filename,
           version=version,
           ncpus=ncpus,
           output_name=output_name,
           input_directory=input_directory)

os.remove(mol_for_nrg_filename)
os.chdir(input_directory)

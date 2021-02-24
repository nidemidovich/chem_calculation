import os
import sys
from subprocess import call

from ase.calculators.gamess_us import GAMESSUS
from ase.io import read
from openbabel import pybel

class Job:
    
    # Optimize calc
    opt_calc = GAMESSUS(basis=dict(gbasis='sto', ngauss=3), 
                        contrl=dict(scftyp='rhf', runtyp='optimize'),
                        statpt=dict(opttol=0.0001, nstep=20))

    # Energy calc
    energy_calc = GAMESSUS(basis=dict(gbasis='sto', ngauss=3),
                           contrl=dict(scftyp='rhf', runtyp='energy'))
    
    @staticmethod
    def run_gamess(name, input_file, version, 
                   ncpus, output_name, input_directory):
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
    
    def __init__(self, i, smiles, ncpus, version, path_to_gamess):
        self.smiles = smiles.rstrip()
        self.ncpus = ncpus
        self.version = version
        self.path_to_gamess = path_to_gamess
        self.name = str(i)
        
    def smi_to_mol(self):
        smiles = pybel.readstring('smi', self.smiles)
        smiles.make3D()
        smiles.write(format='mol', filename=self.name+'.mol', overwrite=True)
        
    def create_inp_file(self):
        mol = read(self.name+'.mol')
        Job.opt_calc.write_input(atoms=mol)
        self.mol_for_opt_filename = self.name + '_for_opt.inp'
        os.rename('gamess_us.inp', self.mol_for_opt_filename)
        self.mol_for_opt_name = os.path.splitext(self.mol_for_opt_filename)[0]
    
    def geometry_optimization(self):
        self.output_name = self.name+'_opt.gamout'
        input_directory = os.getcwd()
        os.rename(self.mol_for_opt_filename,
                  os.path.join(self.path_to_gamess, self.mol_for_opt_filename))
        os.chdir(self.path_to_gamess)
        print(self.mol_for_opt_filename, os.getcwd())

        Job.run_gamess(name=self.mol_for_opt_name,
                       input_file=self.mol_for_opt_filename,
                       version=self.version,
                       ncpus=self.ncpus,
                       output_name=self.output_name,
                       input_directory=input_directory)
        
        os.remove(self.mol_for_opt_filename)
        os.chdir(input_directory)
        
    def get_optimized_coords(self):
        mol = next(pybel.readfile('gamout', self.output_name))
        mol.write(format='mol', filename='coords_opt.mol', overwrite=True)
        
    def create_inp_file_for_nrg_calc(self):
        mol = read('coords_opt.mol')
        Job.energy_calc.write_input(atoms=mol)
        self.mol_for_nrg_filename = self.name + '_for_nrg.inp'
        os.rename('gamess_us.inp', self.mol_for_nrg_filename)
        self.mol_for_nrg_name = os.path.splitext(self.mol_for_nrg_filename)[0]
        os.remove('coords_opt.mol')
        
    def nrg_calc(self):
        self.output_name = self.name+'_nrg.gamout'
        input_directory = os.getcwd()
        os.rename(self.mol_for_nrg_filename, 
                  os.path.join(self.path_to_gamess, self.mol_for_nrg_filename))
        os.chdir(self.path_to_gamess)

        Job.run_gamess(name=self.mol_for_nrg_name,
                       input_file=self.mol_for_nrg_filename,
                       version=self.version,
                       ncpus=self.ncpus,
                       output_name=self.output_name,
                       input_directory=input_directory)

        os.remove(self.mol_for_nrg_filename)
        os.chdir(input_directory)
        
    def run(self):
        self.smi_to_mol()
        self.create_inp_file()
        self.geometry_optimization()
        self.get_optimized_coords()
        self.create_inp_file_for_nrg_calc()
        self.nrg_calc()
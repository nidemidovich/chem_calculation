import glob
import os
import sys
from subprocess import call

from ase.calculators.gamess_us import GAMESSUS
from ase.io import read
from openbabel import pybel

class Job:
    
    i = 0
    
    @staticmethod
    def run_gamess(name, input_file, temp, path,
                   version, ncpus, output_name):
        cwd = os.getcwd()
        os.chdir(os.path.expanduser('~/'))
        files = glob.glob('gamess-files/'+temp+'/'+name+'*')
        for f in files:		
            try:
                os.remove(f)
            except:
                pass
        os.chdir(cwd)
        call('{0} {1} {2} {3} {4}>{5}'.format(path, input_file, version, ncpus, temp, output_name), 
             shell=True)
        print('\n\n')
    
    def __init__(self, smiles, temp, ncpus, version, path, iterations):
        self.smiles = smiles.rstrip()
        self.temp = temp
        self.ncpus = ncpus
        self.version = version
        self.path = path
        Job.i += 1
        self.name = str(Job.i)
        self.iterations = iterations
        
    def smi_to_mol(self):
        smiles = pybel.readstring('smi', self.smiles)
        smiles.make3D(steps=int(self.iterations))
        smiles.write(format='mol', filename=self.name+'.mol', overwrite=True)
        
    def create_inp_file(self):
        mol = read(self.name+'.mol')
        opt_calc = GAMESSUS(basis=dict(gbasis='sto', ngauss=3), 
                            contrl=dict(scftyp='rhf', runtyp='optimize'),
                            statpt=dict(opttol=0.0001, nstep=20))
        opt_calc.write_input(atoms=mol)
        self.mol_for_opt_filename = self.name + '_for_opt.inp'
        os.rename('gamess_us.inp', self.mol_for_opt_filename)
        self.mol_for_opt_name = os.path.splitext(self.mol_for_opt_filename)[0]
    
    def geometry_optimization(self):
        self.output_name = self.name+'_opt.gamout'
        Job.run_gamess(name=self.mol_for_opt_name,
                       input_file=self.mol_for_opt_filename,
                       temp=self.temp,
                       path=self.path,
                       version=self.version,
                       ncpus=self.ncpus,
                       output_name=self.output_name)
        
    def get_optimized_coords(self):
        mol = next(pybel.readfile('gamout', self.output_name))
        mol.write(format='mol', filename='coords_opt.mol', overwrite=True)
        
    def create_inp_file_for_nrg_calc(self):
        mol = read('coords_opt.mol')
        energy_calc = GAMESSUS(basis=dict(gbasis='sto', ngauss=3),
                               contrl=dict(scftyp='rhf', runtyp='energy'))
        energy_calc.write_input(atoms=mol)
        self.mol_for_nrg_filename = self.name + '_for_nrg.inp'
        os.rename('gamess_us.inp', self.mol_for_nrg_filename)
        self.mol_for_nrg_name = os.path.splitext(self.mol_for_nrg_filename)[0]
        os.remove('coords_opt.mol')
        
    def nrg_calc(self):
        self.output_name = self.name+'_nrg.gamout'
        Job.run_gamess(name=self.mol_for_nrg_name,
                       input_file=self.mol_for_nrg_filename,
                       temp=self.temp,
                       path=self.path,
                       version=self.version,
                       ncpus=self.ncpus,
                       output_name=self.output_name)
        
    def run(self):
        self.smi_to_mol()
        self.create_inp_file()
        self.geometry_optimization()
        self.get_optimized_coords()
        self.create_inp_file_for_nrg_calc()
        self.nrg_calc()
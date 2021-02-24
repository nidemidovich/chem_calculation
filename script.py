import sys
import json

from calc import Job

# Get config
with open('config.json') as f:
	config = json.load(f)

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

with open(input_file) as f:
    i = 0
    for smiles in f:
        i += 1
        print(smiles)
        job = Job(i, smiles, ncpus, version, path_to_gamess)
        job.run()
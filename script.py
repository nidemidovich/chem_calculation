import sys
import json

from calc import Job

# Get config
with open('config.json') as f:
	config = json.load(f)

path = config['path']
version = config['version']

# Get arguments from console
if len(sys.argv) <= 1:
	input_file = input('Enter path to input: ')
else:
	input_file = sys.argv[1]

try:
	temp = sys.argv[2]
except:
	temp = 'temp'

try:
	ncpus = sys.argv[3] 
except:
	ncpus = 1

with open(input_file) as f:
    for smiles in f:
        print(smiles)
        job = Job(smiles, temp, ncpus, version, path)
        job.run()
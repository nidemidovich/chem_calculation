import multiprocessing
import sys
import json

from calc import Job

# Get config
with open('config.json') as f:
	config = json.load(f)

path = config['path']
version = config['version']
iterations = config['iterations']

# Get arguments from console
input_file = sys.argv[1]

try:
	temp = sys.argv[2]
except:
	temp = 'temp'

try:
	ncpus = sys.argv[3] 
except:
	ncpus = '1'

processes = []
with open(input_file) as f:
	for smiles in f:
		print(smiles)
		job = Job(smiles, temp, ncpus, version, path, iterations)
		p = multiprocessing.Process(target=job.run)
		p.start()
		processes.append(p)

for process in processes:
	process.join()
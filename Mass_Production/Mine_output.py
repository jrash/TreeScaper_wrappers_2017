#! /usr/bin/env python

import os
import glob
import shutil 
import itertools
import re 

# Open file for recording output
outOutFile = open('Compiled_info.out', 'w' )

# Pattern to find in file
pattern = 'Tacos'

# Header. Eventually pull from file, in case order changes. For now I'm lazy.
header = "Eat,Name,Tips,Trees,Cloud_RF,Network,Model,Communities,Plateau_low,Plateau_highq,Lambda,Time,Rooted,Weighted,Fixed_Lambda_Cov,High_Freq,Low_Freq,Fixed_Lambda_Aff,Distance_Metric,Affinity_Transformation\n"

outOutFile.write(header)

# Make sure you are in directory that you are calling script from 
mainDir = os.getcwd()

# Iterate through folders
for t in glob.glob('*.nex'):
	# Split off name, dependant on what your iterating handle files are
	repName = t.split(".")[0]
	outFile = repName+".out"
	# Create outfile path 
	outFilePath = os.path.join(mainDir,repName,outFile)
	with open(outFilePath, "r") as file:
		for line in file:
			if pattern in line:
				outOutFile.write(line)

outOutFile.close()
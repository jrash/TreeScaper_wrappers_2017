#! /usr/bin/env python

import os
import glob
import shutil 
import itertools
import re 

pbs = "job_runCLV.pbs"
# Make sure you are in directory that you are calling script from 
mainDir = os.getcwd()

for t in glob.glob('*.nex'):
	# Split off name, dependant on what your iterating handle files are
	repName = t.split(".")[0]
	# Create path 
	dirPath = os.path.join(mainDir,repName)
	os.chdir(dirPath)
	newPbs = str(pbs.split(".")[0])+"_"+repName+".pbs"
	os.system("qsub %s" % (newPbs))
	os.chdir(mainDir)
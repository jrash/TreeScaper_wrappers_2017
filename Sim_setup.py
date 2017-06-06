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
	# Make directory 
	if not os.path.exists(dirPath):
		os.mkdir(dirPath)	
	# Copy nexus to folder
	os.system("cp %s %s" % (t,dirPath))
	# Copy pbs to folder
	newPbs = str(pbs.split(".")[0])+"_"+repName+".pbs"
	os.system("cp %s %s" % (pbs,dirPath+"/"+newPbs))
	os.chdir(dirPath)
	with open(newPbs, "r+") as f:
			filedata = f.read()
			filedata = re.sub("tacocat",str(repName), filedata)
			f.seek(0)
			f.write(filedata)
			f.truncate()
	os.chdir(mainDir)

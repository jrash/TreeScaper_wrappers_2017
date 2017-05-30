#!/usr/bin/env python

# Eventually want to just add this to the other wrapper python file. 
# for now, just find a way to output the right # 
#Usage: autoFindlambda.py
#Inputs:
# - automatic output file. Ex: all_trees_CovAuto.out
# - probably want to specify cov vs aff output files?

from __future__ import division
import numpy



inOutFile = 'all_trees_CovAuto.out'


def autoFindlambda(inOutFile):
	largePlateau=[]
	allPlateaus=[]
	pattern = 'The found plateaus are:'
	with open(inOutFile , 'r' ) as file:
		for line in file:
			if pattern in line:
				# This pulls both the largest plateau, and a line containing all lambdas
				# The second occurance isnt parse properly in this script. But keeping it for parsing later if needed.
				a = (next(file,'').strip())
				b = a.split(",")[0].strip("[")
				c = a.split(",")[1].strip("]")
				# Store both lines for later use
				allPlateaus.append(a)
				# Add first and second numbers
				largePlateau.append(b)
				largePlateau.append(c)
	# Pull out largest plateau and get lambda in middle of it. 
	plat = [float(largePlateau[0]), float(largePlateau[1])]
	manLamdba = numpy.mean(plat)
	return manLamdba

'''
def main():
	outFile = sys.argv[1]

if __name__=='__main__':
	main()
'''
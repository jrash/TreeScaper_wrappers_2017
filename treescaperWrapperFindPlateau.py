#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Usage: treescaperWrapperFindPlateau.py [CLV path] [model] [network] [rooted]
#[model] can be CNM/CPM/ERNM/NNM

#***Check the ./CLVTreescaper settings within the script.  You might want to run ./CLVTreescaper with different options.  
# Make sure the numbering of all the indices is correct. The developers have changed between starting at 1 and starting at 0.  
# This includes the affinityCommunities script that is imported into this script.

#useful if you have found the plateau with the automatic search function of the treescaper GUI.  If you enter the lambda values where the plateau was found for both affinity and covariance matrices, you will get all the useful output of treescaperWrapperV2.py.  See treescaperWrapperV2.py for usage and output.

#output files

# ['type'] can be Covariance or Affinity
# ['treeset'] tree set name

#['treeset']_['type']WholeCommunity_results.out: community results over the whole range of lambda values
#['treeset']_CovPlateauCommunity.out: community structure of the plateau
#['treeset']_comKey.out: key showing you which bipartitions are in which communities
# AffinityCom[number].nex: a nexus file of the trees in an affinity community
# AffinityCom[number].nex.con: consensus tree of an affinity community
# AffinityCom[number].nex.con.pdf: pdf of consensus tree of an affinity community

import re
import os
import sys
import glob
import numpy as np
import fnmatch
from collections import Counter
import filecmp
from AffinityCommunities import affinityCommunityConsensus
import dendropy
import StringIO
import heapq
from shutil import move, copyfile
from tempfile import mkstemp


def make_list(pre_ls,convert):
	# Turns tab delimited line into a list. 
	pre_ls = pre_ls.split("\t")
	ls = [i for i in pre_ls if i != "\n"]
	if convert == "int":
		ls = [int(i) for i in ls]
	if convert == "float":
		ls = [float(i) for i in ls]
	#ls = ls[1:]
	#print("list - "+str(ls))
	return ls


def reg_ex_match(file, pattern):
	# Returns the first match of a reg ex search
	file.seek(0)
	for line in file:
		m = pattern.match(line)
		if m:
			return m.group(1)

def autoFindlambda(inOutFile):
	# Currently pulls out largest plateau found in auto run. 
	# Would be useful to also pull out all possibly plateaus. 
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
	manLamdba = np.mean(plat)
	return manLamdba


def mode_function(lst):
	# Returns a list of all the possible plateaus
	# Not sure what input is. 

	counterLst = Counter(lst)
	_,val = counterLst.most_common(1)[0]
	modeLS = [x for x,y in counterLst.items() if y == val]
	lstNoMode = [x for x in lst if x not in modeLS] #remove all the values that are the mode of the list
	if len(counterLst.most_common(2)) == 2: # find the next mode of the list
		counterLstNoMode = Counter(lstNoMode)
		_,val2 = counterLstNoMode.most_common(1)[0]
		if val - val2 >= 2:
			#if the largest plateau is bigger than the second largest by one increment, than its possible that the true size of the second largest
			#plateau is bigger than the largest. The second largest plateau could encompass the space between the lambda value immediately below and immediately above its lambda range.  The second largest plateau could
			#grow by 1.9999... increment and largest plateau could not grow at all. But if the largest plateau is bigger than the second largest by two or more increments, you are sure that plateau is the largest
			#****newly added feature, not used for any of my results
			return modeLS
		else:
			return [x for x,y in counterLst.items() if y == val or y == val2]
	else:
		return modeLS
	
	
def get_plateau(clvPath, treeSet, treeSetTrunc, type, model, rooted, plateauLambda):
	# Get info shit and output some files

	if type == "Covariance":

		# Open out file from manual run
		comFile = open("%s_CovCommunity.out" %  treeSet, 'r')
		pattern = re.compile('Number of communities: (\d+)')
		coms = int(reg_ex_match(comFile, pattern))
		comKey = open("%s_comKey.out" %  treeSetTrunc, 'w')

		# Makes a key to decode the bipartitions in each community
		for i in range(1, coms+1): 
			# Get nodes for each community
			pattern = re.compile('Community '+str(i)+' includes nodes: (.+)')
			comStr = reg_ex_match(comFile, pattern)
			comLS = comStr.split(",")
			print('Community '+str(i)+' includes nodes: '+str(comLS))
			comLS = filter(None, comLS)
			comKey.write("Com %s:\n" % i)
			# For each node
			for j in comLS:
				# Bipartitions count from 1 while nodes count from 0.
				j = int(j)+1
				pattern = re.compile("bipartition "+str(j)+" : ([0-1]+?), appear times: ([0-9]+?)$")
				# Some file copying, temporary
				fh, absPath = mkstemp()
				copyfile("%s_CovCommunity.out" %  treeSet, absPath)
				comTempFile = open(absPath,'r')
				# Grab bipartitions
				for line in comTempFile:
					m = pattern.match(line)
					if m:
						# Grab binary code for bipartition
						bipart = m.group(1)
						print("Bipartion "+str(j)+" : "+str(bipart))
						# Grab number of times bipartition occurs
						freq = m.group(2)
						bipartLS = []
						# Create list from binary code
						for c in bipart:
							bipartLS.append(c)
						# Create a list of taxa included in each bipartition
						indices = [x+1 for x, y in enumerate(bipartLS) if y == '1']
						# Get taxon names
						for k in indices:
							pattern2 = re.compile("(.+) , "+str(k))
							taxon = reg_ex_match(comFile, pattern2)
							indices[indices.index(k)] = taxon
				# Close temp file
				comTempFile.close()
				os.close(fh)
				# Write taxon and frequency for bipartitions in communit
				comKey.write("%s %s\n" % (indices,freq))
			comKey.write("\n")
		comKey.close()
#
	if type == "Affinity":
		os.system("%s -trees -f %s -w 0 -r %s -o Community -t Affinity -cm %s -lm manu -dm URF -am Exp -lp %s -ln 1 " % (clvPath, treeSet, rooted, model, plateauLambda)+\
		" > %s_AffPlateauCommunity.out" %  treeSetTrunc)#outputs community structure for current lambda values

	return plateauLambda


def main():
	clvPath = sys.argv[1]
	inNexus = sys.argv[2]
	model = sys.argv[3]
	network = sys.argv[4]
	rooted = sys.argv[5]

	# See treescaperWrapperKnownPlateau.py for a script that can add a translate block via dendropy. Not used here because it was glitching on supermike. 

	treeSet=str(inNexus)
	treeSetIndex = treeSet.find(".")
	treeSetTrunc = treeSet[:treeSetIndex]
	os.system("echo 'hi'")
	
	if network == 'Covariance':

		# Run automatic plateau finder
		print("Running manual with lambda = %s. Log file: %s_CovCommunity.out" %  (plateau, treeSet))
		os.system("%s -trees -f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm auto -hf .95 -lf .05" % (clvPath, treeSet, rooted, model)+\
	 	" > %s_CovAuto.out" %  treeSet)

	 	# Get middle of largest plateau
		outFile = "%s_CovAuto.out" %  treeSet
		plateau = autoFindlambda(outFile)

		# Run manual plateau
		# Outputs community structure for current lambda values
		print("Running manual with lambda = %s. Log file: %s_CovCommunity.out" %  plateau, treeSet)
	 	os.system("%s -trees -f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm manu -lp %s -ln 1 -hf .95 -lf .05" % (clvPath, treeSet, rooted, model, plateau)+\
	 	" > %s_CovCommunity.out" %  treeSet)

	 	# Get output file from automatic run
	 	cmCar=glob.glob('%s*_Covariance_Matrix_*community_auto_results.out' % (treeSetTrunc))

	 	# Change name, might want to turn this into cp instead of mv
	 	os.system("mv %s %s_CovWholeCommunity_results.out" % (str(cmCar[0]), treeSetTrunc))
	 	print("Create extra output files")
	 	plateauLambda = get_plateau(clvPath, treeSet, treeSetTrunc, "Covariance", model, rooted, plateau)

	 	print("platLambd"+str(plateauLambda))

	if network == 'Affinity':

		# Run automatic plateau finder
		os.system("%s -trees -f %s -ft Trees -w 0 -r %s -o Community -t Affinity -cm %s -lm auto -dm URF -am Exp" % (clvPath, treeSet, rooted, model)+\
	 	" > %s_AffAuto.out" %  treeSet)

	 	# Get middle of largest plateau
		outFile = "%s_AffAuto.out" %  treeSet
		plateau = autoFindlambda(outFile)

		# Re-run Treescaper with manual plateau

		os.system("%s -trees -f %s -w 0 -r %s -o Community -t Affinity -cm %s -lm manu -dm URF -am Exp -lp %s -ln 0 " % (clvPath, treeSet, rooted, model, plateau)+\
		" > %s_AffCommunity.out" %  treeSet)

		print("affinity plateau "+str(plateau))

		# change file name. 
		aCar=glob.glob('%s*_Affinity-*community_auto_results.out' % (treeSetTrunc))

		os.system("mv %s %s_AffWholeCommunity_results.out" % (str(aCar[0]), treeSetTrunc))
#
		plateauLambda = get_plateau(clvPath, treeSet, treeSetTrunc, "Affinity", model, rooted, plateau)
		print("platLambd "+str(plateauLambda))
#
		if plateauLambda:
			affinityCommunityConsensus(clvPath, treeSet, model, plateauLambda, rooted)
	print(type(treeSetTrunc))
 	os.system("sumtrees.py -r -o %s.con %s" % (treeSetTrunc, inNexus))
	os.system("cat ./SeqSim/FigTreeBlock.txt >> %s.con" % (treeSetTrunc))
 	#os.system("/Applications/FigTree/FigTree_v1.4.3/bin/figtree -graphic PDF all_trees.con all_trees.pdf")

	os.system("echo 'treescaperWrapperKnownPlateau.py %s\n' >> commands.txt" % ' '.join(sys.argv[1:]))


	#os.system("cat RAxML_bestTree.G1_gene* > RAxML_allTree.G1.tre")
if __name__=='__main__':
	main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

#### Input

# Usage: treescaperWrapperFindPlateau.py [CLV path] [inputTrees.nex] [model] [network] [rooted] 
# [model] can be CNM/CPM/ERNM/NNM
# [network] can be Covariance/Affinity
# See main() for additional CLVTreeScaper options


#### Output 

#['inputTrees']_['network']WholeCommunity_results.out: community results over the whole range of lambda values when auto plateau is run.
#['inputTrees']_CovCommunity.out: community structure of the plateau with fixed lambda
#['inputTrees']_comKey.out: key showing you which bipartitions are in which communities
# AffinityCom[number].nex: a nexus file of the trees in an affinity community
# AffinityCom[number].nex.con: consensus tree of an affinity community
# AffinityCom[number].nex.con.pdf: pdf of consensus tree of an affinity community
# There are additional output files not described here yet. Please feel free to fill in! 

#### Additional Notes
# Double check indexing of nexus file vs trees assigned to Affinity communities. This script should be able to handle nexus counting from 0 or from 1. Will break if Treescaper counts from 1. 
# affinityCommunities script has been incorperated into this script


import re
import os
import sys
import glob
import numpy as np
import fnmatch
import time
from collections import Counter
import filecmp
import dendropy
import StringIO
import heapq
from shutil import move, copyfile
from tempfile import mkstemp

startTime = time.time()

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

def edit_treeset(treeFileEditPath):
	# Add comment blocks that number each tree with the indices used by TreeScaper. Pulled from AffinityCommunities.py

	treeFileEdit = open(treeFileEditPath, 'r')

	lineNum = 1
	# make a temp file
	fh, absPath = mkstemp()
	tempFile = open(absPath,'w')
	for line in treeFileEdit:
		if line.find('[&U]') != -1: # Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
			if line.find('['+str(lineNum)+']') != -1:
				tempFile.write(line)
			else:
 				tempFile.write(line.replace('=','['+str(lineNum)+']='))
			lineNum += 1
		elif line.find('[&R]') != -1: # Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
			if line.find('['+str(lineNum)+']') != -1:
				tempFile.write(line)
			else:
 				tempFile.write(line.replace('=','['+str(lineNum)+']='))
			lineNum += 1
		else:
			tempFile.write(line)
	# close temp file
	tempFile.close()
	os.close(fh)
	treeFileEdit.close()
	# Remove original file
	os.remove(treeFileEditPath)
	# Move new file
	move(absPath, treeFileEditPath)
	return lineNum

def affinityCommunityConsensus(clvPath,treeFile,model,plateau,rooted):
	# Best function ever
	# Parses communities and associate trees. Creates a nexus file and consensus tree for each community. 
	
	if model not in ('CPM','ERNM','CNM','NNM'):
		print "invalid model choose CPM, ERNM, CNM, or NNM"
		model = raw_input('Enter Model: ')

	print("Plateau lambda: "+str(plateau))

	# Get number of communities from manual output file
	comFile = open( "%s_AffCommunity.out" % (treeFile) , 'r' )
	pattern = re.compile('Number of communities: (\d+)')
	coms = int(reg_ex_match(comFile, pattern))
	print("Number of communities: "+str(coms))

	totalTrees = edit_treeset(treeFile)
	treeFile = open(treeFile,'r')

	# Make a file to count the frequency and relative frequency of trees in each community
	treeCountFile = open("AffinityCommunitiesTreeCount.txt", 'w')

	# Get the tree indices from each community
	for i in range(1, coms+1):
		pattern = re.compile('Community '+str(i)+' includes nodes: (.+)')
		comStr = reg_ex_match(comFile, pattern)
		comLs = comStr.split(",")
		print("Trees in community "+str(i)+" : "+str(comLs))
		comLs = filter(None, comLs)
		# Pull out these trees from the original trees. Make a nexus file containing the trees for each community
		comTreeSetStr = 'AffinityCom'+str(i)+'.nex'
		comTreeSet = open(comTreeSetStr,'w')
		treeCount = 0
		treeFile.seek(0)
		for line in treeFile:
			# Print translate block to new community#.nex file
			if line.find('[&U]') == -1 & line.find('[&R]') == -1:
				comTreeSet.write(line)
			else:
				if line.find('[0]') != -1:
					# Assuming that Treescaper counts from 0 and Nexus file counts from 0.
					for j in comLs:
						#print("j: "+str(j))
						if line.find('['+str(j)+']') != -1:
							print("line: "+str(line))
							comTreeSet.write(line)
							treeCount += 1
				else:
					# Adjust for Nexus file counting from 1. 
					for j in comLs:
						k = int(j) + 1
						#print("j: "+str(j))
						#print("k: "+str(k))
						if line.find('['+str(k)+']') != -1:
							#print("line: "+str(line))
							comTreeSet.write(line)
							treeCount += 1
		# Write frequency and relative frequency of trees in the community
		treeCountFile.write("%s\t%s\t%.2f%% of trees\n" % (comTreeSetStr, str(treeCount),100*(treeCount/totalTrees)))
		comTreeSet.close()
		# Make a consensus tree of the affinity community
		comTreeConStr = comTreeSetStr + ".con"
		if rooted == '0':
			os.system("sumtrees.py -r --unrooted -o %s %s" % (comTreeConStr,comTreeSetStr))
		if rooted == '1':
			os.system("sumtrees.py -r --rooted -o %s %s &> dendropy_%s.out" % (comTreeConStr,comTreeSetStr,comTreeSetStr))
		os.system("cat ./SeqSim/FigTreeBlock.txt >> %s" % (comTreeConStr))
		#Make a pdf of the consensus tree
		#os.system("figtree -graphic PDF %s %s.pdf" % (comTreeConStr, comTreeConStr))
	
def parse_output(clvPath, treeSet, treeSetTrunc, type, model, rooted, plateauLambda):
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
	if type == "Affinity":
		affinityCommunityConsensus(clvPath, treeSet, model, plateauLambda, rooted)


def main():
	clvPath = sys.argv[1]
	inNexus = sys.argv[2]
	model = sys.argv[3]
	network = sys.argv[4]
	rooted = sys.argv[5]
	startTime = time.time()

	# Manually edit below values if needed. 
	w = 0
	# Covariance
	ln_c = 1
	hf = 0.95
	lf = 0.05
	# Affinity
	ln_a = 0
	dm = "URF"
	am = "Exp"
	# Info:
	#w: Indicate whether trees are weighted. Options are: ‘1’: weighted. ‘0’: 	unweighted
	#ln: Specify a fixed value of λ−. Must be between 0 and 1. Used when -lniv is zero.
	#hf: Frequency upper bound. A number between 0 and 1. Nodes with frequencies above this value are ignored.
	#lf: Frequency lower bound. A number between 0 and 1. Nodes with frequencies below this value are ignored.
	#dm: Indicates the distance metric. Options are: ‘URF’: Unweighted 	Robinson-Foulds distance. ‘RF’: Weighted Robinson-Foulds distance. ‘Mat’: 	Matching distance. ‘SPR’: Subtree-Prune-Regraft
	#am: Indicates the distance to affinity transformation. Options are: ‘Rec’: Reciprocal. ‘Exp’: Exponential

	# See treescaperWrapperKnownPlateau.py for a script that can add a translate block via dendropy. Not used here because it was glitching on supermike. 

	# Get some strings for naming files
	treeSet=str(inNexus)
	treeSetIndex = treeSet.find(".")
	treeSetTrunc = treeSet[:treeSetIndex]

	# Open file for keeping track of time
	timeFile = open( "%s_%s_Time.out" % (inNexus, network) , 'w' )


	os.system("echo 'Hello There'")
	
	if network == 'Covariance':

		# Run automatic plateau finder
		print("Running automatic. Log file: %s_CovAuto.out" %  treeSet)
		startTime2 = time.time()
		os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Covariance -cm %s -lm auto -hf %s -lf %s" % (clvPath, treeSet, w, rooted, model, hf, lf)+\
	 	" > %s_CovAuto.out" %  treeSet)
		endTime2 = time.time()

	 	# Get middle of largest plateau
		outFile = "%s_CovAuto.out" %  treeSet
		plateau = autoFindlambda(outFile)

		# Run manual plateau
		# Outputs community structure for current lambda values
		print("Running manual with lambda = %s. Log file: %s_CovCommunity.out" %  (plateau, treeSet))
		startTime1 = time.time()
	 	os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Covariance -cm %s -lm manu -lp %s -ln %s -hf %s -lf %s" % (clvPath, treeSet, w, rooted, model, plateau, ln_c, hf, lf)+\
	 	" > %s_CovCommunity.out" %  treeSet)
	 	endTime1 = time.time()
	 	# Get output file from automatic run
	 	cmCar=glob.glob('%s*_Covariance_Matrix_*community_auto_results.out' % (treeSetTrunc))

	 	# Change name, might want to turn this into cp instead of mv
	 	os.system("mv %s %s_CovWholeCommunity_results.out" % (str(cmCar[0]), treeSetTrunc))

	 	print("Parse output into useful information")
	 	startTime3 = time.time()
	 	parse_output(clvPath, treeSet, treeSetTrunc, "Covariance", model, rooted, plateau)
	 	endTime3 = time.time()

	if network == 'Affinity':

		# Run automatic plateau finder
		print("Running automatic. Log file: %s_AffAuto.out" %  treeSet)
		startTime2 = time.time()
		os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Affinity -cm %s -lm auto -dm %s -am %s" % (clvPath, treeSet, w, rooted, model, dm, am)+\
		" > %s_AffAuto.out" %  treeSet)
		endTime2 = time.time()

		# Get middle of largest plateau
		outFile = "%s_AffAuto.out" %  treeSet
		plateau = autoFindlambda(outFile)

		# Run Treescaper with manual plateau
		print("Running manual with lambda = %s. Log file: %s_AffCommunity.out" %  (plateau, treeSet))
		startTime1 = time.time()
		os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Affinity -cm %s -lm manu -dm %s -am %s -lp %s -ln %s " % (clvPath, treeSet, w, rooted, model, dm, am, plateau, ln_a)+\
		" > %s_AffCommunity.out" %  treeSet)
		endTime1 = time.time()

		# Get output file from automatic run
		aCar=glob.glob('%s*_Affinity-*community_auto_results.out' % (treeSetTrunc))

		# Change name, might want to turn this into cp instead of mv
		os.system("mv %s %s_AffWholeCommunity_results.out" % (str(aCar[0]), treeSetTrunc))

		print("Parse output into useful information")
		startTime3 = time.time()
		parse_output(clvPath, treeSet, treeSetTrunc, "Affinity", model, rooted, plateau)
		endTime3 = time.time()


	# If you've screwed something up...

	if (network != 'Affinity') and (network != 'Covariance'):
		print("Check spelling and order of input values")


	# Make consensus tree for inNexus
	print("Building consensus tree for input file. Log file: dendropy_%s.out" %  treeSetTrunc)
 	startTime4 = time.time()
 	os.system("sumtrees.py -r -o %s.con %s &> dendropy_%s.out" % (treeSetTrunc, inNexus,inNexus))
 	endTime4 = time.time()

 	# Add FigTree block to file, make PDF
	os.system("cat ./SeqSim/FigTreeBlock.txt >> %s.con" % (treeSetTrunc))
 	os.system("./bin/figtree -graphic PDF all_trees.con all_trees.pdf")
 	#os.system("cat RAxML_bestTree.G1_gene* > RAxML_allTree.G1.tre")

	os.system("echo 'treescaperWrapperKnownPlateau.py %s\n' >> commands.txt" % ' '.join(sys.argv[1:]))

	endTime = time.time()

	# Calculate time for each section, write to file. 

	time1 = "Manual_CLVTreeScaper:\t"+str(round(endTime1 - startTime1, 5))
	time2 = "Automatic_CLVTreeScaper:\t"+str(round(endTime2 - startTime2, 5))
	time3 = "File_Parsing:\t"+str(round(endTime3 - startTime3, 5))
	time4 = "SumTree_inNexus:\t"+str(round(endTime4 - startTime4, 5))
	timeAll = "Total_time:\t"+str(round(endTime - startTime, 5))
	
	timeFile.write("%s\n%s\n%s\n%s\n%s" % (time1, time2, time3, time4, timeAll))
	timeFile.close()

if __name__=='__main__':
	main()

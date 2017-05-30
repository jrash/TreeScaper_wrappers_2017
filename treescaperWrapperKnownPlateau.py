#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Usage: treescaperWrapperKnownPlateau.py [CLV path] [model] [plateau] [network] [rooted]
#[model] can be CNM/CPM/ERNM/NNM

#***Check the ./CLVTreescaper settings within the script.  You might want to run ./CLVTreescaper with different options.  
# Make sure the numbering of all the indices is correct. The developers have changed between starting at 1 and starting at 0.  
# This includes the affinityCommunities script that is imported into this script.

#useful if you have found the plateau with the automatic search function of the treescaper GUI.  If you enter the lambda values where the plateau was found for both affinity and covariance matrices, you will get all the useful output of treescaperWrapperV2.py.  See treescaperWrapperV2.py for usage and output.

#output files

# ['type'] can be Covariance or Affinity
# ['treeset'] tree set name
# ['plateauLambda'] the lambda value were the plateau was found

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
import dendropy
import StringIO
import heapq
from shutil import move, copyfile
from tempfile import mkstemp


def make_list(pre_ls,convert):
	pre_ls = pre_ls.split("\t")
	ls = [i for i in pre_ls if i != "\n"]
	if convert == "int":
		ls = [int(i) for i in ls]
	if convert == "float":
		ls = [float(i) for i in ls]
	# ls = ls[1:]
	# print("list - "+str(ls))
	return ls

def reg_ex_match(file, pattern):
	# Returns the first match of a reg ex search

	file.seek(0)
	for line in file:
		m = pattern.match(line)
		if m:
			return m.group(1)

def edit_treeset(treeFileEditPath):
	# Add comment blocks that number each tree with the indices used by TreeScaper. Pulled from AffinityCommunities.py

	treeFileEdit = open(treeFileEditPath, 'r')

	lineNum = 1
	# make a temp file
	fh, absPath = mkstemp()
	tempFile = open(absPath,'w')
	for line in treeFileEdit:
		if line.find('[&U]') != -1: # Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
 			tempFile.write(line.replace('=','['+str(lineNum)+']='))
			lineNum += 1
		elif line.find('[&R]') != -1: # Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
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
	plateau = sys.argv[4]
	network = sys.argv[5]
	rooted = sys.argv[6]

	'''
	# Use to input a translate block with dendropy. The wrong version of dendropy actually deletes the translate block. So use careully. 

	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, inNexus):
			if rooted == '1':
				tlst = dendropy.TreeList.get_from_path(inNexus, "nexus", rooting='force-rooted')
			else:
				tlst = dendropy.TreeList.get_from_path(inNexus, "nexus", rooting='force-unrooted')

	tlst.write_to_path("all_trees.pre.nex",'nexus', simple=True, translate_tree_taxa=True)
	
	with open("all_trees.nex", "w") as fout:
		with open("all_trees.pre.nex", "r") as fin:
			for line in fin:
				fout.write(re.sub('END;\n', 'END;',line))

	os.system("rm all_trees.pre.nex")
	'''

	treeSet=str(inNexus)
	treeSetIndex = treeSet.find(".")
	treeSetTrunc = treeSet[:treeSetIndex]
	os.system("echo 'hi'")
	
	if network == 'Covariance':

		# Run manual plateau
		# Outputs community structure for current lambda values
		print("Running manual with lambda = %s. Log file: %s_CovCommunity.out" %  (plateau, treeSet))
	 	os.system("%s -trees -f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm manu -lp %s -ln 1 -hf .95 -lf .05" % (clvPath, treeSet, rooted, model, plateau)+\
	 	" > %s_CovCommunity.out" %  treeSet)

	 	# Run automatic plateau finder
	 	print("Running automatic. Log file: %s_CovAuto.out" %  treeSet)
		os.system("%s -trees -f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm auto -hf .95 -lf .05" % (clvPath, treeSet, rooted, model)+\
	 	" > %s_CovAuto.out" %  treeSet)

	 	# Get output file from automatic run
	 	cmCar=glob.glob('%s*_Covariance_Matrix_*community_auto_results.out' % (treeSetTrunc))

	 	# Change name, might want to turn this into cp instead of mv
	 	os.system("mv %s %s_CovWholeCommunity_results.out" % (str(cmCar[0]), treeSetTrunc))


	 	print("Parse output into useful information")
	 	parse_output(clvPath, treeSet, treeSetTrunc, "Covariance", model, rooted, plateau)


	if network == 'Affinity':

		# Run Treescaper with manual plateau
		print("Running manual with lambda = %s. Log file: %s_AffCommunity.out" %  (plateau, treeSet))
		os.system("%s -trees -f %s -ft Trees -w 0 -r %s -o Community -t Affinity -cm %s -lm manu -dm URF -am Exp -lp %s -ln 0 " % (clvPath, treeSet, rooted, model, plateau)+\
		" > %s_AffCommunity.out" %  treeSet)

		# Run automatic plateau finder
		print("Running automatic. Log file: %s_AffAuto.out" %  treeSet)
		os.system("%s -trees -f %s -ft Trees -w 0 -r %s -o Community -t Affinity -cm %s -lm auto -dm URF -am Exp" % (clvPath, treeSet, rooted, model)+\
		" > %s_AffAuto.out" %  treeSet)

		# Get output file from automatic run
		aCar=glob.glob('%s*_Affinity-*community_auto_results.out' % (treeSetTrunc))

		# Change name, might want to turn this into cp instead of mv
		os.system("mv %s %s_AffWholeCommunity_results.out" % (str(aCar[0]), treeSetTrunc))

		print("Parse output into useful information")
		parse_output(clvPath, treeSet, treeSetTrunc, "Affinity", model, rooted, plateau)

 	os.system("sumtrees.py -r -o %s.con %s &> dendropy_%s.out" % (treeSetTrunc, inNexus,inNexus))
	os.system("cat ./SeqSim/FigTreeBlock.txt >> %s.con" % (treeSetTrunc))
 	#os.system("/Applications/FigTree/FigTree_v1.4.3/bin/figtree -graphic PDF all_trees.con all_trees.pdf")

	os.system("echo 'treescaperWrapperKnownPlateau.py %s\n' >> commands.txt" % ' '.join(sys.argv[1:]))


	#os.system("cat RAxML_bestTree.G1_gene* > RAxML_allTree.G1.tre")
if __name__=='__main__':
	main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Usage: treescaperWrapperV2.py [model] [interval] [network]

#***Check the CLVtreescaper settings within the script.  You might want to run CLVtreescaper with different options.  Make sure the
#numbering of all the indices is correct. The developers have changed between starting at 1 and starting at 0.  This includes
#the affinityCommunities script that is imported into this script.

#runs community detection on the tree set with range of lambda values specified. Finds the plateau of communities detected and prints to ouput files.
#the input tree set be a newick file named "all_trees.new" The lambda ranges and increments may need to be changed for affinity and covariance community
#detection. These sections are marked in the main() function.

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
import numpy as np
import fnmatch
from collections import Counter
import filecmp
from AffinityCommunities import affinityCommunityConsensus
from treescaperWrapper import mode_function, reg_ex_match
import dendropy
import StringIO
import heapq
from shutil import move, copyfile
from tempfile import mkstemp


def make_list(ls,convert):
#takes a string of elements seperated by tabs.  splits the string into a list.  allows user to select what format the elements are converted to.
	ls = ls.split("\t")
	ls.remove("\n")
	if convert == "int":
		ls = [int(i) for i in ls]
	if convert == "float":
		ls = [float(i) for i in ls]
	#remove first column, which does not provide useful values
	ls = ls[1:]
	return ls

def get_plateau(treeSet, treeSetTrunc, type, model, rooted):
#finds plateau of community detection results.  If a single plateau isnt found, returns None

	#open the community results file
	if type == "Covariance":
		comResults = open("%s_CovWholeCommunity_results.out" % treeSetTrunc)
	if type == "Affinity":
		comResults = open("%s_AffWholeCommunity_results.out" % treeSetTrunc)
	#make a list of the label of each unique community structure produced by each lambda value.  This number starts at 0 and increases by 1 every time a unique community structure is found.
	labelLS = comResults.readline()
	labelLS = make_list(labelLS, "int")
	print labelLS

	#make a list of the lambda values
	lambdaLS = comResults.readline()
	lambdaLS = make_list(lambdaLS,"float")
	print lambdaLS

	#make a list of community numbers produced by each lambda value
	comNumLS = comResults.readline()
	conNumLS = make_list(comNumLS, "int")
	print conNumLS
	m = max(labelLS)
	#assuming that the range of community structures went from one community for each node and one community for all nodes, remove the label values
	#that correspond to these community structures.
	labelLStrim= [i for i in labelLS if i != m and i != 0]
	print labelLStrim

	#also remove the second largest labels because those correspond to one community for every unique topology
	if type == "Affinity":

		m = max(labelLStrim)
		print m

		labelLStrim= [i for i in labelLStrim if i != m]
		print labelLStrim

	#if using traditional search for plateau, call mode_function
	if labelLStrim != []:
		plateauLabel = mode_function(labelLStrim)
	#if using automatic search and feeding the results, use the only result in the list
	else:
		plateauLabel = labelLS[0]
	print plateauLabel

	#if there is only one largest plateau, output plateau community structure.  if the communities are covariance communities, make a key for the bipartitions
	if len(plateauLabel) == 1:

		#find the lambda value that corresponds to the plateau community structure
		plateauIndex = labelLS.index(plateauLabel[0])
		plateauLambda = lambdaLS[plateauIndex]
		print plateauLambda


		if type == "Covariance":
			#save plateau community structure to file
			os.system("CLVTreeScaper -trees "+\
			"-f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm manu -lp 1 -ln %s -hf .95 -lf .05" % (treeSet, rooted, model, plateauLambda)+\
			" > %s_CovPlateauCommunity.out" %  treeSetTrunc)

			#get the number of communities
			comFile = open("%s_CovPlateauCommunity.out" %  treeSetTrunc, 'r')
			pattern = re.compile('Number of communities: (\d+)')
			coms = int(reg_ex_match(comFile, pattern))
			#make a key for the bipartitions
			comKey = open("%s_comKey.out" %  treeSetTrunc, 'w')

			for i in range(1, coms+1):
				#make a list of all the nodes in the community
				pattern = re.compile('Community '+str(i)+' includes nodes: (.+)')
				comStr = reg_ex_match(comFile, pattern)
				comLS = comStr.split(",")
				print comLS
				comLS = filter(None, comLS)
				comKey.write("Com %s:\n" % i)
				#find the bipartition that each nodes corresponds to
				for j in comLS:
					#j = int(j)+1
					pattern = re.compile("bipartition "+str(j)+" : ([0-1]+?), appear times: ([0-9]+?)$")
					#make a temperary copy of "*_CovPlateauCommunity.out", because it needs to be searched in two different ways
					fh, absPath = mkstemp()
					copyfile("%s_CovPlateauCommunity.out" %  treeSetTrunc, absPath)
					comTempFile = open(absPath,'r')
					#find the bipartition in the file
					for line in comTempFile:
						m = pattern.match(line)
						if m:
							#save the bipartition and its frequency
							bipart = m.group(1)
							freq = m.group(2)
							#make a list of all the taxa on one side of the biparition.  The position of each binary digit corresponds to the numbering
							# of the taxa in the output’s taxa list.  all the taxa on one side of an internal branch are assigned a “1”
							bipartLS = []
							for c in bipart:
								bipartLS.append(c)
							indices = [x+1 for x, y in enumerate(bipartLS) if y == '1']
							for k in indices:
								pattern2 = re.compile("(.+) , "+str(k))
								taxon = reg_ex_match(comFile, pattern2)
								indices[indices.index(k)] = taxon
					 #close temp file
					comTempFile.close()
					os.close(fh)

					comKey.write("%s %s\n" % (indices,freq))
				comKey.write("\n")
			comKey.close()










		if type == "Affinity":
			#save plateau community structure to file
			os.system("CLVTreeScaper -trees "+\
			"-f %s -w 0 -r %s -o Community -t Affinity -cm %s -lm manu -dm URF -am Exp -lp %s -ln 1 " % (treeSet, rooted, model, plateauLambda)+\
			" > %s_AffPlateauCommunity.out" %  treeSetTrunc)


		return plateauLambda



def main():
#convert newick formatted tree set to a nexus file in a format TreeScaper will read.

	#convert newick formatted tree set to a nexus file in a format TreeScaper will read
# 	tlst = dendropy.TreeList.get_from_path("all_trees.new", "newick", as_rooted=True)
# 	tlst.write_to_path("all_trees.pre.nex",'nexus')
# 	os.system("paup PaupFile")

	#save name of tree set without extension
	treeSet="all_trees.nex"
	treeSetIndex = treeSet.find(".")
	treeSetTrunc = treeSet[:treeSetIndex]
	model = sys.argv[1]
	interval = sys.argv[2]
	network = sys.argv[3]
	rooted = sys.argv[4]

	#find covariance community structures over a range of lambda values.

	#***These values may need to be adjusted to find the plateau****
	if network == 'Covariance':
		os.system("CLVTreeScaper -trees "+\
		"-f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm manu -lp 1 -lns .5 -lne 1.1 -lniv %s -hf .95 -lf .05" % (treeSet,rooted, model,interval)+\
		" > %s_CovCommunity.out" %  treeSet)

		#save the results to a file that won't be overwritten when TreeScaper is run again
		os.system("mv %s_Covariance\ Matrix_community_results.out %s_CovWholeCommunity_results.out" % (treeSetTrunc,treeSetTrunc))

		#search the community results for a single plateau. None is returned if one isnt found
		plateauLambda = get_plateau(treeSet, treeSetTrunc, "Covariance", model)

	#find affinity community structures over a range of lambda values

	#***These values may need to be adjusted to find the plateau****
	if network == 'Affinity':
		os.system("CLVTreeScaper -trees "+\
		"-f %s -ft Trees -w 0 -r %s -o Community -t Affinity -cm %s -lm manu -dm URF -am Exp -ln 1 -lps 0 -lpe 1 -lpiv %s" % (treeSet,rooted, model,interval)+\
		" > %s_AffCommunity.out" %  treeSet)
		os.system("mv %s_Affinity-URF_community_results.out %s_AffWholeCommunity_results.out" % (treeSetTrunc, treeSetTrunc))

		#search the community results for a single plateau. None is returned if one isnt found
		plateauLambda = get_plateau(treeSet, treeSetTrunc, "Affinity", model, rooted)

		#If there is a single plateau, find consensus tree of each affinity community
		if plateauLambda:
			affinityCommunityConsensus(treeSet,model,plateauLambda,rooted)

	#find consensus tree of tree set
 	os.system("sumtrees.py -r -o all_trees.con all_trees.nex")
	#make pdf of consensus tree
 	os.system("cat ~/Documents/BrownLabJash/Programs/bin/SeqSim/FigTreeBlock.txt >> all_trees.con")
 	os.system("figtree -graphic PDF all_trees.con all_trees.pdf")

	os.system("echo 'treescaperWrapperV2.py %s\n' >> commands.txt" % ' '.join(sys.argv[1:]))




	#os.system("cat RAxML_bestTree.G1_gene* > RAxML_allTree.G1.tre")
if __name__=='__main__':
	main()

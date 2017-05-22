#!/usr/bin/env python
##Make a nexus file of the trees in each affinity community. Make consensus tree of each affinity comunity


##Usage AffinityCommunities.py 'Path/To/Treeset' Model Plateau

##****The TreeScaper developers have changed what number the indices start from several times.  The may start counting at 0 or 1.  Make sure all the indices used
##in this script are right (for communities and trees).  Adjust the lines searching [&U] for unrooted trees or [&R] for rooted trees according to which type of trees you are using

#Output files

# AffinityCom[number].nex: a nexus file of the trees in an affinity community
# AffinityCom[number].nex.con: consensus tree of an affinity community
# AffinityCom[number].nex.con.pdf: pdf of consensus tree of an affinity community
# AffinityCommunitiesTreeCount.txt: a file to count the frequency and relative frequency of trees in each affinity community

from __future__ import division
import re
import os
import sys
import fnmatch
from shutil import move
from tempfile import mkstemp

def reg_ex_match(file, pattern):
#returns the first match of a reg ex search

	file.seek(0)
	for line in file:
		m = pattern.match(line)
		if m:
			return m.group(1)

def edit_treeset(treeFileEditPath):
#add comment blocks that number each tree with the indices used by TreeScaper

	print treeFileEditPath
	treeFileEdit = open(treeFileEditPath, 'r')


	lineNum = 1
	#make a temp file
	fh, absPath = mkstemp()
	tempFile = open(absPath,'w')
	for line in treeFileEdit:
		if line.find('[&U]') != -1: #Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
 			tempFile.write(line.replace('=','['+str(lineNum)+']='))
			lineNum += 1
		elif line.find('[&R]') != -1: #Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
 			tempFile.write(line.replace('=','['+str(lineNum)+']='))
			lineNum += 1
		else:
			tempFile.write(line)
	 #close temp file
	tempFile.close()
	os.close(fh)
	treeFileEdit.close()
	#Remove original file
	os.remove(treeFileEditPath)
	#Move new file
	move(absPath, treeFileEditPath)
	return lineNum

def affinityCommunityConsensus(treeFile,model,plateau,rooted):

	if model not in ('CPM','ERNM','CNM','NNM'):
		print "invalid model choose CPM, ERNM, CNM, or NNM"
		model = raw_input('Enter Model: ')
	print plateau
	#outputs community structure for current plateau value
	os.system("/Applications/MAC_TreeScaper_v1.0.0_Binary_2016-12-16/MAC_TreeScaper_v1.0.0_Binary_2016-12-16/CLVTreeScaper -trees "+\
	"-f %s -ft Trees -w 0 -r %s -o Community -t Affinity -dm URF -am Exp -cm %s -lm manu -lp %s -ln 0" % (treeFile, rooted, model, plateau)+\
	" > Affinity%s_%s_community.what" %  (model, plateau))

	print("AffinityCom.py")
	print(str(treeFile)+str(model)+str(plateau)+str(rooted))
	#getting number of communities from the output files
	comFile = open('Affinity%s_%s_community.what' % (model, plateau) , 'r' )
	pattern = re.compile('Number of communities: (\d+)')
	coms = int(reg_ex_match(comFile, pattern))
	print coms

	totalTrees = edit_treeset(treeFile)
	treeFile = open(treeFile,'r')

	#make a file to count the frequency and relative frequency of trees in each community
	treeCountFile = open("AffinityCommunitiesTreeCount.txt", 'w')


	#get the tree indices from each community
	for i in range(1, coms+1):
		pattern = re.compile('Community '+str(i)+' includes nodes: (.+)')
		comStr = reg_ex_match(comFile, pattern)
		comLs = comStr.split(",")
		print comLs
		comLs = filter(None, comLs)
		#Pull out these trees from the original trees. Make a nexus file containing the trees for each community
		comTreeSetStr = 'AffinityCom'+str(i)+'.nex'
		comTreeSet = open(comTreeSetStr,'w')
		treeCount = 0
		treeFile.seek(0)
		for line in treeFile:
			if line.find('[&U]') == -1 & line.find('[&R]') == -1:
				comTreeSet.write(line)
			else:
				for j in comLs:
					if line.find('['+str(j)+']') != -1:
						comTreeSet.write(line)
						treeCount += 1
		#Write frequency and relative frequency of trees in the community
		treeCountFile.write("%s\t%s\t%.2f%% of trees\n" % (comTreeSetStr, str(treeCount),100*(treeCount/totalTrees)))
		comTreeSet.close()
		#Make a consensus tree of the affinity community
		comTreeConStr = comTreeSetStr + ".con"
		if rooted == '0':
			os.system("sumtrees.py -r --unrooted -o %s %s" % (comTreeConStr,comTreeSetStr))
		if rooted == '1':
			print "hey"
			os.system("sumtrees.py -r --rooted -o %s %s" % (comTreeConStr,comTreeSetStr))
		os.system("cat /Users/ChatNoir/Projects/TreeScaper/treescaper_scripts_2017_test/SeqSim/FigTreeBlock.txt >> %s" % (comTreeConStr))
		#Make a pdf of the consensus tree
		os.system("figtree -graphic PDF %s %s.pdf" % (comTreeConStr, comTreeConStr))






if __name__=='__main__':
	affinityCommunityConsensus(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])


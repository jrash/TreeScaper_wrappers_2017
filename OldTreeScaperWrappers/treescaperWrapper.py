#!/usr/bin/env python
#Usage: treescaperWrapper.py 'path/to/configFile'

#Takes a config file in the format of "BrownLabJash/Programs/Scripts/treescaperWrapper/configExample.txt
#Moves into a folder for each rep of each branch length/sequence length pair. Runs affinity and covariance
#community detection on each set of bootstrap trees.

#runs community detection on the tree set with range of lambda values specified
#in the config file. If specified, uses three models for community detection: Configuration Null Model,
#Constant Potts Model, Erdos Renyi Null Model. Finds the plateau of communities detected and prints to ouput files

##Output Files

# ['type'] can be Covariance or Affinity
# ['model'] can be CNM,ERNM, or CPM
# ['plateauLambda'] the lambda value were the plateau was found

# ['type']_['model']_LambdaTest.txt: a file for recording the lambda values used and the number of communities found
# ['type']_['model']_community.cat: concatenated file of community detection results across all lambda values
# ['type']['model']_['plateauLambda']_community.txt: community structure of the plateau
# ['type']_modelCompareLambdaTest.txt: community detection results and stats for each community detection model
# ['type']_modelCompareLambdaTest.cat: all the community detection results concatenated and placed in sequence/branch length head directory
# AffinityCom[number].nex: a nexus file of the trees in each affinity community
# AffinityCom[number].nex.con: consensus tree of each affinity community
# AffinityCom[number].nex.con.pdf: pdf of consensus tree of each affinity community



import re
import os
import sys
import numpy as np
import fnmatch
from collections import Counter
import filecmp
from AffinityCommunities import affinityCommunityConsensus


def get_treeset():
#saves treeset to treeFile variable
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*nex.boot.tre'):
			treeFile=file
		elif fnmatch.fnmatch(file, '*run00.boot.tre'):
			treeFile=file
		elif fnmatch.fnmatch(file, '*burn.t'):
			treeFile=file
	return treeFile


def reg_ex_match(file, pattern):
#returns the first match of a reg ex search

	file.seek(0)
	for line in file:
		m = pattern.match(line)
		if m:
			return m.group(1)


def mode_function(lst):
#Returns a list of all the possible plateaus


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

def next_mode_function(lst):
#Returns a list of the next most frequent number(s) in a list

	counter = Counter(lst)
	if len(counter.most_common(2)) == 2:
		_,val = counter.most_common(2)[1]
		return [x for x,y in counter.items() if y == val]

def make_modelCompare_file(comList, modelCompare, treeFile, lambdaRange, nodesEqualsComs, oneCom, model):
#write community detection results and stats to file

	if comList != []:
# 		print mode_function(comList)
		modeStr = ""
		# find plateau(s), if there is more than one, save error flag for the plateau
		for i in mode_function(comList):
			count = comList.count(i)
			size = count*lambdaRange['Step']
			i = str(i)
			modeStr += i + ", "
		plateauComs = modeStr.split(",")
		plateauComs = filter(None, plateauComs)
		plateauComs.remove(' ')
		if len(plateauComs) > 1:
			modeStr = "Step Error"
			comList = []
		#if there is not one community for every bpartition and one community for all bipartitions in the range of communities - save error flag for the plateau
		if nodesEqualsComs == 0 or oneCom == 0:
			modeStr="Range Error" #moved this down from the top
			comList = []
		nextModeStr = ""
		#Range of lambda values with nontrivial communities
		lambdaListRange = len(comList)*lambdaRange['Step']
		#Variance of number of communities in the com list
		comListVar = np.var(comList)
		#find the second largest plateau
		if next_mode_function(comList):
			for i in next_mode_function(comList):
				count2 = comList.count(i)
				size2 = count2*lambdaRange['Step']
				i = str(i)
				nextModeStr += i + ", "
			modelCompare.write("%s\t%s\t%s\t%s\t%s\t%s\t%.3f\n" % (model,modeStr,size,nextModeStr,size2,lambdaListRange,comListVar))
		else:
			modelCompare.write("%s\t%s\t%s\tN/A\t\t%s\t%s\n"  % (model,modeStr,size,lambdaListRange,comListVar))

	else:
# 		print "empty"
		modelCompare.write("%s\tEmpty\n" % model)
	return comList


def com_detect(type, configFile, treeFile, LambdaPcov, seqLen, branchLen, rep):
#runs community detection on the tree set with range of lambda values specified
#in the config file. If specified, uses three models for community detection: Configuration Null Model,
#Constant Potts Model, Erdos Renyi Null Model. Finds the plateau of communities detected and prints to ouput file

	os.system("rm *.pdf")
	os.system("rm *.nex")
	os.system("rm *.nex.con")
    #os.system("rm *.txt")


	modelList = ['CNM','CPM','ERNM']
	rangeLabels = ['From','To','Step']
	#If there are any covariance lambda values in the config file, make a model compare file
	pattern = re.compile('^%s.*?Lambda Values %s:\s*(.*\d+)' % (type, model, rangeLabel))
	if reg_ex_match(configFile, pattern):
		#File for community detection results for each community detection model
		modelCompare = open("%s_modelCompareLambdaTest.txt" % type, 'w')
		modelCompare.write("\n%sbp%srep%s\n\nModel\tPlateau\tSize\t2ndPlateau\tSize\tLambda Value Range\tCommunity Variance\n" %(seqLen, branchLen, rep))

	for model in modelList:

		print "\n\n"+model+"\n\n"

		lambdaRange = dict()
		#if the range values for a model exists in the config file, save these values to a dictionary
		for rangeLabel in rangeLabels:

			pattern = re.compile('^%s %s Lambda Values %s:\s*(.*\d+)' % (type, model, rangeLabel))
			if reg_ex_match(configFile, pattern):
				lambdaRange[rangeLabel] = float(reg_ex_match(configFile, pattern))

		if 'From' in lambdaRange and 'To' in lambdaRange and 'Step' in lambdaRange:

			#make a file for recording the lambda values used and the number of communities found
			modelLambdaFile = open("%s_%s_LambdaTest.txt" % (type, model), 'w')


			diffCom = 0
			comList = []
			lambdaList = []
			coms = 0
			nodes = 1
			comsDiscard = 0
			nodesEqualsComs = 0
			oneCom = 0
			diffCom = 0
			allComList = []
			comList = []
			lambdaList = []
			#Makes lambda values go from highest to lowest, instead of lowest to highest for covariance communities.
			#Numbers of covariance communities starts at high numbers at low lambda values and goes to low numbers at high lambda values.
			#The inverse is true for affinity communities. This for simiplicity, so covariance and affinity community detection
			#can stop  when the number of communities = number of bipartitions.
			if type == 'Covariance':

				lambdaRange['Temp'] = lambdaRange['From']
				lambdaRange['From'] = lambdaRange['To']
				lambdaRange['To'] = lambdaRange['Temp']
				lambdaRange['Step'] *= -1
				print lambdaRange['Step']
				print lambdaRange
			for Lambda in np.arange(lambdaRange['From'],lambdaRange['To'],lambdaRange['Step']):

				if coms != nodes:
					lambdaList.append(Lambda)
					#outputs community structure for current lambda values
					if type == 'Covariance':
						os.system("CLVTreeScaper -trees "+\
						"-f %s -w 0 -r 0 -o Community -t Covariance -cm %s -lp %s -ln %s -hf .90 -lf .10" % (treeFile, model, LambdaPcov, Lambda)+\
						" > %s%s_%s_community.out" %  (type, model, Lambda))
					#outputs community structure for current lambda values ***rooted used because unrooted distance matrix is erroneous
					if type == 'Affinity':
						os.system("CLVTreeScaper -trees "+\
						"-f %s -w 1 -r 1 -o Community -t Affinity -dm URF -cm %s -lp %s -ln 1" % (treeFile, model, Lambda)+\
						" > %s%s_%s_community.out" %  (type, model, Lambda))

					#getting number of nodes and communities from the output files
					comFile = open('%s%s_%s_community.out' % (type, model, Lambda) , 'r' )
					pattern = re.compile('network has (\d+) nodes')
					nodes = int(reg_ex_match(comFile, pattern))
					pattern = re.compile('Number of communities: (\d+)')
					coms = int(reg_ex_match(comFile, pattern))
					allComList.append(coms)

					if Lambda == lambdaRange['From']:
						modelLambdaFile.write("Number of Nodes: %s\n\n" % nodes)
						modelLambdaFile.write("Lambda\tCommunities\n")

					#gets the community structure for the prior lambda value
					if Lambda != lambdaRange['From']:

						LambdaBefore = Lambda - lambdaRange['Step']
						print LambdaBefore
						comFileBefore = open('%s%s_%s_community.out' % (type, model, LambdaBefore), 'r' )
						pattern = re.compile('Number of communities: (\d+)')
						comsBefore = int(reg_ex_match(comFileBefore, pattern))

						if comsBefore != coms:
							#sets indicator of the same community number but different community structure back to 0 if the number of communities changes
							diffCom = 0

							#Save the number of communities right before the number of communities equals the number of nodes.  These community numbers will be discarded, because
							# the largest plateau will always be when there is a separate community for every unique topology
							if type == 'Affinity' and coms == nodes:

								comsDiscard = comsBefore
								print comsDiscard

							# If number of communities has reversed back to a prior community number, set lambdaBefore to the lambda value of the most recent community structure.
							# This is to check if if the community structure has reversed back to the most recent community structure for that community number.
							# Also, save the diffCom of the most recent community structure, so it can be correctly incremented by .001 if the community structures are different
							if coms in comList:

								comBeforeList = [x for x in comList if coms <= x < coms+1]
								diffCom = max(comBeforeList) - coms
								lastOcc = len(comList) - 1 - comList[::-1].index(max(comBeforeList))
								lambdaList[lastOcc] = LambdaBefore
								comFileBefore = open('%s%s_%s_community.out' % (type, model, LambdaBefore), 'r' )


					# makes a list of all the number of non trivial communities for each lambda value.
					if coms !=  nodes and coms != 1:

						#if community number has not occurred before, add to comlist --- edited this, had if lambda = lambda[start] append, if coms != comsbefore and coms not in comslist
						if coms not in comList:
							comList.append(coms)
						#if number of communities is the same as a community number before
						else:
							comMatchCount = 0
							# check to see if the community structure is the same the community structure before
							for i in range(1, coms+1):
								pattern = re.compile('Community '+str(i)+' includes nodes: (.+)')
								comsStr = reg_ex_match(comFile, pattern)
								for j in range(1, comsBefore+1):
									pattern = re.compile('Community '+str(j)+' includes nodes: (.+)')
									comsBeforeStr = reg_ex_match(comFileBefore, pattern)
									if comsStr == comsBeforeStr:
										comMatchCount += 1
							if comMatchCount == coms:
								coms += diffCom
								comList.append(coms)
							else:
								diffCom += .001
								coms += diffCom
								comList.append(coms)

					if coms ==  nodes:
						nodesEqualsComs += 1
					if coms == 1:
						oneCom += 1

					modelLambdaFile.write("%s\t%s\n" % (Lambda, coms))

			print comList
			comListTrim = []

			#Trim out the comumunity numbers that correspond to the number of unique topologies.  This will always be the largest plateau
			# if these community numbers arent trimmed. Also, write community detection results to file and save the plateau community number and lambda value
			if type == 'Affinity' and comList != []:
# 				try:
# 					comsDiscard
# 				except NameError:
# 					comsDiscard = 0
# 					print 'works'
				comList = [x for x in comList if x != comsDiscard]
				print comList
				if comsDiscard == 0:
					comList = []
					print comList
				if comList != []:

					plateauComNum = int(mode_function(comList)[0])
					plateauLambda = str(lambdaList[allComList.index(plateauComNum)])

			comList = make_modelCompare_file(comList, modelCompare, treeFile, lambdaRange, nodesEqualsComs, oneCom, model)

			#write community detection results to file and save the plateau community number and lambda value
			if type == 'Covariance' and comList != []:
				plateauComNum = int(mode_function(comList)[0])
				plateauLambda = str(lambdaList[allComList.index(plateauComNum)])
				#outputs community structure for current lambda values
				os.system("CLVTreeScaper -trees "+\
                "-f %s -w 0 -r 0 -o Community -t Covariance -cm %s -lp %s -ln %s -hf .90 -lf .10" % (treeFile, model, LambdaPcov, plateauLambda)+\
				" > %s%s_%s_community.txt" %  (type, model, Lambda))

			#Get the lambda value for the plateau and find the consensus trees for each community at this lambda value
			if type == "Affinity" and model == "CPM" and comList != []:
				if len(mode_function(comList)) == 1:
					affinityCommunityConsensus(treeFile,model,plateauLambda)

			#concatenate all community detection result files
			os.system("{ find . -name '*community.out' -print0 | xargs -0 cat; } > %s_%s_community.cat" % (type, model))
			#delete community detection result files after concatenation
			os.system("find . -name '*community.out' -print0 | xargs -0 rm")

			modelLambdaFile.close()

		#outputs community structure for current lambda values (NNM model)
		if type == 'Covariance':
			os.system("CLVTreeScaper -trees "+\
			"-f %s -w 0 -r 0 -o Community -t Covariance -cm NNM -hf .90 -lf .10" % treeFile+\
			" > %s_NNM_community.txt" %  type)
		#outputs community structure for current lambda values (NNM model)
		if type == 'Affinity':
			os.system("CLVTreeScaper -trees "+\
			"-f %s -w 1 -r 1 -o Community -t Affinity -dm URF -cm NNM" % treeFile+\
			" > %s_NNM_community.txt" %  type)

	#save NNM community detection results and write them to modelCompare file
	comFile = open('%s_NNM_community.txt' % type , 'r' )
	pattern = re.compile('Number of communities: (\d+)')
	coms = int(reg_ex_match(comFile, pattern))
	modelCompare.write("NNM\t%s\n" % coms)


	modelCompare.close()



def main():
#moves into a folder for each rep of each branch length/sequence length pair. Runs affinity and covariance community detection on each set of bootstrap trees.

	#move into community detection folder and open specified config file
	os.chdir("/home/vestige/Documents/BrownLabJash/RogueTaxaSimulation/ComDetection")
	configFilePath = str(sys.argv[1])
	configFile = open(configFilePath, 'r')

	# Save each set of Sequence lengths, branch lengths and rep numbers in the config file and cd into the specified folders
	pattern = re.compile("^Sequence Lengths:\s*(.*\d+)+")
	sequenceLen = reg_ex_match(configFile, pattern)
	sequenceLenList = []
	for i in sequenceLen.split():
		sequenceLenList.append(i)

	pattern = re.compile("^Branch Length Scales:\s*(.*\d+)+")
	branchLen = reg_ex_match(configFile, pattern)
	branchLenList = []
	for i in branchLen.split():

		branchLenList.append(i)

	pattern = re.compile("Reps:\s*(.*\d+)+")
	reps = reg_ex_match(configFile, pattern)
	repsList = []
	for i in reps.split():
		i = int(i)
		repsList.append(i)
	repsList[1] = repsList[1] + 1

	#save lambda+ value in config file
	pattern = re.compile("Covariance Lambda\+ Value:\s*(\d+)")
	LambdaPcov = reg_ex_match(configFile, pattern)

	#move into each sequence length, branch length, and rep folder specified, unless they dont exist, then skip
	for seqLen in sequenceLenList:
		for branchLen in branchLenList:
			for dir in os.listdir("."):
				#find folder labeled by seqLen and branchLen
				if re.search("%s[a-z]*%s" % (seqLen, branchLen), dir):
					os.chdir("./%s" % dir)
					for rep in range(repsList[0],repsList[1]):
						for dir in os.listdir("."):
							# find rep folder labeled with rep number
							if re.search("(?![0-9]+)[a-z]*%s(?![0-9]+)" % rep, dir):
								os.chdir(dir+"/MLboot")
								treeFile = get_treeset()

								#Save bipartition matrix to file
								os.system("CLVTreeScaper -trees "+\
								"-f %s -w 0 -r 0 -o BipartMatrix -bfm matrix > bipartMatrix.txt" % treeFile)
								#Save unweighted robinsons fold distance matrix to file
								os.system("CLVTreeScaper -trees "+\
								"-f %s -w 0 -r 1 -o Dist -dm URF" % treeFile)
								#Save weighted robinsons fold distance matrix to file
								os.system("CLVTreeScaper -trees "+\
								"-f %s -w 1 -r 1 -o Dist -dm RF" % treeFile)

 								#Covariance community detection
 								com_detect('Covariance', configFile, treeFile, LambdaPcov, seqLen, branchLen, rep)
 								#Affinity community detection
								com_detect('Affinity', configFile, treeFile, LambdaPcov, seqLen, branchLen, rep)

								print os.getcwd()
								os.chdir("../..")

					#concatenate all the community detection stats and place file in sequence/branch length head directory
					os.system("cat rep*/MLboot/Affinity_modelCompareLambdaTest.txt > Affinity_modelCompareLambdaTest.cat")
					os.system("cat rep*/MLboot/Covariance_modelCompareLambdaTest.txt > Covariance_modelCompareLambdaTest.cat")
					os.chdir("..")






if __name__=='__main__':
	main()

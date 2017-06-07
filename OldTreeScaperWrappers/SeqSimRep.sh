#!/bin/bash
#Usage SeqSimRep.sh 'sequence lengths' 'scale of branch lengths' 'rep number from' 'rep number to'
# ./SeqSimRep.sh '1000 5000' '.1 1' 1 5
#Produces 4 folders - sequence length 1000 with .1 and 1 branch length scales and 5000 with .1 and 1 branch length scales
#Each folder contains replicates 1-5

#makes a folder for each sequence length and branch length pair specified.  makes a folder for each rep.  Simulates
#sequences with conflicting signal and then does a mr bayes and garli bootstrap analysis.

#input

#Tree[1/2].phy: the guide trees used to simulate sequence alignments


#output


#SimSeqs/Tree[1/2].seqs.nex: the sequence alignments simulated from either guide tree
#SimSeqs/Tree.cat.seqs.nex: Both sequence alignments concatenated together

#Bayes/Tree.cat.seqs.nex.burn.t: the set of posterior trees from the 5 runs, with burnin removed
#Bayes/Tree.cat.seqs.nex.con.tre: the consensus tree of the above tree set
#Bayes/Tree.con.Garli.${1}bp${2}L.pdf: pdf of the consensus tree

#MLboot/Tree.cat.seqs.nex.boot.tre:  The bootstrap tree set from the garli bootstrap analysis
#MLboot/Tree.cat.seqs.nex.con: the consensus tree of the above tree set
#MLboot/Tree.con.Garli.${1}bp${2}L.pdf: pdf of the consensus tree

#MLbest/Tree.cat.seqs.nex.best.tre: garli maximum likelihood tree

SeqGen () { 
#Uses two guide trees to simulate sequences in Seq-Gen, and concatenates the sequences for each taxon together.
#Analyzes the sequence alignment in Garli (ML tree and bootstrap analysis) and Mr.Bayes.  Makes pdf's of consensus trees

	#simulating two equally sized MSAs using two different guide trees
 	halfSeq=$(($1/2))
 	seq-gen -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l${halfSeq} -s$2 -on  < Tree1.phy > Tree1.seqs.nex
 	seq-gen -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l${halfSeq} -s$2 -on  < Tree2.phy > Tree2.seqs.nex
	#concatenating the MSAs together
 	combineGenes.py Tree1.seqs.nex Tree2.seqs.nex Tree.cat.seqs.prenex 
	#removing the partition block
 	cat Tree.cat.seqs.prenex | head -n 23 > Tree.cat.seqs.nex
 	rm Tree.cat.seqs.prenex
	#adding the Mr Bayes block to the MSA, check this file for specs on how mr bayes was run 
 	cat ~/Documents/BrownLabJash/Programs/Scripts/SeqSim/MBblock.txt >> Tree.cat.seqs.nex
 	mb Tree.cat.seqs.nex
	#move mr bayes results to "Bayes"
 	mv Tree.cat.seqs.nex.* ../Bayes
 	cd ../Bayes
	#run a sequence of Paup commands that fill remove the 25% burnin from all five Mr.Bayes runs, and then concatenate all five runs together 
	#into one large treset.  Check file for specifics.  
 	paup ~/Documents/BrownLabJash/Programs/Scripts/SeqSim/PaupBurn.txt
	#Add a file to the MSA that will specify how the consensus tree will be visualized in FigTree
 	cat ~/Documents/BrownLabJash/Programs/Scripts/SeqSim/FigTreeBlock.txt >> Tree.cat.seqs.nex.con.tre
	#create a pdf of FigTree visualization 
 	figtree -graphic PDF Tree.cat.seqs.nex.con.tre Tree.con.Garli.${1}bp${2}L.pdf
 	cd ../SimSeqs
	#do 100 boostrap replicate analysis in garli.  see garli config file "RTaxaTree.JC69.boot" for details
 	Garli-2.0 ~/Documents/BrownLabJash/Programs/Garli-2.0-IntelOSX-multithread/bin/RTaxaTree.JC69.boot
	#move boostrap results to "MLboot"
 	mv Tree.cat.seqs.nex.* ../MLboot
	cd ../MLboot
	#make a consensus tree of the 100 boostrap replicates
	sumTrees.py -r -o Tree.cat.seqs.nex.con Tree.cat.seqs.nex.boot.tre
	#Add a file to the MSA that will specify how the consensus tree will be visualized in FigTree
	cat ~/Documents/BrownLabJash/Programs/Scripts/SeqSim/FigTreeBlock.txt >> Tree.cat.seqs.nex.con
	#create a pdf of FigTree visualization 
	figtree -graphic PDF Tree.cat.seqs.nex.con Tree.con.Garli.${1}bp${2}L.pdf
 	cd ../SimSeqs
	#do ML tree analysis in garli.  see garli config file "RTaxaTree.JC69" for details
 	Garli-2.0 ~/Documents/BrownLabJash/Programs/Garli-2.0-IntelOSX-multithread/bin/RTaxaTree.JC69
	#move ML tree results to "MLbest"
 	mv Tree.cat.seqs.nex.* ../MLbest

	
}

cd ~/Documents/BrownLabJash/RogueTaxaSimulation/ComDetection

for length in $1 #Enter list of desired sequence lengths
do
    for scale in $2  # Enter list of desired scale of branch lengths
    do 
    	for i in `seq $3 $4`
    	do
      		mkdir ${length}bp${scale}L
    		cd ${length}bp${scale}L
 			cp -r ../TemplateFolder rep$i
        	echo "$length$scale$i"
    		cd rep$i/SimSeqs
    		SeqGen $length $scale
    		cd ../../..
    	done
    done  
done


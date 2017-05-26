#!/bin/bash
rm *.what
rm *.out
rm AffinityCom*.con
rm AffinityCom*.nex
rm AffinityCommunitiesTreeCount.txt
rm all_trees.con

# Try sample commands
#./CLVTreeScaper -trees -f all_trees.nex -ft Trees -w 0 -r 0 -o Community -t Covariance -cm CPM -lm auto -hf .95 -lf .05 > all_trees_CovAuto.out

#./CLVTreeScaper -trees -f all_trees.nex  -ft Trees -w 0 -r 0 -o BipartMatrix -bfm matrix 

#./CLVTreeScaper -trees -f all_trees.nex -ft Trees -w 0 -r 0 -o Community -t Affinity -cm CPM -lm auto -dm URF -am Exp > all_trees_AffAuto.out

#!/bin/bash
# Run this from one folder above your git folder

cd ./TreeScaper_wrappers_2017

new_branch='ggm_May28'


git branch $new_branch
git checkout $new_branch
git push -u origin $new_branch


git checkout master
git pull
git checkout $new_branch
git rebase master

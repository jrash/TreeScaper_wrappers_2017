#!/bin/bash

cd ./TreeScaper_wrappers_2017

new_branch='May31_ggm'


git branch $new_branch
git checkout $new_branch
git push -u origin $new_branch


git checkout master
git pull
git checkout $new_branch
git rebase master

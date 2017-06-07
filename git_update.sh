#!/bin/bash


message="Pull out updated plateau from out file"


echo $message

git add *
git commit -a -m "$message"
git push

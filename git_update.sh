#!/bin/bash


message="Cleaning up folder"


echo $message

git add *
git commit -a -m "$message"
git push

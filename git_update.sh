#!/bin/bash


message="Fussing with output file"


echo $message

git add *
git commit -a -m "$message"
git push

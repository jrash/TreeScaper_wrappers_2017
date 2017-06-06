#!/bin/bash

message="Scripts to set up runs on the cluster"

echo $message

git add *
git commit -a -m "$message"
git push

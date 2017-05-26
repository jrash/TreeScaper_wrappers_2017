#!/bin/bash

message="Fixing CLV path issues"

echo $message

git add *
git commit -a -m "$message"
git push

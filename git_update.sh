#!/bin/bash

message="Fixing issues with lambda being set to 0"

echo $message

git add *
git commit -a -m "$message"
git push

#!/bin/bash

message="Avoid adding repetitive numbers to original .nex file"

echo $message

git add *
git commit -a -m "$message"
git push

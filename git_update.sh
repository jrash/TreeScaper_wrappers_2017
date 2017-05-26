#!/bin/bash

message="add job script"

echo $message

git add *
git commit -a -m "$message"
git push

#!/bin/bash

message="removing output files"

echo $message

git add *
git commit -a -m "$message"
git push

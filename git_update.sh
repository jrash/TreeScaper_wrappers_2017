#!/bin/bash

message="small edits"

echo $message

git add *
git commit -a -m "$message"
git push

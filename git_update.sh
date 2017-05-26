#!/bin/bash

message="edits"

echo $message

git add *
git commit -a -m "$message"
git push

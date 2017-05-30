#!/bin/bash

message="fix new branch script"

echo $message

git add *
git commit -a -m "$message"
git push

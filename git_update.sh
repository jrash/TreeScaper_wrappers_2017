#!/bin/bash

message="delete output file"

echo $message

git add *
git commit -a -m "$message"
git push

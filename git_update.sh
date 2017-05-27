#!/bin/bash

message="fixing input file issues"

echo $message

git add *
git commit -a -m "$message"
git push

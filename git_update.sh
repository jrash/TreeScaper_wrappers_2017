#!/bin/bash

message="Adding figtree CLV and using in scripts"

echo $message

git add *
git commit -a -m "$message"
git push

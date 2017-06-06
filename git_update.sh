#!/bin/bash

message="Output formatting"

echo $message

git add *
git commit -a -m "$message"
git push

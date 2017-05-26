#!/bin/bash

message="deleting extra files"

echo $message

git add *
git commit -a -m "$message"
git push

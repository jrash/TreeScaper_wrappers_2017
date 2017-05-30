#!/bin/bash

message="Adding timekeeping to wrappers"

echo $message

git add *
git commit -a -m "$message"
git push

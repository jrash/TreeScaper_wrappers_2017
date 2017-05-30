#!/bin/bash

message="Commenting and cleaning up wrappers"

echo $message

git add *
git commit -a -m "$message"
git push

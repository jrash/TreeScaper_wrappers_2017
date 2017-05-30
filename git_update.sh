#!/bin/bash

message="Commenting scripts and removing seemingly unnessecary parts"

echo $message

git add *
git commit -a -m "$message"
git push

#!/bin/bash

message="Editing and commenting Affinity.py"

echo $message

git add *
git commit -m "$message"
git push

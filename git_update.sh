#!/bin/bash

message="fix git script"

echo $message

git add *
git commit -a -m "$message"
git push

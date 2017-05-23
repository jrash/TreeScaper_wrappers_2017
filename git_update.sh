#!/bin/bash

message="Works locally. Commit before cluster test"

echo $message

git add *
git commit -m "$message"
git push

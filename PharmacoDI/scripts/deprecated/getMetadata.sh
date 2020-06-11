#!/bin/bash

echo "Cloning the annotation repository..."
git clone https://github.com/BHKLAB-Pachyderm/Annotations.git
echo "Done."

echo "Copying annotation files to metadata directory..."
cp Annotations/* .
echo "Done."
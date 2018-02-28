#!/bin/bash
find . -type f -name "*.gcov" -delete
for filename in $(find ./src -type f -name "*.cxx"); do
    (cd $(dirname $filename) && \
         gcov --branch-probabilities \
              --branch-counts \
              --preserve-paths \
              --relative-only \
              -o . $filename)
done

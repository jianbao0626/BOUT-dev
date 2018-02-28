#!/bin/bash
find . -type f -name "*.gcov" -delete
for filename in $(find . -type f -name "*.gcno"); do
    (cd $(dirname $filename) && \
         gcov --branch-probabilities \
              --branch-counts \
              --preserve-paths \
              --relative-only \
              -o . $filename)
done

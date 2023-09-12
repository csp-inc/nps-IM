#!/usr/bin/env bash

find ./output -mindepth 5 -maxdepth 5 -type d \
    -exec ./src/compile-results.r -d {} \;
# find output -type f -name "mod-summaries.csv" -delete
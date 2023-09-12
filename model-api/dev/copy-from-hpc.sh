#!/usr/bin/env bash

# NOTE: on monsoon first run: find output -type f | grep -i mod-summary.csv | xargs cp {} --parents -t output-tmp
# https://superuser.com/questions/312348/linux-copy-all-files-matching-pattern-from-dir-and-subdirs-into-a-single-dir/312394
scp -r lz62@monsoon.hpc.nau.edu:/scratch/lz62/uplands-ci/output-tmp/output output-tmp

#scp -r lz62@monsoon.hpc.nau.edu:/scratch/lz62/uplands-ci/logs output/
# scp -r lz62@monsoon.hpc.nau.edu:/scratch/lz62/uplands-ci/output/every-model-summary.csv output/
# scp -r lz62@monsoon.hpc.nau.edu:/scratch/lz62/uplands-ci-may/output/every-model-summary.csv output/

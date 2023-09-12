#!/bin/bash

LOGS=output/logs/*.log
mkdir -p output/logs/out-of-memory

for f in $LOGS; do
  if grep -q "Out Of Memory" "$f"; then
    # Some Actions # SomeString was found
    echo "Processing $f file..."
    # mv "$f" output/logs/out-of-memory/
  fi
  # take action on each file. $f store current file name
  # cat $f
done


# for filename in /output/logs/*.log; do
#   echo $filename
#     # for ((i=0; i<=3; i++)); do
#     #     ./MyProgram.exe "$filename" "Logs/$(basename "$filename" .txt)_Log$i.txt"
#     # done
# done


# for file in $(grep -l "*.log" $1/*); do
#     mv $file $2;
# done

#!/usr/bin/env bash

DIR="${1}"

rm -f inits-needed.txt
cd "${DIR}"
find -exec sh -c "cat *.jags 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null 2>/dev/null | sed -n '/~/p' | sed 's/~.*//' | sed 's/\[.*\]//' | sed 's/^\s*//' | sort -u >> tmp-inits-needed.txt" \;
cat tmp-inits-needed.txt | sort -u > inits-needed.txt
rm -f tmp-inits-needed.txt
cd -

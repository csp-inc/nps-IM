#!/usr/bin/env bash

for ISS in $(ls in/ | sed 's/^iss-//' | sed 's/-.*$//' | sort -u); do
    for i in in/iss-${ISS}-*; do
        for j in out/iss-${ISS}-*;do
            echo $i $j;
            diff -rq $i $j;
        done
    done
done

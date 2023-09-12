#!/usr/bin/env bash

WRITE_DIR=$(dirname "${1}")
cat "${1}" | sed -n '/~/p' | sed 's/~.*//' | sed 's/\[.*\]//' | sed 's/^\s*//' | sed 's/[ \t]*$//' | sort -u

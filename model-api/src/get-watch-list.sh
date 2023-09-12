#!/usr/bin/env bash

WRITE_DIR=$(dirname "${1}")
cat "${1}" | sed -n '/<-\|~/p' | sed 's/<-.*//' | sed 's/~.*//' | sed 's/\[.*\]//' | sed 's/^\s*//' | sort -u > "${WRITE_DIR}/00-input/vars-to-watch.txt"

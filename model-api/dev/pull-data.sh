#!/usr/bin/env bash

AZURE="${1}"

mkdir -p data
gsutil -m rsync -r -c gs://nps-ima/data data

mkdir -p config
gsutil -m rsync -r -c gs://nps-ima/config config

if [ "${AZURE}" == 'azure' ]; then  # e.g., ./pull-data.sh azure
    mkdir -p review review/in review/out review/original
    gsutil -m rsync -r -c gs://nps-ima/review-in review/in
    gsutil -m rsync -r -c gs://nps-ima/review-out review/out
    gsutil -m rsync -r -c gs://nps-ima/review-original review/original
fi

#!/usr/bin/env bash

NETWORK=$1
UNIT=$2
FILE=$3

network_park_suffix="${NETWORK,,}-${UNIT,,}"

filelist_suffix="${network_park_suffix}"
if [ "${FILE}" != "" ]; then
	FILE_SANS_EXT="${FILE%%.*}"
  filelist_suffix="${network_park_suffix}-${FILE_SANS_EXT}"
fi

echo "output/hpc/filelist-${filelist_suffix}"

./dev/debugging/batch-analysis.sh "${NETWORK}/${UNIT}/${FILE}" deploy TRUE FALSE run

ls -R -d $PWD/output/hpc/workspace-${filelist_suffix}*/*.RData > output/hpc/filelist-${filelist_suffix}

N_JOBS=$(cat output/hpc/filelist-${filelist_suffix} | wc -l)
sed "s/LINE_COUNT/${N_JOBS}/g" dev/ap-job-array-template.sh > dev/submit-${filelist_suffix}.sh
sed -i "s/FILENAME/${filelist_suffix}/g" dev/submit-${filelist_suffix}.sh
chmod +x dev/submit-${filelist_suffix}.sh
mkdir -p "/scratch/lz62/uplands-ci/logs/${filelist_suffix}"

echo "You will be running ${N_JOBS} jobs...."
echo "Run with: dev/submit-${filelist_suffix}.sh"

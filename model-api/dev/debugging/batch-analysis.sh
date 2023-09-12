#!/usr/bin/env bash

SUBDIR="${1}"  # either 'all' or <network code>, <network code>/<unit code>, etc.
TEST="${2}"  # either 'test' or 'deploy'
MAKE_FILELIST="${3}"  # either TRUE or FALSE
OVERRIDE_DQS="${4}"  # either TRUE or FALSE
RUN="${5}"  # either omit just to list analyses, or 'run' to compile analyses
NUM_CPUS="${6}"

SEARCH_PATH="./config"
if [ "${SUBDIR}" != 'all' ]; then
    SEARCH_PATH="./config/${SUBDIR}"
fi

CPUS=1
if [ "${NUM_CPUS}" != '' ]; then
    CPUS="${NUM_CPUS}"
fi

if [ "${TEST}" == 'test' ]; then
    A_FLAGS="-a 300 -u 500 -n 300 -t TRUE -c ${CPUS} -m ${MAKE_FILELIST} -d ${OVERRIDE_DQS} -l logs/killed_log_paths"
else
    #A_FLAGS="-a 1500 -u 15000 -n 3000 -c ${CPUS} -m ${MAKE_FILELIST}"
    A_FLAGS="-a 5000 -u 15000 -n 5000 -c ${CPUS} -m ${MAKE_FILELIST} -d ${OVERRIDE_DQS} -l logs/killed_log_paths"
fi

# Note: Filenames cannot contain newlines...
mapfile -t my_config < <(find "${SEARCH_PATH}" -name *.yml -not -path "${SEARCH_PATH}/MISC/*" \
    ! -name ".*" \
    ! -name _network-level-attributes.yml \
    ! -name _park-level-attributes.yml)

if [ "${RUN}" == 'run' ]; then
    for c_file in "${my_config[@]}";do
        myscript=$(echo Rscript analysis-pipeline.r -f "${c_file}" "${A_FLAGS}")
        echo "${myscript}"
        $($myscript)
    done
else
    for c_file in "${my_config[@]}";do
        echo Rscript analysis-pipeline.r -f "${c_file}" "${A_FLAGS}"
    done
fi

#!/usr/bin/env bash

RELEASE="${1}"  # or use $(git describe)

rm -f uplands-ci.zip Rplots.pdf
gsutil -m rsync -c -x ".*\.DS_Store$|MISC/|.*-archive/" -r data gs://nps-ima/data

DASH_SUFFIX=""
if [ "${RELEASE}" != '' ]; then  # e.g., ./push-data.sh v0.3.1-ci.0
    zip -r "data-${RELEASE}.zip" data/* -x data/MISC/\* *.DS_Store
    gsutil cp "data-${RELEASE}.zip" gs://nps-ima
    rm "data-${RELEASE}.zip"
    DASH_SUFFIX="-${RELEASE}"
fi

## declare an array variable
declare -a arr=("data" "config")

## now loop through the above array
apt-get update
apt-get install -y zip unzip -qy
for i in "${arr[@]}"
do
   echo "$i"
   gsutil -m rsync -c -x ".DS_Store|MISC/|.*-archive/" -r "$i" "gs://nps-ima/$i"
   zip -r "$i.zip" "$i"/* -x $i/MISC/\* *.DS_Store
   gsutil cp "$i.zip" gs://nps-ima
   rm "$i.zip"
done

zip -r "uplands-ci${DASH_SUFFIX}.zip" . -x \*MISC/\* output/\* *-archive/\* .rstudio/\* rstudio/\* logs/\* sandbox/\* *.DS_Store *.zip Rplots.pdf
gsutil cp "uplands-ci${DASH_SUFFIX}.zip" gs://nps-ima
rm "uplands-ci${DASH_SUFFIX}.zip"

gsutil rm gs://nps-ima/**.DS_Store

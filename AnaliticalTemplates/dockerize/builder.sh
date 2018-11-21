#!/usr/bin/env bash

SOURCE="$1"
EXITPATH="$2"

FILE=${SOURCE##*/}
BASE=${FILE%.*}

echo -e "\e[31mReaching tmp dir...\e[0m"
cd /tmp &&
echo -e "\e[31mCopying package into the container...\e[0m"
cp ${SOURCE} . &&
echo -e "\e[31mUnzipping package...\e[0m"
unzip ${FILE} &&
echo -e "\e[31mAccessing compilation files...\e[0m"
cd ${BASE} &&
echo -e "\e[31mCompiling...\e[0m"
make all &&
echo -e "\e[31mCompressing result files\e[0m"
cd .. && zip -r ${FILE} ${BASE} &&
echo -e "\e[31mMoving final package...\e[0m"
cp ${FILE} ${EXITPATH} &&
echo -e "\e[31mCleaning...\e[0m"
rm -rf /tmp/*


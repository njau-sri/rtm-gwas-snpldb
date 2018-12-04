#!/bin/bash

rm -rf $1
mkdir $1

TARGET=$1/rtm-gwas-snpldb

if [ -z "$RTM_GWAS_VERSION" ]; then
    RTM_GWAS_VERSION=unknown
fi

if [ $1 == "glnx64" ]; then

    g++ *.cpp -o $TARGET -s -O2 -std=c++11 -static -fopenmp \
        -DRTM_GWAS_VERSION=$RTM_GWAS_VERSION

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o $TARGET.exe -s -O2 -std=c++11 -static -fopenmp \
        -DRTM_GWAS_VERSION=$RTM_GWAS_VERSION

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o $TARGET.exe -s -O2 -std=c++11 -static -fopenmp \
        -DRTM_GWAS_VERSION=$RTM_GWAS_VERSION

fi

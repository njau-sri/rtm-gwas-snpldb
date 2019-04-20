#!/bin/bash

rm -rf macos
mkdir macos

TARGET=macos/rtm-gwas-snpldb

if [ -z "$RTM_GWAS_VERSION" ]; then
    RTM_GWAS_VERSION=unknown
fi

g++ src/*.cpp -o $TARGET -O2 -std=c++11 \
    -DRTM_GWAS_VERSION=$RTM_GWAS_VERSION

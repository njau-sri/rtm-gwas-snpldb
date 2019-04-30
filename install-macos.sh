#!/bin/bash

export PATH="/usr/local/opt/llvm@7/bin:$PATH"
export LDFLAGS="-L/usr/local/opt/llvm@7/lib"
export CPPFLAGS="-I/usr/local/opt/llvm@7/include"

rm -rf macos
mkdir macos

TARGET=macos/rtm-gwas-snpldb

if [ -z "$RTM_GWAS_VERSION" ]; then
    RTM_GWAS_VERSION=unknown
fi

clang++ src/*.cpp -o $TARGET -O2 -std=c++11 -fopenmp \
    -DRTM_GWAS_VERSION=$RTM_GWAS_VERSION

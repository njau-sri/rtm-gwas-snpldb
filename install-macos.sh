#!/bin/bash

rm -rf macos
mkdir macos

LIB=/usr/local/lib
LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

g++ *.cpp -o macos/rtm-gwas-snpldb -O2 -std=c++11 $LIBGCC

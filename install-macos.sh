#!/bin/bash

rm -rf macos
mkdir macos

g++ *.cpp -o macos/rtm-gwas-snpldb -O2 -std=c++11

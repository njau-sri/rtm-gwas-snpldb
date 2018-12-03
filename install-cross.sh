#!/bin/bash

rm -rf $1
mkdir $1

if [ $1 == "glnx64" ]; then

    g++ *.cpp -o $1/rtm-gwas-snpldb -s -O2 -std=c++11 -static -fopenmp

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o $1/rtm-gwas-snpldb.exe -s -O2 -std=c++11 -static -fopenmp

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o $1/rtm-gwas-snpldb.exe -s -O2 -std=c++11 -static -fopenmp

fi
#!/bin/bash

# Generates packaging for packagemonkey

rm -f Makefile rpmpackage/libgpr.spec

packagemonkey -n "libgpr" --version "1.03" --cmd --dir "." -l "bsd" -e "Bob Mottram (4096 bits) <bob@robotics.uk.to>" --brief "Library for genetic programming" --desc "Making the inclusion of Genetic Programming easy within any C/C++ application. Genetic programming (GP) is a powerful technique, inspired by the process of natural selection, which can be utilized to automatically discover programs which produce a desired input to output transformation. Both classical tree-based and Cartesian forms of Genetic Programming are supported, including self-modifying variants." --homepage "https://github.com/fuzzgun/libgpr" --repository "https://github.com/fuzzgun/libgpr.git" --section "libs" --categories "Development/ArtificialIntelligence" --cstandard "c99" --compile "-lm -lz -fopenmp" --dependsdeb "gnuplot, libz-dev" --dependsarch "gnuplot, libzip"

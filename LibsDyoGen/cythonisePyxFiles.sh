#!/bin/bash
# cythonise the utils/extractDiags.pyx file into utils/extractSbs.c
# the first argument should be the path to the root of LibsDyogen

PATH_LIBSDYOGEN=${1}
cython ${PATH_LIBSDYOGEN}/utils/extractSbs.pyx
# C compilation of utils/extractSbs.c into utils/extractSbs.so that can be used as a python package
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o ${PATH_LIBSDYOGEN}/utils/extractSbs.so ${PATH_LIBSDYOGEN}/utils/extractSbs.c

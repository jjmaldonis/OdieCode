#!/bin/bash
IFS=:

for p in ${LD_LIBRARY_PATH}; do
    if [ -e ${p}/libimf.so ]; then
    #if [ -e ${p}/mpif.h ]; then
        echo ${p}
    fi
done

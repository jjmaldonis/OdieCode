#!/bin/bash

for f in *.sh
do
    sed -i '5i #$ -V' $f
    sed -i -e 's/\r//' $f
done

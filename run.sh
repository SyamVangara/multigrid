#!/bin/bash

name=$(sed -n "6 p" poisson.in)
#echo -n "Write data to folder > "
#read name
./poisson
mkdir -p data/$name
mv -v *.out data/$name/
cp -r *.in data/$name/
mv -v *.err data/$name/
mv -v *.dat data/$name/



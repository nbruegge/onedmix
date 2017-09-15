#!/bin/bash

#name='fib1'
name='diffusion'

rm ${name}.pyf
rm ${name}.so

echo "Build step 1:"
f2py -m ${name} -h ${name}.pyf test_mod.f90 ${name}.f90
if [ $? != 0 ]; then echo "::: Error! :::"; exit 2; fi

echo "Build step 2:"
f2py -c ${name}.pyf test_mod.f90 ${name}.f90
if [ $? != 0 ]; then echo "::: Error! :::"; exit 2; fi

echo "All done!"

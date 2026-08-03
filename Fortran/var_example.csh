#!/bin/csh

# To compile this code you need to:
# module load OneAPI
# source $SETVARS

ifx var_example.f90 -o var_example.x
./var_example.x

# Burrows-Wheeler transformation and FM index

This repository contains educational implementations of Burrows-Wheeler Tranformation and Ferragina-Manzini index.

This project contains four modules.

- bwt_fm.py - This is main module. Here we read the input file and process it so that all necessary data is initialized
- bwt_build.py - This module is used for creation of Burrows-Wheeler transformation, and suffix array.
- fm_build.py - This module is used for creation of FM index.
- utils.py - This module contains auxilary functions. 

## How to run

Run this project by running *bwt_fm.py* script.

### Test file

You can either use existing default file, or pass your own.
If you choose to run your own file, you have to pass it's path as argument.
Test file should be in conventient format, first line should have some information about processed file, and starting from the second line there should be a string for processing.

### Patterns

There are already three default patterns in this project for default file.
If you choose to run your own file and search for different patterns, you can pass them as arguments in command line after running the *bwt_fm.py* script.

### Factors

After inputing patterns, you will be asked to insert suffix array and tally matrix factor.
Theese steps can be skipped, in which case default values (1,1) will be passed.
In case you want to run optimized solution you should pass two integer values.

## Video presentation

Video presentation in available on YouTube: https://www.youtube.com/watch?v=YYzolQjCGtw

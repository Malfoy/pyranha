#!/usr/bin/env python
import subprocess,sys,struct,os,time,shutil,random

sequence = ""
with open(sys.argv[1]) as infile:  # do not load the whole file, one line at time
    for line in infile:
        string = line.rstrip()
        if ">" not in string:
            sequence += string
        else:
            if sequence != "":
                print sequence
            print string
            sequence = ""
    print sequence

#!/usr/bin/env python

import sys
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-j", help = "Number of jobs to be execute in a bash file", dest = "j", default = 1, type = int)
parser.add_argument("-q", help = "Which nodes to be use", dest = "q", default = "ngs4G")
parser.add_argument("-p", help = "Project ID", dest = "p", default = "MST107119")
parser.add_argument("-n", help = "Job names", dest = "n", default = "job")
parser.add_argument("-c", help = "Number of cpu been used", dest = "c", default = 5, type = int)
args = parser.parse_args()

def print_pbs(output):
    f = open(output + '.sh', 'x')
    f.write("#!/bin/bash\n")
    f.write("#PBS -P %s\n" % args.p)
    f.write("#PBS -W group_list=%s\n" % args.p)
    f.write("#PBS -N %s\n" % output)
    #f.write("#PBS -l select=1:ncpus=%s\n" % args.c)
    #f.write("#PBS -l place=pack\n")
    f.write("#PBS -q %s\n" % args.q)
    f.write("#PBS -o %s\n" % ''.join(["./" + output + ".out"]))
    f.write("#PBS -e %s\n" % ''.join(["./" + output + ".err"]))
    f.write("cd $PBS_O_WORKDIR\n\n")

count = -1
no_output = 1

for line in sys.stdin: 
    if line == "### job\n": 
        count += 1

    if count == -1 and line != "\n":
        file_name = args.n + str(no_output) 
        print_pbs(file_name)
        count += 1
        no_output += 1
    
    if count < args.j:   
        print(line, file = open(file_name + ".sh", 'a'), end = "")
    elif count == args.j:
        count = -1
        print(line, file = open(file_name + ".sh", 'a'))
        os.system("chmod +x %s" % file_name + ".sh")
        os.system(''.join(["qsub ./", file_name, ".sh"]))

if count < args.j and count != -1:
    os.system("chmod +x %s" % file_name + "*.sh")
    os.system(''.join(["qsub ./", file_name, ".sh"]))

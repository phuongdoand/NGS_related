#!/usr/bin/env python

import sys
import glob
#from _future_ import print_function
from argparse import ArgumentParser

parser = ArgumentParser(usage = "run_fastqc.py -d path_to_the_input_file -o path_to_the_output_directory -t number_of_threads")
parser.add_argument("-d", help = "The location of the input data", dest = "d", default="./")
parser.add_argument("-f", help = "The format of the file name", dest = "f", default="*.*")
parser.add_argument("-o", help = "The output directory", dest = "o", default = "./")
parser.add_argument("-e", help = "Uncompressed the output file (boolean)", dest = "e", action = "store_true")
parser.add_argument("-t", help = "Threads", dest = "t", default = 1, type = int)

args = parser.parse_args()
for file_name in glob.glob(args.d + args.f):
    print("fastqc %s \\" % file_name)
    print("-o %s \\" % args.o)
    print("-t %d" % args.t, end = " ")
    print("\\\n--extract" if args.e else "", end = "")
    print("\n### job", end = "\n\n")


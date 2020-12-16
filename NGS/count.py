#!/usr/bin/env python

import glob
import sys
from argparse import ArgumentParser

parser = ArgumentParser(usage = "count.py -d path_to_the_input -o path_to_the_output_directory -t number_of_threads")
parser.add_argument("-d", help = "The location of the input data", dest = "d", default="./")
parser.add_argument("-f", help = "The format of the file name", dest = "f", default="*.*")
parser.add_argument("-o", help = "The output directory", dest = "o", default = "./")
parser.add_argument("-t", help = "Threads", dest = "t", default = 1, type = int)

args = parser.parse_args()
for file_name in glob.glob(args.d + args.f):
    if file_name.endswith(".fq"):
        print("cat %s | grep \"+\" | wc -l" % file_name)


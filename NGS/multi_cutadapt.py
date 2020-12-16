#!/usr/bin/env python

import glob
from argparse import ArgumentParser

parser = ArgumentParser(usage = "trim the adapter ")
parser.add_argument("-d", help = "directory of input data", dest = "d", default = "./")
parser.add_argument("-o", help = "output directory", dest = "o", default = "./")
parser.add_argument("-f", help = "format of the input", dest = "f", default = "*.*")

args = parser.parse_args()
for filename in glob.glob(args.d + args.f):
    file_name = filename[2:]
    print("source activate python36")
    print("cutadapt --overlap 3 -f fastq --minimum-length 17 --quality-cutoff 20 -b TGGAATTCTCGGGTGCCAAGG \\")
    print("--untrimmed-output %s%s \\" % (args.o, (file_name[:-3] + "_fragment.fq")))
    print("--too-short-output %s%s \\" % (args.o, (file_name[:-3] + "_short.fq")))
    print("-o %s%s \\" % (args.o, (file_name[:-3] + "_first.fq")))
    print("%s \\" % filename)
    print("1>>%s 2>>%s" % ((args.o + "std_output"), (args.o + "std_error")))
    print("\n### job", end = "\n\n")

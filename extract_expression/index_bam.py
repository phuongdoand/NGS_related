#!/usr/bin/env python

import glob
from argparse import ArgumentParser

parser = ArgumentParser(usage = "Sort the bam files")
parser.add_argument("-d", help = "The location of the input data", dest = "d", default="./")
parser.add_argument("-f", help = "The format of the file name", dest = "f", default="*.*")
parser.add_argument("-o", help = "The output directory", dest = "o", default = "./")
parser.add_argument("-t", help = "Threads", dest = "t", default = 1, type = int)

args = parser.parse_args()

for filename in glob.glob(args.d + args.f):
    print('source activate python36')
    print('samtools index %s' % filename)
    print('\n### job', end = '\n\n')


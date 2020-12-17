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
    print('samtools sort %s -o %s -@ %s \\' % (filename, args.o + filename[2:-4] + '_sort.bam', args.t))
    print('1>>%s 2>>%s' % (args.o + 'std_output', args.o + 'std_error'))
    print('\n### job', end = '\n\n')

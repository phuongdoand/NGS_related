#!/usr/bin/env python

import glob
from argparse import ArgumentParser
from itertools import zip_longest

parser = ArgumentParser('Use to filter the length of miRNA')
parser.add_argument("-d", help = 'Location of the input data', dest = 'd', default = './')
parser.add_argument("-o", help = 'Location of the output data', dest = 'o', default = './')
parser.add_argument("-f", help = 'The input file type', dest = 'f', default = '*.fq')

args = parser.parse_args()

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

for filename in glob.glob(args.d + args.f):
    f = open(args.o + '/' + filename[2:-9] + '_tripped.fq', 'x')
    with open(filename, 'r') as infile:
        for lines in grouper(infile, 4):
            if (len(lines[1]) >= 16 and len(lines[1]) <= 27):
                for line in lines:
                    f.write(line)
    infile.close()
    f.close()

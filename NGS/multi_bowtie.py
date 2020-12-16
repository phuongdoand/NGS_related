#!/usr/bin/env python

import glob
from argparse import ArgumentParser

parser = ArgumentParser(usage = "Mapping the miRNA")
parser.add_argument("-d", help = "directory of input data", dest = "d", default = "./")
parser.add_argument("-o", help = "output directory", dest = "o", default = "./")
parser.add_argument("-f", help = "format of the input", dest = "f", default = "*.*")
parser.add_argument("-p", help = "number of thread", dest = "p", default = 1)

args = parser.parse_args()
for filename in glob.glob(args.d + args.f):
    file_name = filename[2:]
    print('source activate python36')
    print('bowtie -v 1 --norc -p 4 --un %s \\' % ('/project/GP1/phuongdoan/cryocord/miRNA_rawdata/raw_data/trimmed_adapter/clean/tripped/unmapped/unmapGenome_' + file_name))
    print('--al %s \\' % ('/project/GP1/phuongdoan/cryocord/miRNA_rawdata/raw_data/trimmed_adapter/clean/tripped/mapped/mapGenome_' + file_name))
    print('%s %s \\' % (('/project/GP1/phuongdoan/cryocord/miRNA_rawdata/raw_data/trimmed_adapter/clean/tripped/hsa.index'), '/project/GP1/phuongdoan/cryocord/miRNA_rawdata/raw_data/trimmed_adapter/clean/tripped/' + file_name))
    print('-S %s \\' % ('/project/GP1/phuongdoan/cryocord/miRNA_rawdata/raw_data/trimmed_adapter/clean/tripped/mapped/' + file_name[:-3] + '.sam'))
    print('-p %s \\' % args.p)
    print('1>>%s 2>>%s' % ('/project/GP1/phuongdoan/cryocord/miRNA_rawdata/raw_data/trimmed_adapter/clean/tripped/mapped/' + 'std_output', '/project/GP1/phuongdoan/cryocord/miRNA_rawdata/raw_data/trimmed_adapter/clean/tripped/mapped/' + 'std_error'))
    print("\n### job", end = "\n\n")

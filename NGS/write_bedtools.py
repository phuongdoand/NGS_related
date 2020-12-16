#!/usr/bin/env python

import glob

f = open('multi_bedtools.py', 'x')
f.write('source activate python36\n')
f.write('bedtools multicov -bams \\\n')
for filename in glob.glob('./*.bam'):
    f.write('%s \\\n' % filename[2:])
f.write('-bed ./hsa.gff3 > result.txt')
f.close()

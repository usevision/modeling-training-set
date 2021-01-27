#!/usr/bin/env python3

import sys

prev_chr = ''
prev_end = 0
new_file = open('all_PL.bed', 'w+')
for line in open(sys.argv[1]):
	fields = line.strip('\r\n').split('\t')
	if fields[0] != prev_chr:
		prev_chr = fields[0]
		prev_end = fields[2]
		continue
	else:
		new_file.write('{}\t{}\t{}\n'.format(fields[0], prev_end, fields[1]))
		prev_chr = fields[0]
		prev_end = fields[2]


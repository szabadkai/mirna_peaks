#!/usr/bin/python

from sys import argv

fh = open(argv[1])

#this will proess the 6 type blastn output

def read_blast(line):
	dic=dict()
	for i in zip((	'query id', 'subject id', 'identity', 'alignment length', \
					'mismatches', 'gap opens', 'q. start', 'q. end', 's. start',\
					's. end', 'evalue', 'bit score'),line.split()):
		dic[i[0]]=i[1]
	return dic

for line in fh:
	a=read_blast(line)
	if float(a['identity']) > 50 and int(a['alignment length']) > 15:
		print line,
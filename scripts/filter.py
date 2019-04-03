#!/usr/bin/env python
# filter.py
# Filter a fastq file according to a given regex to match to read sequence
# Copyright (C) Darren J Wilkinson, November 2010, http://tinyurl.com/darrenjw
 
import sys,re
 
def readread(s):
    return [s.readline(),s.readline(),s.readline(),s.readline()]
 
def writeread(sread,s):
    for i in range(4):
        s.write(sread[i])
 
def readfilter(query,s,o):
    sread=readread(s)
    while (sread[0]):
        if (query.search(sread[1])):
            writeread(sread,o)
        sread=readread(s)
 
if __name__=='__main__':
    if (len(sys.argv)!=2):
        print "Usage: python filter.py <regex> < infile.fastq > outfile.fastq"
        exit(1)
    regex=re.compile(sys.argv[1])
    readfilter(regex,sys.stdin,sys.stdout)
     
# eof

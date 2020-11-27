#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:19:10 2020

@author: francescoandreace
"""
import time
import sys
tik = time.time()
one = open(str(sys.argv[1]),"r").read()
two = open(str(sys.argv[2]),"r").read()
outf = str(sys.argv[3])

one_sep = one.split("\n")
two_sep = two.split("\n")
if (len(one_sep) != len(two_sep)):
    sys.exit("The two files have different length. Please use the .1 and .2 paired end reads in fasta format.")

for i in range(len(one_sep)):
    if i%2 == 0:
        one_sep[i] = ">r"+str(int(i/2+1))+".1"
        two_sep[i] = ">r"+str(int(i/2+1))+".2"

with open(outf, "w") as out:
    out.write("\n".join(one_sep))
    out.write("\n".join(two_sep))


out.close()
tok = time.time()

print(tok-tik)
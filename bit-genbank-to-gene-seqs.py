#!/usr/bin/env python

from Bio import SeqIO
import argparse
import re
import sys
import os
from collections import defaultdict


def Capilize(seq):
    outSeq = ""
    for i in seq:
        if i == "a":
            outSeq += "A"

        elif i == "g":
            outSeq += "G"

        elif i == "c":
            outSeq += "C"

        elif i == "t":
            outSeq += "T"
    return outSeq


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def customfilter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls.append(string)
    ls = customfilter(ls, [""])
    return ls


def reverseComplement(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString





parser = argparse.ArgumentParser(description="This script takes a genbank file and returns nucelotide sequences for all genes. Genbank file needs to be Genbank/ENA/DDJB compliant (that is, with genes included [--addgenes]). For version info, run `bit-version`.")

required = parser.add_argument_group('required arguments')

required.add_argument("-i", "--input_gb", help='input Genbank file (e.g. "*.gbk", "*.gb", "*.gbff"). Needs to be Genbank/ENA/DDJB compliant.', action="store", dest="input_gb", required=True)
parser.add_argument("-f", "--output_fasta", help='Output fasta file (default: "clean.faa")', action="store", dest="output_fasta", default="clean.faa")

if len(sys.argv)==1:
  parser.print_help(sys.stderr)
  sys.exit(0)

args = parser.parse_args()

input_gb = open(args.input_gb, "r")

output_fasta = open(args.output_fasta, "w")


contigsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
count = 0
seq = ""
locus = ""
for i in input_gb:
    if re.match(r'LOCUS', i):
        locus = (i.rstrip().split(" ")[7])
        seq = ''
    if re.match(r'ORIGIN', i.rstrip()):
        count += 1
    if count > 0:
        pass
        ls = delim(i.rstrip())
        if len(ls) > 0:
            seq += ("".join(ls[1:len(ls)]))
    if re.match(r'//', i.rstrip()):
        count = 0
        seq = Capilize(seq)
        contigsDict[locus] = seq


input_gb = open(args.input_gb, "r")
redundantList = []
coords = ""
for i in input_gb:
    if re.match(r'LOCUS', i):
        locus = (i.rstrip().split(" ")[7])
    if re.match(r'     gene', i.rstrip()):
        coords = (delim(i.rstrip())[1])
        if re.findall(r'complement', coords):
            start = (coords.split("(")[1].split(")")[0].split("..")[0])
            start = remove(start, ["<"])
            end = (coords.split("(")[1].split(")")[0].split("..")[1])
            end = remove(end, [">"])
            seq = (reverseComplement(contigsDict[locus][int(start)-1:int(end)]))

        else:
            start = (coords.split("..")[0])
            start = remove(start, ["<"])
            end = (coords.split("..")[1])
            end = remove(end, [">"])
            seq = contigsDict[locus][int(start)-1:int(end)]

    if re.findall(r'/locus_tag=', i.rstrip()):
        locusTag = (i.rstrip().split("/locus_tag=")[1])
        locusTag = remove(locusTag, ["\""])
        if locusTag not in redundantList:
            redundantList.append(locusTag)
            output_fasta.write(">" + locusTag + " " + locus + " " + coords + "\n")
            output_fasta.write(seq + "\n")
output_fasta.close()
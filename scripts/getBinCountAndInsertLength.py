#!/usr/bin/env python3

# Copyright (C) 2019 Rishvanth Prabakar
#
# Authors: Rish Prabakar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys, argparse
import pysam
from numpy import mean
from bisect import bisect_right, bisect_left

class Bin:
    def __init__(self, line):
        self.chrom, self.start, self.end = line.strip().split()[:3]
        self.start = int(self.start)
        self.end = int(self.end)
        self.count = 0
        self.insertLen = [1]
    def __str__(self):
        return '%s\t%s\t%s\t%d\t%f' % \
            (self.chrom, self.start, self.end, self.count, mean(self.insertLen))
    def increment(self):
        self.count += 1
    def addInsertLen(self, length):
        self.insertLen.append(length)


def loadChromInfo(chromInfoFile):
    chromSizes = {}
    chromNames = []
    for i in open(chromInfoFile):
        chrom, chromLen = i.strip().split()[:2]
        chromLen = int(chromLen)
        chromNames.append(chrom)
        chromSizes[chrom] = chromLen
    chromNames.sort()

    chromOffsets = {}
    totalLength = 0
    for i in chromNames:
        chromOffsets[i] = totalLength
        totalLength += chromSizes[i]

    return chromOffsets, totalLength


def loadBins(chromOffsets, binsFile):
    bins = [Bin(i) for i in open(binsFile)]
    binSorter = [(chromOffsets[bins[i].chrom] + bins[i].start, i)
                    for i in range(len(bins))]
    binSorter.sort()
    # print(binSorter)
    binPosAbs = [i[0] for i in binSorter]
    binLookup = [i[1] for i in binSorter]
    # print(binPosAbs)
    # print(binLookup)
    return bins, binPosAbs, binLookup


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam-file', dest='samFile',
                        required=True, help='paired end sam file')
    parser.add_argument('-b', '--infile-bins', dest='inBinsFile',
                        required=True, help='bins file in 3-col BED format')
    parser.add_argument('-c', '--infile-chroms', dest='inChromFile',
                        required=True, help='2 col file: chrom and chromLen')
    parser.add_argument('-l', '--max-insert-len', dest='maxInsertLen', type=int,
                        required=True, help='maximum insert length')
    parser.add_argument('-o', '--outfile-bins', dest='outBinsFile',
                        required=True, help='output file for bin counts')
    args = parser.parse_args()


    chromOffsets, genomeSize = loadChromInfo(args.inChromFile)
    # print(chromOffsets)
    # print(genomeSize)
    bins, binsPosAbs, binLookup = loadBins(chromOffsets, args.inBinsFile)

    chromsWithBins = dict([(i.chrom, True) for i in bins])
    # print(chromsWithBins)

    totalReads = 0
    samFile = pysam.AlignmentFile(args.samFile, 'r', check_sq=False)
    for read in samFile.fetch():
        if (read.is_proper_pair
            and read.is_read1
            and read.reference_name in chromsWithBins
            and abs(read.template_length) <= args.maxInsertLen):
            # print (read.query_name, read.reference_start,
            #    read.next_reference_start, abs(read.template_length))
            readPosAbs = chromOffsets[read.reference_name] + \
                            read.reference_start 
            # print(readPosAbs)
            binIdx = binLookup[bisect_right(binsPosAbs, readPosAbs) - 1]
            bins[binIdx].increment()
            bins[binIdx].addInsertLen(abs(read.template_length))
            totalReads += 1
    readsPerBin = float(totalReads)/len(bins)

    out = open(args.outBinsFile, 'w')
    for i in range(len(bins)):
        print('%s\t%f' % (str(bins[i]), bins[i].count/readsPerBin), file=out)
    

    samFile.close()
    out.close()

if __name__ == '__main__':
    main()

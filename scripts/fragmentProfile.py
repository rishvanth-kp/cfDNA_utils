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
        self.insertLen = []
        self.shortFragCount = 0
        self.longFragCount = 0
    def __str__(self):
        meanInsertLen = 0
        if (len(self.insertLen) > 0):
           meanInsertLen = mean(self.insertLen)
        fragmentRatio = 0
        if (self.shortFragCount > 0 and self.longFragCount > 0):
            fragmentRatio = self.shortFragCount / self.longFragCount
        return '%s\t%s\t%s\t%d\t%f\t%d\t%d\t%f' % \
                    (self.chrom,
                    self.start,
                    self.end,
                    len(self.insertLen),
                    meanInsertLen,
                    self.shortFragCount,
                    self.longFragCount,
                    fragmentRatio)
    def addInsertLen(self, length):
        self.insertLen.append(length)
    def incrementShortFrag(self):
        self.shortFragCount += 1
    def incrementLongFrag(self):
        self.longFragCount += 1


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
    binPosAbs = [i[0] for i in binSorter]
    binLookup = [i[1] for i in binSorter]

    return bins, binPosAbs, binLookup


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam-file', dest='samFile',
                        required=True, help='paired end sam file')
    parser.add_argument('-b', '--infile-bins', dest='inBinsFile',
                        required=True, help='bins file in 3-col BED format')
    parser.add_argument('-c', '--infile-chroms', dest='inChromFile',
                        required=True, help='2 col file: chrom and chromLen')
    parser.add_argument('-l', '--lower-len-bound', dest='lowerLenBound',
                        type=int, default=100, required=False,
                        help='minimum insert length')
    parser.add_argument('-u', '--upper-len-bound', dest='upperLenBound',
                        type=int, default=220, required=False,
                        help='maximum insert length')
    parser.add_argument('-m', '--mid-len', dest='midLen', type=int, default=150,
                        required=False, help='insert length mid-point')
    parser.add_argument('-o', '--outfile-bins', dest='outBinsFile',
                        required=True, help='output file for fragment profile')
    args = parser.parse_args()


    chromOffsets, genomeSize = loadChromInfo(args.inChromFile)
    bins, binsPosAbs, binLookup = loadBins(chromOffsets, args.inBinsFile)

    chromsWithBins = dict([(i.chrom, True) for i in bins])

    totalReads = 0
    samFile = pysam.AlignmentFile(args.samFile, 'r', check_sq=False)
    for read in samFile.fetch():
        if (read.is_proper_pair
            and read.is_read1
            and read.reference_name in chromsWithBins
            and abs(read.template_length) <= args.upperLenBound
            and abs(read.template_length) >= args.lowerLenBound):

            readPosAbs = chromOffsets[read.reference_name] + \
                            read.reference_start
            binIdx = binLookup[bisect_right(binsPosAbs, readPosAbs) - 1]
            bins[binIdx].addInsertLen(abs(read.template_length))
            if abs(read.template_length) <= args.midLen:
                bins[binIdx].incrementShortFrag()
            else:
                bins[binIdx].incrementLongFrag()

    out = open(args.outBinsFile, 'w')
    for i in range(len(bins)):
        print('%s' % (str(bins[i])), file=out)


    samFile.close()
    out.close()

if __name__ == '__main__':
    main()

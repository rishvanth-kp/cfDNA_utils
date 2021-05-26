#!/usr/bin/env python3

# Copyright (C) 2021 Rishvanth Prabakar
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
from collections import Counter
from numpy import average

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in-file', dest='inFile',
                        required=True, help='paired-end sam file')
    parser.add_argument('-b', '--in-bed-file', dest='inBedFile',
                        required=True, help='bed file')
    parser.add_argument('-q', '--min-mapq', dest='minMapq', type=int,
                        default=0, required=False, help='minimum MAPQ')
    parser.add_argument('-o', '--out-file', dest='outFile',
                        required=True, help='output file')
    parser.add_argument('-s', '--stats-file', dest='statsFile',
                        required=True, help='output stats file')
    args = parser.parse_args()

    inSamFile = pysam.AlignmentFile(args.inFile, 'r')

    maxLen = 0
    insertSzList = []

    for regions in open(args.inBedFile):
        regions = regions.strip().split()
        insertSz = Counter();
        for read in inSamFile.fetch(contig=regions[0],
                        start=int(regions[1]), end=int(regions[2])):
            # Assuming the mapper sets the 'properly paired' flag only
            # for FR reads, and that the TLEN for one alignment is positive
            # and the other is negative.
            # Dangling mates are not tracked, this could potentially lead to
            # missing alignments when one read of a pair is aligned inside a
            # region and the other is aligned outside.
            if (not read.is_supplementary and
                read.mapping_quality >= args.minMapq and
                read.is_proper_pair):
                if (read.template_length > 0):
                    insertSz[read.template_length] += 1
                    if (read.template_length > maxLen):
                        maxLen = read.template_length
        insertSzList.append(insertSz)

    insertSz = Counter();
    for read in inSamFile.fetch():
        if (not read.is_supplementary and
            read.mapping_quality >= args.minMapq and
            read.is_proper_pair):
            if (read.template_length > 0):
                insertSz[read.template_length] += 1
                if (read.template_length > maxLen):
                    maxLen = read.template_length
    insertSzList.append(insertSz)

    out = open(args.outFile, 'w')
    header = 'len'
    for regions in open(args.inBedFile):
        regions = regions.strip().split()
        header += ',%s:%s-%s' % (regions[0], regions[1], regions[2])
    header += ',all'
    print(header, file=out)

    for i in range(maxLen):
        line = str(i)
        for insertSz in insertSzList:
            line += ',%d' % (insertSz[i])
        print(line, file=out)
    out.close()


    statsOut = open(args.statsFile, 'w')
    i = 0
    header = 'region,count,mean,min,max'
    print(header, file=statsOut)
    for regions in open(args.inBedFile):
        regions = regions.strip().split()
        line = '%s:%s-%s' % (regions[0], regions[1], regions[2])
        line += ',%d' % (sum(insertSzList[i].values()))
        line += ',%f' % (average(list(insertSzList[i].keys()),
                         weights=list(insertSzList[i].values())))
        line += ',%d' % (min(insertSzList[i].values()))
        line += ',%d' % (max(insertSzList[i].values()))
        print(line, file=statsOut)
        i += 1

    line = 'all'
    line += ',%d' % (sum(insertSzList[i].values()))
    line += ',%f' % (average(list(insertSzList[i].keys()),
                     weights=list(insertSzList[i].values())))
    line += ',%d' % (min(insertSzList[i].values()))
    line += ',%d' % (max(insertSzList[i].values()))
    print(line, file=statsOut)

    statsOut.close()

if __name__ == '__main__':
    main()

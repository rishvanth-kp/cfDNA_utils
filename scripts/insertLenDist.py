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

# Forward-reverse alignments
# Flag 83:
#   read paired (0x1)
#   read mapped in proper pair (0x2)
#   read reverse strand (0x10)
#   first in pair (0x40)
# Flag 163:
#   read paired (0x1)
#   read mapped in proper pair (0x2)
#   mate reverse strand (0x20)
#   second in pair (0x80)
def isFR(r1Flag, r1Tlen, r2Flag, r2Tlen):
    if (r1Flag == 83 and r1Tlen < 0) and (r2Flag == 163 and r2Tlen >= 0):
        return True
    elif (r1Flag == 163 and r1Tlen >= 0) and (r2Flag == 83 and r2Tlen < 0):
        return True
    elif (r1Flag == 99 and r1Tlen >= 0) and (r2Flag == 147 and r2Tlen < 0):
        return True
    elif (r1Flag == 147 and r1Tlen < 0) and (r2Flag == 99 and r2Tlen >= 0):
        return True
    else:
        return False

# Reverse-forward alignments
# Flags 81 and 161: Similar to 83 and 163, but with 'read mapped in
# proper pair (0x2)' unset.
def isRF(r1Flag, r1Tlen, r2Flag, r2Tlen):
    if (r1Flag == 81 and r1Tlen >= 0) and (r2Flag == 161 and r2Tlen < 0):
        return True
    elif (r1Flag == 161 and r1Tlen < 0) and (r2Flag == 81 and r2Tlen >= 0):
        return True
    elif (r1Flag == 97 and r1Tlen < 0) and (r2Flag == 145 and r2Tlen >= 0):
        return True
    elif (r1Flag == 144 and r1Tlen >= 0) and (r2Flag == 97 and r2Tlen < 0):
        return True
    else:
        return False

# Forward-forward alignments
# Flag 65:
#   read paired (0x1)
#   first in pair (0x40)
# Flag 129:
#   read paired (0x1)
#   second in pair (0x80)
def isFF(r1Flag, r1Tlen, r2Flag, r2Tlen):
    if (r1Flag == 65) and (r2Flag == 129):
        return True
    elif (r1Flag == 129) and (r2Flag == 65):
        return True
    else:
        return False

# Reverse-reverse alignments
# Flag 113:
#   read paired (0x1)
#   read reverse strand (0x10)
#   mate reverse strand (0x20)
#   first in pair (0x40)
# Flag 177:
#   read paired (0x1)
#   read reverse strand (0x10)
#   mate reverse strand (0x20)
#   second in pair (0x80)
def isRR(r1Flag, r1Tlen, r2Flag, r2Tlen):
    if (r1Flag == 113) and (r2Flag == 177):
        return True
    elif (r1Flag == 177) and (r2Flag == 113):
        return True
    else:
        return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--in-file', dest='inFile',
                        required=True, help='paired-end sam file')
    parser.add_argument('-q', '--min-mapq', dest='minMapq', type=int,
                        default=0, required=False, help='minimum MAPQ')
    parser.add_argument('-o', '--out-file', dest='outFile',
                        required=True, help='output file')
    parser.add_argument('-s', '--stats-file', dest='statsFile',
                        required=True, help='output stats file')
    args = parser.parse_args()

    inSamFile = pysam.AlignmentFile(args.inFile, 'r')

    lonelyMates = {}
    readCount = 0
    otherCount = 0
    flagStat = Counter()
    frInsertSz = Counter()
    rfInsertSz = Counter()
    ffInsertSz = Counter()
    rrInsertSz = Counter()

    for read in inSamFile.fetch():
        readCount += 1
        if (not read.is_supplementary and
            read.mapping_quality >= args.minMapq):
            if read.query_name in lonelyMates:
                r1 = read
                r2 = lonelyMates.pop(read.query_name)
                flagStat[r1.flag] += 1
                flagStat[r2.flag] += 1
                if isFR(r1.flag, r1.template_length,
                        r2.flag, r2.template_length):
                    frInsertSz[abs(r1.template_length)] += 1

                elif isRF(r1.flag, r1.template_length,
                          r2.flag, r2.template_length):
                    rfInsertSz[abs(r1.template_length)] += 1

                elif isFF(r1.flag, r1.template_length,
                          r2.flag, r2.template_length):
                    ffInsertSz[abs(r1.template_length)] += 1

                elif isRR(r1.flag, r1.template_length,
                          r2.flag, r2.template_length):
                    rrInsertSz[abs(r1.template_length)] += 1

                else:
                    otherCount += 1

            else:
                lonelyMates[read.query_name] = read

    out = open(args.outFile, 'w')
    for i in range(max(frInsertSz.keys())):
        line = str(i)
        line += ',%d,%d,%d,%d' % (frInsertSz[i], rfInsertSz[i],
                 ffInsertSz[i], rrInsertSz[i])
        print(line, file=out)
    out.close()

    statsOut = open(args.statsFile, 'w')

    statsOut = close()

    print(len(lonelyMates))
    print(flagStat)
    print(readCount)
    print(max(frInsertSz.keys()))
    print(max(rfInsertSz.keys()))
    print(max(ffInsertSz.keys()))
    print(max(rrInsertSz.keys()))
    print(otherCount)

if __name__ == '__main__':
    main()

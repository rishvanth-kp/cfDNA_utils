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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--in-sam-file', dest='inSamFile',
                        required=True, help='paired-end sam file')
    parser.add_argument('-l', '--min-insert-len', dest='minInsertLen', type=int,
                        required=True, help='min insert len to output')
    parser.add_argument('-u', '--max-insert-len', dest='maxInsertLen', type=int,
                        required=True, help='max insert len to output')
    parser.add_argument('-o', '--out-sam-file', dest='outSamFile',
                        required=True, help='output sam file')
    args = parser.parse_args()

    inSamFile = pysam.AlignmentFile(args.inSamFile, 'r', check_sq=False)
    outSamFile = pysam.AlignmentFile(args.outSamFile, 'w', template=inSamFile)

    for read in inSamFile.fetch():
        if (abs(read.template_length) >= args.minInsertLen
            and abs(read.template_length) <= args.maxInsertLen):
            outSamFile.write(read)

if __name__ == '__main__':
    main()

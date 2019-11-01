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

import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample-name', dest='sampleName',
                        required=True, help='sample name')
    parser.add_argument('-m', '--mapped-stat', dest='inMappedStat',
                        required=True, help='all mapped flagstat')
    parser.add_argument('-u', '--unique-stat', dest='inUniqueStat',
                        required=True, help='uniquely mapped flagstat')
    args = parser.parse_args()

    mappedStatFile = open(args.inMappedStat, 'r')
    uniqueStatFile = open(args.inUniqueStat, 'r')

    mappedStats = mappedStatFile.readlines()
    uniqueStats = uniqueStatFile.readlines()

    totalReads = int(mappedStats[0].split()[0])
    mappedReads = int(mappedStats[4].split()[0])
    uniqueReads = int(uniqueStats[4].split()[0])

    print ('%s\t%d\t%d\t%f\t%d\t%f' %
            (args.sampleName,
            totalReads,
            mappedReads,
            mappedReads / totalReads,
            uniqueReads,
            uniqueReads / totalReads))

    mappedStatFile.close()
    uniqueStatFile.close()

if __name__ == '__main__':
    main()

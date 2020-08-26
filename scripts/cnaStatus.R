#!/usr/bin/env Rscript

# Copyright (C) 2020 Rish Prabakar
#
# Authors: Rish Prabakar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

main <- function() {
  
  args = commandArgs(trailingOnly=TRUE)
  if (length(args) != 5) {
    stop("cnaStatus.R <cna.csv> <chrom-regions.csv> <gain-thresold>
          <loss-threshold> <outfile.csv>", call.=FALSE)
  }

  # col 1-3 (header): chr (chrom), chr pos (chrompos), 
  #   absolute chr pos (abspos).
  # col 4- (header): one column per sample (sample name)
  # rows: segmented bin values, one row per bin
  cna <- read.csv(args[1], header = TRUE)
  # col 1-4 (no header): chr, chr start, chr end, region name
  # one row per region. 
  regions <- read.csv(args[2], header = FALSE)
  names(regions) <- c("chr", "chr.start", "chr.end", "name")
  # gain and loss threshold 
  gain.th <- as.numeric(args[3])
  loss.th <- as.numeric(args[4])
  # csv file to output gain/loss status for each region
  outfile <- args[5]

  print(dim(cna))
  print(gain.th)
  print(loss.th)  
  n.samples <- ncol(cna) - 3

  for (i in 1 : nrow(regions)) {
    print(regions[i,])
  }

}

main()

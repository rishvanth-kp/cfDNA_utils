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
          <loss-threshold> <outfile-prefix>", call.=FALSE)
  }

  # col 1-3 (header): chr (chrom), chr pos (chrompos),
  #   absolute chr pos (abspos).
  # col 4- (header): one column per sample (sample name)
  # rows: segmented bin values, one row per bin
  cna <- read.csv(args[1], header = TRUE)
  # col 1-4 (no header): region name, chr, chr start, chr end
  # one row per region. chr is encoded as 1-24 to keep it compatible
  # with the segmented csv file from the CNV app
  # (using a csv file and encoding chromosomes are 1-24 are both bad
  # designs. Ideally, a bed file should be used and chromosomes should 
  # be encoded using the standard UCSC format.)
  regions <- read.csv(args[2], header = FALSE)
  names(regions) <- c("name", "chrom", "chrom.start", "chrom.end")
  # gain and loss threshold
  gain.th <- as.numeric(args[3])
  loss.th <- as.numeric(args[4])
  # outfile prefix
  outfile.prefix <- args[5]

  n.samples <- ncol(cna) - 3
  sample.names <- names(cna[,4:ncol(cna)])

  out.median <- c()
  out.status <- c()
  for (i in 1 : nrow(regions)) {
    cna.in.region <- cna[which(cna$chrom == regions[i,]$chrom &
                               cna$chrompos >= regions[i,]$chrom.start &
                               cna$chrompos < regions[i,]$chrom.end),
                               4:ncol(cna)]

    region.median <- sapply(cna.in.region, median)
    
    if (nrow(cna.in.region) == 0) {
      region.status <- rep("NA", times=n.samples) 
    }
    else {
      region.status <- rep("neutral", times=n.samples)
      region.status[region.median >= gain.th] <- "gain"
      region.status[region.median <= loss.th] <- "loss"
    }
  
    out.median <- rbind(out.median, region.median)
    out.status <- rbind(out.status, region.status)
  }

  out.status <- as.data.frame(out.status, row.names = FALSE)
  names(out.status) <- sample.names
  out.status <- cbind(regions, out.status)
  write.csv(out.status, quote = FALSE, row.names = FALSE,
            file = sprintf("%s_status.csv", outfile.prefix))

  out.median <- as.data.frame(out.median, row.names = FALSE)
  names(out.median) <- sample.names
  out.median <- cbind(regions, out.median)
  write.csv(out.median, quote = FALSE, row.names = FALSE,
            file = sprintf("%s_median.csv", outfile.prefix))
}

main()

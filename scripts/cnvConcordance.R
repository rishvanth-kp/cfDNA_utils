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
  if (length(args) != 4) {
    stop ("cnaConcordance.R <cna-tumor> <cna-rb>
            <gain-threshold> <loss-threshold>", call.=FALSE)
  }

  # two column csv file with chrom and segmented values
  # First row has the sample name
  cna.tumor <- read.csv(args[1], header=T)
  cna.rb <- read.csv(args[2], header=T)
  gain.thresh <- as.numeric(args[3])
  loss.thresh <- as.numeric(args[4])

  names(cna.tumor) <- c("chrom", "seg")
  names(cna.rb) <- c("chrom", "seg")

  cna.tumor <- cna.tumor[cna.tumor$chrom != 23,]
  cna.tumor <- cna.tumor[cna.tumor$chrom != 24,]
  cna.rb <- cna.rb[cna.rb$chrom != 23,]
  cna.rb <- cna.rb[cna.rb$chrom != 24,]

  print(dim(cna.tumor))
  print(dim(cna.rb))

  cna.ratio <- (cna.tumor$seg / cna.rb$seg)
  concordance <- (length(cna.ratio[cna.ratio >= loss.thresh &
                    cna.ratio <= gain.thresh]) / length(cna.ratio))
  concordance <- concordance * 100

  print(concordance)
}

main()

#!/usr/bin/env Rscript

# Copyright (C) 2021 Rish Prabakar
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

library(reshape2)
library(ggplot2)

main <- function () {

  args = commandArgs(trailingOnly=TRUE)
  if (length(args) != 2) {
    stop("plotBinnedInsertLenDist.R <insert-len-csv> <sample-name>",
          call.=FALSE)
  }

  inserts <- read.csv(args[1], header=TRUE)
  sample <- args[2]

  # The last col had the distribution of all the reads
  # Not plotting it to avoid dwarfing the others
  inserts.melt <- melt(inserts[,1:ncol(inserts)-1], id.vars="len")
  ggplot(data=inserts.melt) +
    geom_line(mapping=aes(x=len, y=value, color=variable),
      size=0.3, alpha=0.5) +
    labs(x="Length", y="Read count", color="Region",
      title=sample) +
    theme_bw()
  ggsave(sprintf("%s_alt_ins_sz.pdf", sample), height=4, width=6)


  for (i in 2:ncol(inserts)) {
    inserts[,i] <- inserts[,i]/sum(inserts[,i])
  }
  inserts.melt <- melt(inserts, id.vars="len")
  ggplot(data=inserts.melt) +
    geom_line(mapping=aes(x=len, y=value, color=variable),
      size=0.3, alpha=0.5) +
    labs(x="Length", y="Normalized read count", color="Region",
      title=sample) +
    theme_bw()
  ggsave(sprintf("%s_alt_ins_sz_norm.pdf", sample), height=5, width=6)

}

main ()

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

main <- function () {

  args = commandArgs(trailingOnly=TRUE)
  if (length(args) != 5) {
    stop("insertLenDist <picard-insert-len> <min-len> <max-len>
          <sample-name> <pdf-name>", call.=FALSE)
  }

  inserts <- read.table(args[1], skip=11)
  names(inserts) <- c("len", "reads")
  min.len <- as.numeric(args[2])
  max.len <- as.numeric(args[3])
  sample.name <- args[4]
  pdf.name <- args[5]

  max.insert.len <- max(inserts$len)
  if (max.len > max.insert.len) {
    max.len <- max.insert.len
  }

  inserts <- inserts[inserts$len >= min.len & inserts$len <= max.len,]

  max.y <- max(inserts$reads/sum(inserts$reads))
  print(max.y)

  pdf(pdf.name, width=4, height=4, useDingbats=FALSE)
  par(pin=c(2,2))
  plot(inserts$len, inserts$reads / sum(inserts$reads), type="n",
       xlim=c(min.len, max.len), ylim=c(0, max.y),
       xlab="Fragment length (bp)", ylab="Normalized read count",
       main=sample.name)
  lines(inserts$len, inserts$reads / sum(inserts$reads), col="red")
  dev.off()
}

main ()

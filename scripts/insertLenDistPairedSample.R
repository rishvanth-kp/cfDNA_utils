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
    stop("insertLenDist <picard-insert-len-1> <picard-insert-len-2>
          <min-len> <max-len> <sample-name>", call.=FALSE)
  }

  inserts.s1 <- read.table(args[1], skip=11)
  names(inserts.s1) <- c("len", "reads")
  inserts.s2 <- read.table(args[2], skip=11)
  names(inserts.s2) <- c("len", "reads")
  min.len <- as.numeric(args[3])
  max.len <- as.numeric(args[4])
  sample.name <- args[5]

  max.insert.len <- max(inserts.s1$len, inserts.s2$len)
  if (max.len > max.insert.len) {
    max.len <- max.insert.len
  }

  inserts.s1 <- inserts.s1[inserts.s1$len >= min.len &
                           inserts.s1$len <= max.len,]
  if (min(inserts.s1$len) > min.len) {
    inserts.s1 <- rbind(c(min.len, 0), inserts.s1)
  }
  if (max(inserts.s1$len) < max.len) {
    inserts.s1 <- rbind(inserts.s1, c(max.len, 0))
  }

  inserts.s2 <- inserts.s2[inserts.s2$len >= min.len &
                           inserts.s2$len <= max.len,]
  if (min(inserts.s2$len) > min.len) {
    inserts.s2 <- rbind(c(min.len, 0), inserts.s2)
  }
  if (max(inserts.s2$len) < max.len) {
    inserts.s2 <- rbind(inserts.s2, c(max.len, 0))
  }


  max.y <- max(inserts.s1$reads/sum(inserts.s1$reads),
               inserts.s2$reads/sum(inserts.s2$reads))

  pdf(sprintf("%s.pdf", sample.name), width=4, height=4, useDingbats=FALSE)
  par(pin=c(2,2))
  plot(inserts.s1$len, inserts.s1$reads / sum(inserts.s1$reads), type="n",
       xlim=c(min.len, max.len), ylim=c(0, max.y),
       xlab="Fragment length (bp)", ylab="Normalized read count",
       main=sample.name)
  lines(inserts.s1$len, inserts.s1$reads / sum(inserts.s1$reads), col="red")
  lines(inserts.s2$len, inserts.s2$reads / sum(inserts.s2$reads), col="blue")
  dev.off()
}

main ()

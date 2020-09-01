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
  if (length(args) != 4) {
    stop("insertLenDistMerged <sample-list.csv> <min-len> <max-len>
          <sample-name>", call.=FALSE)
  }

  # 2 column headerless csv file. 1st col: sameple type.
  # 2nd col: path to picard output
  sample.list <- read.csv(args[1], header=FALSE)
  names(sample.list) <- c("type", "path")
  # min and max insert length to plot
  min.len <- as.numeric(args[2])
  max.len <- as.numeric(args[3])
  sample.name <- args[4]

  sample.types <- unique(sample.list$type)
  inserts.mean <- data.frame(len=seq(min.len, max.len))
  for (i in sample.types) {
    inserts.type <- data.frame(len=seq(min.len, max.len))
    sample.count <- 0
    for (j in 1 : nrow(sample.list)) {
      if (sample.list[j,]$type == i) {
        inserts <- read.table(sample.list[j,]$path, skip=11)
        names(inserts) <- c("len", "read.count")

        inserts <- inserts[inserts$len >= min.len & inserts$len <= max.len,]
        inserts$read.count <- inserts$read.count / sum(inserts$read.count)

        inserts.type <- merge(inserts.type, inserts, by="len", all=TRUE)
        names(inserts.type)[ncol(inserts.type)] <- sprintf("%s.%d",
                                                  i, sample.count)

        sample.count <- sample.count + 1
      }
    }

    inserts.type[is.na(inserts.type)] <- 0
    type.mean <- rowMeans(inserts.type[,2:ncol(inserts.type)])

    inserts.mean <- cbind(inserts.mean, type.mean)
    names(inserts.mean)[ncol(inserts.mean)] <- i

  }

  max.y <- max(inserts.mean[,2:ncol(inserts.mean)])
  line.color <- c("red", "blue", "black", "gray45")
  pdf(sprintf("%s.pdf", sample.name), width=4, height=4, useDingbats=FALSE)
  par(pin=c(1.5,1.5))
  plot(inserts.mean[,1:1], inserts.mean[,2:2], type="n",
       xlim=c(min.len, max.len), ylim=c(0, max.y),
       xlab="Fragment length (bp)", ylab="Normalized read count",
       main=sample.name)
  for (i in 1 : length(sample.types)) {
    lines(inserts.mean[,1:1], inserts.mean[,(i+1):(i+1)], col=line.color[i])
  }
  dev.off()

}

main()

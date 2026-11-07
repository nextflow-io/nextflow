#!/usr/bin/env Rscript

# Contract fixture for native external R scripts.
#
# This file should be valid R as written. Nextflow injects an idiomatic runtime
# object from the task sidecar; this fixture uses S4-style access by analogy
# with Snakemake while leaving the final API open to implementation.

count_records <- function(fasta) {
    lines <- readLines(fasta)
    sum(startsWith(lines, ">"))
}

reads <- nextflow@input[["reads"]]
summary <- nextflow@output[["summary"]]
context <- nextflow@output[["context"]]
prefix <- nextflow@params[["prefix"]]

writeLines(sprintf("%s\trecords\t%s", prefix, count_records(reads)), summary)
writeLines(c(
    "language=r",
    sprintf("reads=%s", reads),
    sprintf("summary=%s", summary),
    sprintf("prefix=%s", prefix),
    sprintf("cpus=%s", nextflow@cpus),
    sprintf("attempt=%s", nextflow@attempt)
), context)

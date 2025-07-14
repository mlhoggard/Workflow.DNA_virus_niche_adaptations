#!/usr/bin/env Rscript

# MIT License
# Copyright (c) 2025 Michael Hoggard
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Load arguments
args = commandArgs(trailingOnly=TRUE)

# Load libraries
library(readr)
library(dplyr)
library(qvalue)

# Load file, calculate qvalues, calculate fdr adjusted p-values, sort by p.adj
## Filter with parameters: score >= 0.9 & pvalue <= 0.05 & p.adj <= 0.1
## export as tsv (overwrite original)
result <- read_tsv(args[1]) %>%
    mutate(., qvalue = tryCatch(qvalue(.$pvalue, pi0.meth="bootstrap")$qvalues, error=function(e) "error")) %>%
    mutate(., p.adj = p.adjust(.$pvalue, method="fdr")) %>%
    arrange(., p.adj) %>%
    filter(., (score >= 0.9 & pvalue <=0.05 & p.adj <= 0.1)) %>%
    write_tsv(., args[1])

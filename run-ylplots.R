source("R/rearrangement_plot.R")
library(stringr)
library(parallel)

make_plot <- function(in_row, out_dir, sample_annot, chrom_sizes) {
    fn_bedpe <- in_row[2]
    fn_cn <- in_row[3]
    this_sample <- in_row[1]
    annot <- na.omit(sample_annot[sample_annot$CASM.sample.name == this_sample, c("chr", "start", "end", "Guide")])
    if (!file.exists(fn_bedpe)) {
        print(paste0("file not found ", fn_bedpe))
        return()
    }
    if (!file.exists(fn_cn)) {
        print(paste0("file not found ", fn_cn))
        return()
    }
    bedpe <- read.table(fn_bedpe, header = F, sep = "\t", stringsAsFactors = F)
    cn <- read.table(fn_cn, header = F, sep = "\t", stringsAsFactors = F,
        colClasses = c("factor", "numeric", "numeric", "numeric"))
    chr_lens <- read.table(chrom_sizes, header = F, sep = "\t",
        row.names = 1, colClasses = c("character", "character", "numeric"))
    temp <- rownames(chr_lens)
    chr_lens <- chr_lens[, 2]
    names(chr_lens) <- temp
    pdf(paste0(out_dir, "/", this_sample, ".pdf"), h = 4, w = 10)
    for (c in paste0("chr", c(seq(1, 22), "X", "Y"))) {
        plot_rearrangements(bedpe = bedpe, chrs = c, chr_lens = chr_lens, cn_bedgraph = cn,
        cn_win_size = 1e4, cn_cex = 0.3, ref = "hg38", ideogram = T, annot = annot)
    }
    dev.off()
    warnings()
}

args <- commandArgs(trailingOnly = TRUE)
input_sample_list <- args[1]
output_dir <- args[2]
guide_list <- args[3]
cut_sites_file <- args[4]
chrom_sizes <- args[5]

sample_cut_sites <- read.csv(guide_list, na.strings = c(" - "))
cut_sites <- read.csv(cut_sites_file)
sample_annot <- merge(sample_cut_sites, cut_sites[cut_sites$Mismatches. == 0, ], by.x = "Guide", by.y = "Cut.Site", all.x = TRUE)
sample_annot[c("chr", "start", "end")] <- str_split_fixed(sample_annot$Location, "(:|-)", 3)
sample_annot <- sample_annot[c("CASM.sample.name", "chr", "start", "end", "Guide")]
sample_annot$start <- as.integer(sample_annot$start)
sample_annot$end <- as.integer(sample_annot$end)

inputs <- read.table(input_sample_list, header = F, sep = "\t")
apply(inputs, 1, make_plot, output_dir, sample_annot, chrom_sizes)

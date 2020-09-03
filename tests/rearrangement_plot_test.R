library(rearrangement.plot)

fn_bedpe <- '/lustre/scratch116/casm/cgp/pancancer/downloads/donors/BOCA-UK/CGP_donor_1691124/sanger/fc96f0df-ef04-411e-e040-11ac0d484291_SV_4_0_12/fc96f0df-ef04-411e-e040-11ac0d484291_vs_fc96f0df-ef01-411e-e040-11ac0d484291.annot.bedpe'
fn_cn <- '/lustre/scratch116/casm/cgp/pancancer/downloads/donors/BOCA-UK/CGP_donor_1691124/sanger/fc96f0df-ef04-411e-e040-11ac0d484291_SV_4_0_12/intermediates/fc96f0df-ef04-411e-e040-11ac0d484291.ngscn.abs_cn.bg.gz'
zz_cn <- gzfile(fn_cn, 'rt')  
bedpe <- read.table(fn_bedpe, header = F, sep = "\t", stringsAsFactors = F)
cn    <- read.table(zz_cn,    header = F, sep = "\t", stringsAsFactors = F, colClasses = c("factor", "numeric", "numeric", "numeric"))

chr_lens <- t(as.vector(read.table('~/dev/params/hg19.chrom_sizes.txt', header=F, sep="\t", row.names=1, colClasses = c("character","NULL", "numeric")))

this_sample <- 'test'
fpath_out <- '~/dev/testplot'
pdf(paste0(fpath_out,".pdf"), h = 4, w = 10)
plot_rearrangements(bedpe = bedpe, chrs=c(1), chr_lens=chr_lens, cn_bedgraph = cn, cn_win_size = 1e4, cn_cex = 0.3)
dev.off()


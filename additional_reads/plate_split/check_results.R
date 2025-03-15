library(snpR);library(ggplot2)

setwd("additional_reads/plate_split/")
meta <- read.table("sample_metadata.txt", header = T)
region_table <- data.frame(pop = 1:14,
                           region = c(rep("Shasta-Trinity", 2),
                                      rep("SixRivers",5),
                                      rep("Plumas", 3),
                                      rep("Shasta-Trinity", 3),
                                      "Mendicino"))
meta$region <- region_table[match(meta$pop, region_table$pop), 2]
# ngsadmix <- plot_structure(x = ".qopt", k = 1:5, reps = 10, qsort_K = 2,
#                            clumpp_path = "C://usr/bin/CLUMPP.exe", 
#                            clumpp.opt = "large.k.greedy", facet = as.character(meta$region), 
#                            facet.order = c("Shasta-Trinity", "SixRivers", "Plumas", "Mendicino"))

# setwd("..")
genotypes <- data.table::fread("genotypes.geno.gz")
rmcol <- ncol(genotypes)
genotypes <- genotypes[,-..rmcol]
snp.meta <- genotypes[,1:2]
colnames(snp.meta) <- c("scaffold", "position")
dat <- import.snpR.data(genotypes[,-c(1:2)], snp.meta = snp.meta, sample.meta = meta)

dat <- dat[region = -c("Mendicino")]
dat <- filter_snps(dat, hwe = .05, fwe_method = "holm", hwe_facets = "pop", remove_garbage = .25, min_ind = .5, min_loci = .5)

dat <- calc_pairwise_fst(dat, c("pop"), verbose = TRUE)
plot_pairwise_fst_heatmap(dat, c("pop"), mark_sig = .05)

dat <- calc_pi(dat, c("pop", "region"))
dat <- calc_ho(dat, c("pop", "region"))
dat <- calc_he(dat, c("pop", "region"))
dat <- calc_fis(dat, c("pop", "region"))
dat <- calc_hs(dat, c("pop", "region"))
dat <- calc_private(dat, c("pop", "region"))
res.means <- get.snpR.stats(dat, c("region"), c("ho", "pi", "fis", "private", "hs", "he"))
res.means$weighted.means
ggplot(res.means$sample, aes(x = region, y = hs)) + geom_violin() + geom_jitter(height = 0, width = .2)


pca <- plot_clusters(dat, "region")

saveRDS(list(NGSplot = ngsadmix$plot, statistics = res.means, fst_matrix = fstm, private_alleles = pa, snpRdata = dat), "results_for_cody.RDS")

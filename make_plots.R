library(ggplot2); library(snpR); library(plotly); library(sf); library(ggnewscale)
# PCA
mat <- read.table("IBS_mat.ibsMat")
meta <- read.table("sample_metadata.txt", header = T)
# scale and center
mat <- scale(mat)
mat[is.nan(mat)] <- 0

## remove NAs
# good.samps <- !is.na(mat[,1])
# mat <- mat[good.samps,]
# mat <- mat[,good.samps]
# meta <- meta[good.samps,]
# row.nas <- rowSums(is.na(mat))
# row.nas <- !as.logical(row.nas)
# mat <- mat[row.nas, row.nas]
# meta <- meta[row.nas,]
## pca
pca_r <- svd(mat)
colnames(pca_r$u) <- paste0("PC", 1:ncol(pca_r$u))
pca <- cbind(meta, pca_r$u)

# bind more metadata info
pop_group_tab <- data.frame(pop = 1:14,
                            location = c(rep("Shasta-Trinity", 2),
                                         rep("Six Rivers", 5),
                                         rep("Plumas", 3),
                                         rep("Shasta-Trinity", 3),
                                         "Mendocino"))
pca$location <- pop_group_tab[match(pca$pop, pop_group_tab$pop),2]

loadings <- (pca_r$d^2)/sum(pca_r$d^2)
loadings <- round(loadings * 100, 2)

ggplot(pca, aes(x = PC1, y = PC2, color = as.factor(location))) + 
  geom_point() + theme_bw() +
  scale_color_viridis_d() +
  scale_shape_discrete(solid = TRUE)


## umap
umap_r <- umap::umap(mat, n_epochs = 1000, verbose = T)
colnames(umap_r$layout) <- paste0("PC", 1:ncol(umap_r$layout))
umap_r$layout <- cbind(umap_r$layout, meta, location = pca$location)
umap_r$layout$location <- factor(umap_r$layout$location, c("Shasta-Trinity", "Six Rivers", "Plumas", "Mendocino"))

umap_p <- ggplot(as.data.frame(umap_r$layout), aes(x = PC1, y = PC2, color = as.factor(location))) + 
  geom_point(size = 4) + theme_bw() +
  khroma::scale_color_okabeito() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)) +
  guides(color = guide_legend(title = "Location"))

ggsave("paper/Figure_4.pdf", umap_p, device = "pdf", width = 15, height = 11)


# missing value imputation with missMDA
# nb <- estim_ncpPCA(mat,ncp.max=10, method.cv = "Kfold")
# ipca <- imputePCA(mat, nb$ncp)
# ## pca
# pca_r <- prcomp(ipca$completeObs)
# pca <- as.data.frame(pca_r$x) #grab the PCA vectors.
# pca <- cbind(pca, meta)
# loadings <- (pca_r$sdev^2)/sum(pca_r$sdev^2)
# loadings <- round(loadings * 100, 2)
# ggplot(pca, aes(x = PC1, y = PC2, color = as.factor(pop))) + theme_bw() +
#   scale_color_viridis_d() + geom_text(aes(label = pop)) + 
#   ggplot2::xlab(paste0("PC1 (", loadings[1], "%)")) + 
#   ggplot2::ylab(paste0("PC2 (", loadings[2], "%)"))
# ## umap
# umap_r <- umap::umap(ipca$completeObs, n_epochs = 1000, verbose = T)
# colnames(umap_r$layout) <- paste0("PC", 1:ncol(umap_r$layout))
# umap_r$layout <- cbind(umap_r$layout, meta)
# umap_r$layout$missing <- rowSums(is.na(mat))
# 
# ggplot(as.data.frame(umap_r$layout), aes(x = PC1, y = PC2, color = as.factor(pop))) + theme_bw() +
#   scale_color_viridis_d() + geom_text(aes(label = pop)) +
#   ggplot2::xlab("Dim 1") + ggplot2::ylab("Dim 2")


# NGSrelate
setwd("additional_reads/plate_split/NGSadmix/")
meta <- read.table("../sample_metadata.txt", header = T)
meta$location <- pop_group_tab[match(meta$pop, pop_group_tab$pop),2]
ap <- khroma::color(palette = "okabeito")
ap <- as.character(ap(4))[c(2, 1, 3, 4)]
ngs <- plot_structure(".qopt", meta$location, k = 1:4, 
                      reps = 10, clumpp.opt = "large.k.greedy", clumpp = T, 
                      alt.palette = ap, 
                      facet.order = c("Shasta-Trinity", "Six Rivers", "Plumas", "Mendocino"), 
                      t.sizes = c(20, 20, 18))
ngs <- ngs$plot
ngs <- ngs + xlab("Location") +
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

leg <- ggpubr::get_legend(ngs)
# ngs <- ngs + guides(fill = "none", color = "none")


# prep cv storage
cv_storage <- expand.grid(K = 1:10, rep = 1:10)
cv_storage$est_ln_prob <- 0

# run for each k and rep
for(j in 1:10){
  for(i in 1:10){
    # run structure
    outfile <- paste0("full_K", i, "_r", j, ".log")
    
    # grab ln summary data
    suppressWarnings(lndat <- readLines(outfile))
    sline <- grep("best like", lndat)
    lndat <- lndat[sline]
    lndat <- as.numeric(gsub(" after.+", "", gsub("^.+=", "", lndat)))
    cv_storage[which(cv_storage$K == i & cv_storage$rep == j),"est_ln_prob"] <- lndat
  }
}

cv_storage <- as.data.table(cv_storage)

evanno <-cv_storage[, mean(est_ln_prob), by = K]
colnames(evanno)[2] <- "mean_est_ln_prob"
evanno$lnpK <- NA
evanno$lnppK <- NA
evanno$deltaK <- NA
evanno$sd_est_ln_prob <- cv_storage[, sqrt(stats::var(est_ln_prob)), by = K][[2]]
evanno$lnpK[-1] <- evanno$mean_est_ln_prob[-1] - evanno$mean_est_ln_prob[-nrow(evanno)]
evanno$lnppK[-nrow(evanno)] <- abs(evanno$lnpK[-nrow(evanno)] - evanno$lnpK[-1])
# evanno$deltaK[-c(1, nrow(evanno))] <- abs((evanno$mean_est_ln_prob[-1][-1] - 
#                                             2*evanno$mean_est_ln_prob[-c(1, nrow(evanno))] +
#                                             evanno$mean_est_ln_prob[-nrow(evanno)][-(nrow(evanno) - 1)])/
#                                           evanno$sd_est_ln_prob[-c(1, nrow(evanno))])
evanno$deltaK[-c(1, nrow(evanno))] <- abs(evanno$lnppK)[-c(1, nrow(evanno))]/evanno$sd_est_ln_prob[-c(1, nrow(evanno))] # no reason to resolve for ln''(K)

infs <- which(is.infinite(evanno$deltaK))
if(length(infs) > 0){
  evanno$deltaK[infs] <- NA
  warning(paste0("For some values of K (", paste0((k)[infs], collapse = ", "), "), all reps had the same estimated ln(likelihood). Since calculating deltaK involves dividing by the standard deviation of the ln(likelihood) estimates across reps, this will return 'Inf', and have thus been assigned a deltaK of NA."))
}

evplot <- ggplot(evanno[2:9], aes(x = K, y = deltaK)) + geom_point() + geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 2:9) +
  ylab(expression(Delta*K)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

# F6 <- gridExtra::grid.arrange(evplot + ggtitle("A") + theme(title = element_text(size = 24)), 
#                         ngs + ggtitle("B") + theme(plot.margin = unit(c(.1, .1, .1, 2), "cm"), title = element_text(size = 24)), 
#                         widths = c(.5, 1))


# NGS relate map
coords <- c(40.962944, -122.795056,
            40.985278, -122.773750,
            41.850472, -123.907667,
            41.922139, -123.867833,
            41.892861, -123.864389,
            41.872944, -123.885417,
            41.856278, -123.901306,
            40.013028, -121.126556,
            40.013278, -121.126306,
            40.001722, -120.978694,
            41.313472, -122.421833,
            41.406528, -122.522222,
            41.385639, -122.535306,
            39.253472, -123.746528)
coords <- matrix(coords, ncol = 2, byrow = TRUE)
coords <- as.data.frame(coords)
coords$pop <- as.character(1:14)
colnames(coords)[1:2] <- c("lat", "long")
psf <- sf::st_as_sf(as.data.frame(coords), coords = c("long", "lat"))
psf <- sf::`st_crs<-`(psf, "EPSG:4326")
pop <- as.character(meta$pop)
assn <-plot_structure(".qopt", pop, k = 3, 
                              reps = 10, clumpp.opt = "large.k.greedy", clumpp = T, 
                              alt.palette = ap, 
                              facet.order = as.character(1:14), 
                              t.sizes = c(20, 20, 18))

background <- rnaturalearth::ne_states(iso_a2 = "US", returnclass = "sp")
background <- sf::st_as_sf(background)
background <- background[background$name %in% c("California", "Nevada", "Oregon"),]
extent <- sf::st_bbox(background)
extent["ymin"] <- 38
extent["ymax"] <- 42.5
extent["xmax"] <- -119
background <- sf::st_crop(background, extent)
psf <- sf::st_transform(psf, sf::st_crs(background))

# veg data, from (Ecological Subsections: Landcover)[https://data.fs.usda.gov/geodata/edw/datasets.php?dsetCategory=biota]
veg <- st_read("../../../shpfiles/S_USA.NationalLandcoverSubsections/")
veg <- st_transform(veg, st_crs(background))
sf::sf_use_s2(FALSE)
veg <- st_crop(veg, extent)
veg$`% Forest Cover` <- veg$EVERGREEN1 + veg$DECID_FO_1 + veg$MIXED_FO_1



map <- plot_structure_map(assn, k = 3, facet = "pop", pop_coordinates = psf,
                                    layers = list(geom_sf(data = veg, aes(fill = `% Forest Cover`), color = NA, lwd = 0, alpha = .5),
                                                  khroma::scale_fill_batlow(reverse = TRUE),
                                                  # guides(fill = guide_legend(title = "% Forest Cover")),
                                                  geom_sf(data = background, fill = NA, lwd = 1, color = "black"),
                                                  ggnewscale::new_scale_fill()),
                                    radius_scale = .1,
                                    alt.palette = ap[1:3], crop = FALSE,
                                    scale_bar = list(text_cex = 1, height = unit(.03, "native")),
                                    compass = list(style = ggspatial::north_arrow_orienteering, location = "br"), 
                                    pop_names = FALSE, ask = FALSE) +
  scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_y_continuous(expand = c(0.05, 0.05))
map <- map + theme(legend.text = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 16),
            legend.title = element_text(size = 16))


Figure_1 <- gridExtra::grid.arrange(umap_p  + ggtitle("A") + theme(plot.title = element_text(size = 20)), 
                                    map + ggtitle("B")+ theme(plot.title = element_text(size = 20)), 
                                    ngs + ylab("Ancestry Proportion") + ggtitle("C")+ theme(plot.title = element_text(size = 20)),
                                    layout_matrix = matrix(c(1, 3, 2, 3), ncol = 2))
ggsave("../../../Figure1.pdf", Figure_1, "pdf", height = 11, width = 15)

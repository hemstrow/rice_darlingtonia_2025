library(readr)

# import data
barcodes <- read_delim("RAD barcodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
plate_barcodes <- as.data.frame(read_delim("NEBNext_barcodes.txt", "\t"))
meta <- read_delim("darlingtonia_metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE, col_names = F)
meta[is.na(meta)] <- "EMPTY"
bamlist <- read.table("gbamlist.txt", stringsAsFactors = F)
plate_index <- rbind(c("1", "NEBNext23"),
                     c("2", "NEBNext22"),
                     c("3", "NEBNext21"))

# pull out barcode and plate from bamlist
bamtable <- data.frame(Plate = substr(bamlist[,1], 9, 14),
                       Index = substr(bamlist[,1], 21, 28),
                       ord = 1:nrow(bamlist))
bamtable$Plate <- plate_barcodes[match(bamtable$Plate, plate_barcodes[,2]),1]
bamtable$Plate <- plate_index[match(bamtable$Plate, plate_index[,2]),1]

# merge
colnames(meta) <- c("Plate", "DNAA", "Well", "pop", "PopSampNum")
sample.meta <- merge(meta, barcodes, by = "Well", sort = F, all.x = T)
sample.meta <- merge(sample.meta, bamtable, by = c("Plate", "Index"), sort = F)
sample.meta <- sample.meta[order(sample.meta$ord),]
sample.meta <- sample.meta[,c(5,6)]
write.table(sample.meta,"sample_metadata.txt", F, F, "\t", row.names = F, col.names = T)

# missing
good.meta <- sample.meta[sample.meta$pop != "EMPTY",]
bamlist.ord <- bamlist
bamlist.ord$ord <- 1:nrow(bamlist)
good.meta <- merge(good.meta, bamlist.ord, "ord")
write.table(good.meta$V1, "ogbamlist.txt", row.names = F, col.names = F, quote = F)

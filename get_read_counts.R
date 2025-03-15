library(readr)

#=========import read data===========
rs <- read.table("read_counts.txt", header = F, stringsAsFactors = F)
colnames(rs) <- c("reads", "bamfile")
rs$reads <- rs$reads/4

rsA <- rs[grep("_RA_", rs$bamfile),]
rsB <- rs[grep("_RB_", rs$bamfile),]
add_rs_inf <- function(rs){
  rs$Plate_seq <- substr(rs$bamfile, 9, 14)
  rs$Index <- substr(rs$bamfile, 21, 28)
  if(any(grepl("_RA_", rs$bamfile))){
    colnames(rs)[1] <- "RA_reads"
  }
  else{
    colnames(rs)[1] <- "RB_reads"
    rs$bamfile <- gsub("_RB_", "_RA_", rs$bamfile)
  }
  
  return(rs)
}
rsA <- add_rs_inf(rsA)
rsB <- add_rs_inf(rsB)
rs <- merge(rsA, rsB, by = c("bamfile", "Index", "Plate_seq", "Index")) # note, should always be the same number of reads for both!
rs$reads <-rs$RA_reads
rs$RA_reads <- NULL
rs$RB_reads <- NULL


#==========import metadata===========
barcodes <- read_delim("RAD barcodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
plate_barcodes <- as.data.frame(read_delim("NEBNext_barcodes.txt", "\t", stringsAsFactors = F))
meta <- read_delim("darlingtonia_metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE, col_names = F)
meta$X5[is.na(meta$X5)] <- "EMPTY"
bamlist <- read.table("gbamlist.txt", stringsAsFactors = F)
plate_index <- rbind(c("1", "NEBNext23"),
                     c("2", "NEBNext22"),
                     c("3", "NEBNext21"))

# pull out barcode and plate from bamlist
bamtable <- data.frame(Plate = substr(bamlist[,1], 9, 14),
                       Index = substr(bamlist[,1], 21, 28),
                       bamfile = bamlist[,1],
                       ord = 1:nrow(bamlist))
bamtable$Plate_ID <- plate_barcodes[match(bamtable$Plate, plate_barcodes[,2]),1]
bamtable$Plate_seq <- bamtable$Plate
bamtable$Plate <- plate_index[match(bamtable$Plate_ID, plate_index[,2]),1]

# merge
colnames(meta) <- c("Plate", "DNAA", "Well", "pop", "PopSampNum")
sample.meta <- merge(meta, barcodes, by = "Well", sort = F, all = T)
sample.meta <- merge(sample.meta, bamtable, by = c("Plate", "Index"), sort = F, all = T)
sample.meta <- sample.meta[order(sample.meta$ord),]
sample.meta$bamfile <- gsub(".sort.flt.bam", ".fastq", sample.meta$bamfile)

#=========merge=======

rs_full <- merge(rs, sample.meta, by = c("bamfile", "Index", "Plate_seq"))
colnames(rs_full)[1] <- "fastq"

rs_full <- rs_full[,c(1, 7, 6, 4, 8, 9)]
write.table(rs_full, "Darlingtonia_read_counts_per_sample.txt", quote = F, col.names = T, row.names = F)

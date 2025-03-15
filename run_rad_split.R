# parameters to change!

args <- commandArgs(TRUE)
read_dir <- as.character(args[1]) # directory with reads
r1 <- as.character(args[2]) # r1 file
r2 <- as.character(args[3]) # r2 file
barcodes <- as.character(args[4]) # barcodes file, can be a full pathsuch as /me/my_projects/barcodes.txt
out_prefix <- as.character(args[5]) # note, this can includes paths, such as /me/my_projects/rad_split

leading <- as.numeric(args[6]) # number of leading bases, such as "GG"
err_leading <- as.numeric(args[7]) # number of allowable errors in leading bases
err_barcode <- as.numeric(args[8]) # number of allowable errors in barcode proper

chunk_size <- as.numeric(args[9]) # Number of reads to process simultaneously during demultiplexing. Easy enough to hand this from a shell script, just use a commandArgs call. If NULL, no memory limit.



#=================set parameters, get data=============
source("match_barcodes.R")
setwd(read_dir)
r1 <- readLines(r1) 
r2 <- readLines(r2)
barcodes <- read.table(barcodes,stringsAsFactors = F)

#==================demultiplex=============
# initialize and strip down r1/r2
r1 <- process_reads(r1, barcodes, leading, chunk_size)
r2 <- process_reads(r2, barcodes, leading, chunk_size)


# determine which reads are matchable to barcodes
## candidates are missing less than err_barcode in the barcodes, are un-ambiguous between the two reads,
## are un-ambigous compared to barcodes, and must have less than err_leading in the leading in the better
## barcode misses ok
good_bar_miss_r1 <- which(r1$matches$err_barcode <= err_barcode)
good_bar_miss_r2 <- which(r2$matches$err_barcode <= err_barcode)

## leading misses ok
good_lead_miss_r1 <- which(r1$matches$err_leading <= err_leading)
good_lead_miss_r2 <- which(r2$matches$err_leading <= err_leading)

## no ambiguous barcode assignment
good_unam_r1 <- which(!is.na(r1$matches$match))
good_unam_r2 <- which(!is.na(r2$matches$match))

## good across all categories
all_goods_r1 <- intersect(intersect(good_bar_miss_r1, good_lead_miss_r1), good_unam_r1)
all_goods_r2 <- intersect(intersect(good_bar_miss_r2, good_lead_miss_r2), good_unam_r2)

## not ambiguous as to which is is RA and which is RB (only good for one of the two categories!)
not_unambig_r1 <- which(all_goods_r1 %in% all_goods_r2)
r1_reads <- all_goods_r1
if(length(not_unambig_r1) > 0){
  r1_reads <- all_goods_r1[-not_unambig_r1]
}

not_unambig_r2 <- which(all_goods_r2 %in% all_goods_r1)
r2_reads <- all_goods_r2
if(length(not_unambig_r2) > 0){
  r2_reads <- all_goods_r2[-not_unambig_r2]
}



#===============save============
interleve_reads <- function(r){
  tlines <- c(r$id, r$seq, rep("+", length(r$seq)), r$qual)
  ord <- rep(1:length(r$seq), each = 4) + (0:3) * length(r$seq)
  tlines <- tlines[ord]
}

bar_sub <- substr(barcodes[,1], (leading + 1), nchar(barcodes[1,1]))

for(i in 1:nrow(barcodes)){
  treads_r1_RA <- which(r1$matches$match[r1_reads] == bar_sub[i])
  treads_r2_RA <- which(r2$matches$match[r2_reads] == bar_sub[i])
  
  tRA <- c(interleve_reads(list(id = r1$id[r1_reads][treads_r1_RA],
                                seq = r1$seq[r1_reads][treads_r1_RA],
                                qual = r1$qual[r1_reads][treads_r1_RA])),
           interleve_reads(list(id = r2$id[r2_reads][treads_r2_RA],
                                seq = r2$seq[r2_reads][treads_r2_RA],
                                qual = r2$qual[r2_reads][treads_r2_RA])))
  
  
  tRB <- c(interleve_reads(list(id = r2$id[r1_reads][treads_r1_RA],
                                seq = r2$seq[r1_reads][treads_r1_RA],
                                qual = r2$qual[r1_reads][treads_r1_RA])),
           interleve_reads(list(id = r1$id[r2_reads][treads_r2_RA],
                                seq = r1$seq[r2_reads][treads_r2_RA],
                                qual = r1$qual[r2_reads][treads_r2_RA])))
  
  if(ncol(barcodes) == 2){
    match_to_barcode <- barcodes[i,2]
  }
  else{
    match_to_barcode <- barcodes[i,1]
  }
  writeLines(tRA, paste0(out_prefix, match_to_barcode, "_RA.fastq"), sep = "\n")
  writeLines(tRB, paste0(out_prefix, match_to_barcode, "_RB.fastq"), sep = "\n")
}
cat("Finished. Recovered ", sum(length(r1_reads), length(r2_reads)), "out of ", nrow(r1$matches), "total reads.\n")

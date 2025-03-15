library(data.table)
# read in metadata
meta <- fread("Darlingtonia_read_counts_per_sample.txt", header = T)
collection <- as.data.table(openxlsx::read.xlsx("Darlingtonia Collection Data.xlsx"))
collection <- collection[-1,]
collection[,Date := c("7/8/2017",
                      "7/8/2017",
                      "7/21/2017",
                      "7/21/2017",
                      "7/21/2017",
                      "7/21/2017",
                      "7/21/2017",
                      "8/10/2017",
                      "8/10/2017",
                      "8/10/2017",
                      "9/9/2017",
                      "9/9/2017",
                      "9/9/2017",
                      "10/7/2017")] # since dates read in wierd
collection[,Population := as.numeric(Population)]

comb <- merge(meta, collection, by.x = "pop", by.y = "Population")

#==========biosample========
# prepare biosample info
comb[,isolate := paste0("S",1:nrow(comb))]
comb[,sample_name := paste0(gsub(" ", "_", Location), "_", isolate)]
comb[,organism := "Darlingtonia Californica"]
comb[,collection_date := lubridate::mdy(gsub("/", ":", Date))]
comb[,env_broad_scale := "not collected"]
comb[,env_local_scale := "not collected"]
comb[,env_medium := "not collected"]
comb[,geo_loc_name := paste0("USA: California, ", Location)]
comb[,isol_growth_condt := "not applicable"]
comb[,lat := gsub("N.+", "", Coordinates)]
comb[,long := gsub(", ", "", gsub(" W", "", gsub(".+N", "", Coordinates)))]
arc2dec <- function(x){
  s <- strsplit(gsub(" ", "", x), "[^0-9|\\.]")
  s <- t(as.data.frame(s))
  s <- dplyr::mutate_all(as.data.frame(s), as.numeric)
  s <- s[,1] + s[,2]/60 + s[,3]/(60^2)
  return(round(s, 5))
}
comb[,lat_lon := paste0(arc2dec(lat), " N ",
                         arc2dec(long), " W")]
kcols <- c("sample_name", "organism", "isolate", 
           "collection_date", "env_broad_scale",
           "env_local_scale", "env_medium",
           "geo_loc_name", "isol_growth_condt",
           "lat_lon")
biosamp <- comb[,..kcols]
fwrite(biosamp, "biosample.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

#=======SRA metadata=========
SRA <- comb
fastqs <- list.files("D://Cluster_storage/cody/fastqs/", "^SOMM.+\\.fastq$")
RA <- sort(fastqs[grep("_RA_", fastqs)])
RB <- sort(fastqs[grep("_RB_", fastqs)])
fastqs <- lapply(list(RA, RB), function(x){
  res <- data.table(fastq = x)
  res[,plate := gsub("_R.+", "", gsub("SOMM371_", "",  fastq))]
  res[,barcode := gsub(".+_GG", "", gsub("TGCAGG.fastq", "", fastq))]
})
fastqs <- merge(fastqs[[1]], fastqs[[2]], by = c("plate", "barcode"))
SRA <- merge(fastqs, SRA, by.x = "fastq.x", by.y = "fastq")
SRA[,library_ID := paste0(plate, "_", barcode)]
SRA[,title := "RADseq of Darlingtonia Californica cuttings"]
SRA[,library_strategy := "RAD-Seq"]
SRA[,library_source := "GENOMIC"]
SRA[,library_selection := "Restriction Digest"]
SRA[,library_layout := "paired"]
SRA[,platform := "ILLUMINA"]
SRA[,instrument_model := "Illumina HiSeq 2500"]
SRA[,design_description := "RAD-seq libraries using Sbf1"]
SRA[,filetype := "fastq"]
SRA[,filename := fastq.x]
SRA[,filename2 := fastq.y]

kcols <- c("sample_name", "library_ID", "title", "library_strategy",
           "library_source", "library_selection", "library_layout",
           "platform", "instrument_model", "design_description",
           "filetype",
           "filename",
           "filename2")
fwrite(SRA[,..kcols], "SRA_meadata.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


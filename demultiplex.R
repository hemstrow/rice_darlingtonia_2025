r1 <- "test_1.fq"
r2 <- "test_2.fq"
names <- "names.txt"
mismatches <- 1
output <- "%_.fq"


lines <- R.utils::countLines(r1)
gc()

names <- readLines(names)

demultiplex <- function(r1, r2, names, mismatches, lines){
  cut <- .01
  i <- 1
  while(i <= lines){
    if(i/lines >= cut){cat(floor(i/lines*100), "% completed.\n"); cut <- cut + .01;}
    index <- scan(r1, skip = i - 1, what = "character", nlines = 1, sep = ":", quiet = T)
    index <- index[length(index)]
    dists <- adist(index, names)
    
    min.dists <- which(dists == min(dists))
    
    if(length(min.dists) == 1 & dists[min.dists][1] <= mismatches){
      t.hit <- names[min.dists]
      write(scan(r1, nlines = 4, skip = i - 1, sep = "\n", what = "character", quiet = T), gsub("%", paste0(t.hit, "_R1"), output),
            sep = "\n", append = T)
      write(scan(r2, nlines = 4, skip = i - 1, sep = "\n", what = "character", quiet = T), gsub("%", paste0(t.hit, "_R2"), output),
            sep = "\n", append = T)
    }
    i <- i + 4
  }
}

demultiplex(r1, r2, names, mismatches, lines)

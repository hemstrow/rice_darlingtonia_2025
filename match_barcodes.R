best_match <- function(gseq, lookup, leading = 2, barcode_length){

  # function to count the number of matches between sequence and barcode:
  match_count_seq <- function(seqs, bars, not.matches, save_mem = T){
      # works by repping the barcodes and comparing
      # each to the available barcodes.
      split_bars <- rep(bars, each = length(not.matches))
      split_bars <- strsplit(split_bars, "")
      split_seq <- strsplit(seqs[not.matches], "")
      
      split_bars <- unlist(split_bars)
      split_seq <- unlist(split_seq)
      
      # get the number of matches between each barcode and each sequence
      split_match <- split_seq == split_bars
      split_match <- matrix(split_match, nrow = nchar(bars[1]))
      split_match_sum <- matrixStats::colSums2(split_match)
      
      split_match_sum <- matrix(split_match_sum, nrow = length(not.matches), ncol = length(bars))
      rownames(split_match_sum) <- seqs[not.matches]
      colnames(split_match_sum) <- bars
      return(split_match_sum)
  }

  # find the exact matches
  nseq <- length(gseq)
  out <- data.frame(match = character(nseq), err_barcode = numeric(nseq), err_leading = numeric(nseq), stringsAsFactors = F)
  exact.matches <- which(gseq %in% lookup)
  if(length(exact.matches) > 0){
    out[exact.matches, 1] <- substr(lookup[match(gseq[exact.matches], lookup)], (leading + 1), nchar(gseq[1]))
    out[exact.matches, 2] <- 0
    out[exact.matches, 3] <- 0
  }

  # for not exact matches, get the count matched vs each barcode at all and leading bases.
  not.matches <- (1:nrow(out))[-exact.matches]
  if(length(not.matches) > 0){

    sub_seq <- substr(gseq, leading + 1, nchar(gseq)[1])
    sub_lookup <- substr(lookup, leading + 1, nchar(gseq)[1])
    all_bases <- match_count_seq(sub_seq, sub_lookup, not.matches)

    # get the best matching barcode for each:
    out[not.matches, 1] <- colnames(all_bases)[max.col(all_bases,ties.method="first")]
    out[not.matches, 2] <- (barcode_length - leading) - matrixStats::rowMaxs(all_bases)

    # note anywhere where multiple barcodes had equally good matches (to reject regardless)
    multi.match <- which(matrixStats::rowSums2(all_bases == out[not.matches, 2]) > 1)
    out[not.matches, 1][multi.match] <- NA

    # find counts for first two
    matching_full_barcodes <- lookup[match(out[,1][not.matches][-multi.match], sub_lookup)] # put the leading bases back on
    first_two <- unlist(strsplit(substr(gseq[not.matches][-multi.match], 1, leading), ""))
    compare <- unlist(strsplit(substr(matching_full_barcodes, 1, leading), ""))
    compare <- matrix(first_two == compare, nrow = leading, ncol = length(gseq[not.matches][-multi.match]))
    out[not.matches, 3][-multi.match] <- leading - matrixStats::colSums2(compare)
  }

  out$err_total <- out$err_barcode + out$err_leading

  return(out)
}

seperate_parts_add_meta <- function(r, barcodes, leading, barcode_length){
  barcode_length <- nchar(barcodes[1])
  r_id <- r[seq(1, length(r), by = 4)]
  r_seq <- r[seq(2, length(r), by = 4)]
  r_qual <- r[seq(4, length(r), by = 4)]
  r_barcode <- substr(r_seq, 1, barcode_length)
  r_matches <- best_match(r_barcode, barcodes, leading = leading, barcode_length)

  return(list(id = r_id, seq = r_seq, qual = r_qual, matches = r_matches))
}


process_reads <- function(r, barcodes, leading, chunk_size = 1000){
  barcode_length <- nchar(barcodes[1,1])
  
  # do this in chunks to avoid using a ton of memory if requested
  if(is.null(chunk_size)){
    r <- seperate_parts_add_meta(r, barcodes[,1], leading = leading, barcode_length = barcode_length)
  }
  else{
    bind_parts <- function(r_out){
      r_comb <- list(id = unlist(purrr::map(r_out, "id")),
                     seq = unlist(purrr::map(r_out, "seq")),
                     qual = unlist(purrr::map(r_out, "qual")),
                     matches = dplyr::bind_rows(purrr::map(r_out, "matches")))
      return(r_comb)
    }
    
    # loop through chunks
    r_chunks <- vector("list", ceiling((length(r)/4)/chunk_size))
    start <- 1
    cat("Beginning demultiplex.\n\tChunk:\n")
    for(i in 1:length(r_chunks)){
      cat("\t\t", i, "of", length(r_chunks), "\n")
      end <- ((start - 1) + chunk_size*4)
      end <- min(end, length(r))
      r_chunks[[i]] <- seperate_parts_add_meta(r[start:end], barcodes[,1], leading = leading, barcode_length = barcode_length)
      start <- end + 1
    }
    
    # bind
    r <- bind_parts(r_chunks)
    
    return(r)
  }
}

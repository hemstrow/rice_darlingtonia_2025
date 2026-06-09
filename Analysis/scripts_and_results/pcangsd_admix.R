args <- commandArgs(TRUE)
infile <- as.character(args[1])
outfile <- as.character(args[2])
max_K <- as.character(args[3])
par <- as.numeric(args[4])

cat("Beginning run. Parms:\n\tInfile:", infile, "\n\tOutfile prefix:", outfile, "\n\tmax K:", max_K)

cmd <- paste0("pcangsd", 
              " -b ", infile,
              " -t ", par, 
              " -o ", outfile,
              " --inbreed-samples")
system(cmd)

for(i in 2:max_K){
  cat("Running K = ", i, "\n")
  cmd <- paste0("pcangsd", 
                " -b ", infile,
                " -e ", i - 1,
                " -t ", par, 
                " -o ", paste0(outfile, "K_", i),
                " --admix")
  system(cmd)
}

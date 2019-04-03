
library(dada2); packageVersion("dada2")
# Filename parsing
path <- "/Users/Clara/Projects/osd/data/miseq/ERR771024" # directory containing your demultiplexed fastq files
filtpath <- "/Users/Clara/Projects/osd/res/dada2/miseq/ERR771024" # Filtered files go into the filtered/ subdirectory
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              truncLen=300, maxEE=1, truncQ=11, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
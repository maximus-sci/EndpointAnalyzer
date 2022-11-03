#set data location and analysis options
dataDir    <- "I:/Bed"
outDir     <- NULL   #will default to ./dataDir/results
refSeq     <- "I:/Cx3cr1_11.fa"
chromSizes <- "I:/Cx3cr1_11.chrom.sizes"


  
runOpts <- list(quality  = 2,            #reads with quality scores below this number will be removed. 1 or 0 are multimappers
                qualOnly = T,            #if T, only the quality filter is used. If F, region and size filters are also used
                exclude1 = c(18,1077),   #reads starting and ending within this region will be removed
                exclude2 = c(1662,2720), #reads starting and ending within this region will be removed
                minReadSize = 500,       #reads shorter than this will be excluded. To remove filter, set it to 0
                maxReadSize = 2733)      #reads longer than this will be excluded. To remove filter, set it to Inf




{
#load packages and set up parallel processing
  install.packages(c("dplyr","stringr","readr","seqinr","doParallel","foreach"))
  
  require(dplyr)
  require(stringr)
  require(readr)
  require(seqinr)
  require(doParallel)
  require(foreach)
  registerDoParallel(cores=detectCores() - 1)
  
  
#setup functions and directory structure
  if(file.exists(paste0(outDir,"/Convert_NormBedgraph_to_bigwig.sh"))){stop("You have already run this script in this directory. Remove previous results to continue")}
  if(length(outDir)==0){outDir <- paste0(dataDir,"/results")}
  dir.create(outDir)
  refSeq <- seqinr::read.fasta(file = refSeq, as.string = T)
  
  checkSequence <- function(bedData){
    bedData %>%
      as.data.frame()
    if (length(unique(bedData$seqnames)) > 1) {stop("Input bed file has reads from multiple reference sequences")}
    if (unique(bedData$seqnames) != getName(refSeq)) {stop(paste0("Input bed file does not match",getName(refSeq)))}
    if (max(bedData$end) > getLength(refSeq)) {stop("Annotation error, bed file has reads outside bounds of reference sequence")}
  }
  
  dir.create(paste0(outDir,"/filtered_bed"))
  dir.create(paste0(outDir,"/endpoint_bedgraph"))
  dir.create(paste0(outDir,"/bigWig"))


#gather data 
  inBeds <- list.files(dataDir,pattern="*.bed",ignore.case = T,full.names = T)
  if(length(inBeds)==0){stop(paste0("No input bed files found in: ",dataDir))
    } else {message(paste0(length(inBeds)," bed files located in input directory"))}


#filter bed files and dump to files
  bedDat <- list()
  maxLen <- c()
  for (i in 1:length(inBeds)){
    inFile     <- rtracklayer::import(inBeds[i], format="bed") %>% as.data.frame()
    checkSequence(inFile)
    inSize     <- nrow(inFile)
    maxLen     <- append(maxLen,nchar(inBeds[i]))
    bedDat[[i]]<- inFile %>%
      dplyr::filter(score >= runOpts$quality)
    if(!runOpts$qualOnly){
      bedDat[[i]]  <- bedDat[[i]] %>%
        dplyr::filter(!(start >= runOpts$exclude1[1] & end <= runOpts$exclude1[2])) %>%
        dplyr::filter(!(start >= runOpts$exclude2[1] & end <= runOpts$exclude2[2])) %>%
        dplyr::filter(width > runOpts$minReadSize) %>%
        dplyr::filter(width < runOpts$maxReadSize)}
    names(bedDat)[i] <- gsub(".bed","",gsub("^.*/", "", inBeds[i]))
    if (nchar(inBeds[i]) < max(maxLen)){
      gap <- paste0(rep(" ",max(maxLen)-nchar(inBeds[i])), collapse = "")
      } else {gap <- ""}
    message(paste0("For ",inBeds[i]," kept",gap,"\t",nrow(bedDat[[i]])," out of ",inSize," lines (",round((nrow(bedDat[[i]])/inSize)*100,digits = 1),"%)"))
  }
  
  filtBeds <- list()
  for (i in 1:length(bedDat)){
    write_tsv(bedDat[[i]], file = paste0(outDir,"/filtered_bed/BothStrands_",names(bedDat)[i],".filtered.bed"))
    filtBeds[[paste0("Both_",names(bedDat)[i])]] <- bedDat[[i]]
    write_tsv(bedDat[[i]] %>%
                dplyr::filter(strand == "+"),
              file = paste0(outDir,"/filtered_bed/PlusStrand_",names(bedDat)[i],".filtered.bed"))
    filtBeds[[paste0("Plus_",names(bedDat)[i])]] <- bedDat[[i]] %>%
      dplyr::filter(strand == "+")
    write_tsv(bedDat[[i]] %>%
                dplyr::filter(strand == "-"),
              file = paste0(outDir,"/filtered_bed/MinusStrand_",names(bedDat)[i],".filtered.bed"))
    filtBeds[[paste0("Minus_",names(bedDat)[i])]] <- bedDat[[i]] %>%
      dplyr::filter(strand == "-")
  }
  message("\n************************\nFiles filtered, beginning endpoint analysis\n")


#perform endpoint analysis
  bgDat <- list()
  bgDat <- foreach (i = 1:length(filtBeds), .packages = c("tibble","seqinr","hash")) %dopar%{
    bg <- tibble(seqnames = getName(refSeq), start = seq(1:(getLength(refSeq)))-1, end = seq(1:(getLength(refSeq))), depth = 0)
    h <- hash(keys = seq_len(getLength(refSeq)), values = 0)
    for (j in seq_len(nrow(filtBeds[[i]]))){
      h[[as.character(filtBeds[[i]][j,]$start)]] <- h[[as.character(filtBeds[[i]][j,]$start)]] + 1
      h[[as.character(filtBeds[[i]][j,]$end)]]   <- h[[as.character(filtBeds[[i]][j,]$end)]] + 1
    }
    for (k in seq_len(getLength(refSeq))){
      bg[k,4] <- h[[as.character(k)]]
    }
    bg
    }

  bgDat<-setNames(bgDat, names(filtBeds))
  message("Endpoint analysis complete, normalizing bedGraph files")
  
  #normalize bedGraph files
  normBg <- list()
  for (i in seq_len(length(bgDat))){
    totalReads <- bgDat[[i]] %>%
      summarise(sum(depth)) %>%
      pull(`sum(depth)`)
    normBg[[names(bgDat)[i]]] <- bgDat[[i]] %>%
      mutate(depth = (depth+0.01)/totalReads)
  }
  
  message("Bedgraphs normalized, writing data to disk")
  
  invisible(lapply(names(normBg), function(x){
    write_tsv(normBg[[x]], file = paste0(outDir,"/endpoint_bedgraph/",x,".norm.filt.bg"))
  }))
  
  #if (file.exists(paste0(outDir,"/Convert_NormBedgraph_to_bigwig.sh"))){file.remove(paste0(outDir,"/Convert_NormBedgraph_to_bigwig.sh"))}
  for (i in 1:length(names(normBg))){
    write(paste0("home/meganfr/bedGraphToBigWig ",paste0(outDir,"/endpoint_bedgraph/",names(normBg)[i],".norm.filt.bg "),chromSizes," ",
                 paste0(outDir,"/bigWig/",names(normBg)[i],".norm.filt.bw")),
          file = paste0(outDir,"/Convert_NormBedgraph_to_bigwig.sh"), append = T)
  }
  
  message("To get bigWigs, run:")
  message(paste0(". ", outDir,"/Convert_NormBedgraph_to_bigwig.sh"))


#write the log
  
  if (file.exists(paste0(outDir,"/",Sys.Date(),"_logFile.txt"))){file.remove(paste0(outDir,"/",Sys.Date(),"_logFile.txt"))}
  
  write(paste0("Date of run: ",Sys.time()),
        file = paste0(outDir,"/",Sys.Date(),"_logFile.txt"), append = T)
  write(paste0("DataInput: ",dataDir),
        file = paste0(outDir,"/",Sys.Date(),"_logFile.txt"), append = T)
  write(paste0("Output: ",outDir,"\r\rOptions used:"),
        file = paste0(outDir,"/",Sys.Date(),"_logFile.txt"), append = T)
  for (i in seq_len(length(runOpts))){
    write(paste0("\t",names(runOpts)[i]," : ",runOpts[[i]]),
          file = paste0(outDir,"/",Sys.Date(),"_logFile.txt"), append = T)}
  }

# normR fine enrichment
#
# helmuth 2017-05-30
#
runMean <- function(x, n=5) {
  cx <- cumsum(as.numeric(x))
  return((cx[n:length(cx)] - c(0, cx[1:(length(cx)-n)]))/n)
}
binMean <- function(x, n=5) {
  rmean <- runMean(x, n)
  bmean <- rmean[seq(1,length(rmean),n)]
  xl <- length(x)
  if (xl%%n==0) { #x is even in n
    return(bmean)
  } else { #x is odd in n -> add mean of rest of list
    return(c(bmean, mean(x[(xl-xl%%n+1):xl])))
  }
}
getEnrichmentFine <- function(chip.bampath, ctrl.bampath, enrichR.fit,
                              bs.fine = 25, scale.pseudo = T, 
                              method = c("counts", "coverage"),
                              bgIdx = 1, fgIdx = 2, procs = 1, ...) {
  method <- match.arg(method)

  #Get fit information
  require(normr)
  gr.init <- getRanges(enrichR.fit)
  cts.init <- getCounts(enrichR.fit)
  postB.init <- getPosteriors(enrichR.fit)[,bgIdx]
  
  #Count in new binsize
  require(bamsignals)
  gr <- GRanges(seqnames(seqinfo(gr.init)), IRanges(1, seqlengths(seqinfo(gr.init))))
  gr <- gr[seqnames(gr) == 'chrX']
  if (method == "counts") { #1bp counting
    cts.fine <- mclapply(c(ctrl.bampath, chip.bampath), bamProfile, 
                         gr=gr, binsize=bs.fine, ss=F, mc.cores = procs, ...)
    cts.fine <- lapply(cts.fine, function(x) unlist(as.list(x@signals)))
  } else { #count fragments overlapping bins
    cts.fine <- mclapply(c(ctrl.bampath, chip.bampath), function(f) {
      prf <- bamCoverage(f, gr, ...)
      prf.bins <- binMean(unlist(as.list(prf)), bs.fine)
      return(prf.bins)
    }, mc.cores = procs)
  }
  
  #Get model-derived pseudo counts and scale if required
  ps <- mclapply(cts.init, function(x) {
          sum(postB.init * x)/sum(postB.init)
        }, 
        mc.cores = procs)
  if (scale.pseudo) {
    ps <- lapply(ps, "/", width(gr.init)[1]/bs.fine)
  }
  
  #Calculate fine enrichR enrichment
  fc <- log((cts.fine[[2]] + ps[[2]])/(cts.fine[[1]] + ps[[1]]))
  regul <- log(ps[[1]]/ps[[2]])
  stdrz <- log(enrichR.fit@theta[fgIdx]/(1-enrichR.fit@theta[fgIdx]) * 
               (1-enrichR.fit@theta[bgIdx])/enrichR.fit@theta[bgIdx])
  e <- (fc + regul)/stdrz
  e[e < 0] = 0; e[e > 1] = 1; e[is.na(e) | is.infinite(e)] = 0
  e <- e*100 #scale for histone density

  #Construct a genomic ranges with scores as enrichment
  require(GenomicRanges)
  gr.fine <- tileGenome(seqinfo(gr.init), tilewidth=bs.fine,cut.last.tile.in.chrom=T)
  gr.fine = gr.fine[seqnames(gr.fine) == 'chrX']
  gr.fine$score <- e
 
  return(gr.fine)
}


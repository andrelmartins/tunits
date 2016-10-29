#
# TUnit HMM Decoding (generic version)
#

blocks.to.coords <- function(blocks, step) {
  apply(blocks, 2, function(col) {
    # col = [start, end]
    c((col[1] - 1) * step, col[2] * step) # [bed start, bed end)
  })
}

coords.to.bed6 <- function(chrom, strand, coords) {
  return(data.frame(chrom = chrom, start = coords[1,], end = coords[2,], name = "", score = 0, strand = strand))
}

decode.simple <- function(chrom, strand, step, path, states) {
  blocks = path.blocks.qhmm(path, states)
  # 2 x N matrix (rows: start, end; 1-based closed intervals)

  # map back if minus strand
  if (strand == '-') {
    L = length(path)
    starts = rev(L - blocks[2,] + 1)
    ends = rev(L - blocks[1,] + 1)
    blocks = rbind(starts, ends)
  }

  coords = blocks.to.coords(blocks, step) 
  # 2 x N matrix
  
  return(coords.to.bed6(chrom, strand, coords))
}

decode.transform <- function(chrom, strand, step, path, states, positions) {
  blocks = path.blocks.qhmm(path, states)
  # 2 x N matrix (rows: start, end; 1-based closed intervals)
  
  # use positions to remap blocks back into un-transformed space
  blocks.remapped = rbind(positions[blocks[1,]], positions[blocks[2,]]) # TODO: consider adding a +1 here, since we accepted this distance as true ...
  # 2 x N matrix (rows: start, end; 1-based closed intervals)
  
  # for minus strand, need to reverse ranges
  if (strand == '-')
    blocks.remapped = rbind(rev(blocks.remapped[2,]), rev(blocks.remapped[1,]))
  
  # get actual coordinates
  coords = blocks.to.coords(blocks.remapped, step) 
  # 2 x N matrix

  return(coords.to.bed6(chrom, strand, coords))
  #return(coords.to.bed6(chrom, strand, blocks.remapped))
}

#' Decode dataset using given HMM
#'
#' @param hmm QHMM instance to use for viterbi decoding
#' @param dataset dataset produced by load.dataset or load.data
#' @param states vector of state numbers that count as transcribed
#' @param covar.lst list of covariate tables for each sequence in dataset
#' @param min.len if not NA, then predictions smaller than min.len bp are culled
#' @return bed6 data.frame with predictions
#' @export
decode.dataset <- function(hmm, dataset, states, covar.lst = NULL, min.len = NA) {
  bed = NULL
  N = length(dataset)
  for (i in 1:N) {
    chrom = dataset[[i]]$chrom
    strand = dataset[[i]]$strand
    step = dataset[[i]]$step
    
    path = viterbi.qhmm(hmm, dataset[[i]]$data, covars = covar.lst[[i]], missing = dataset[[i]]$missing)
    
    if (is.null(dataset[[i]]$positions))
      bed = rbind(bed, decode.simple(chrom, strand, step , path, states))
    else
      bed = rbind(bed, decode.transform(chrom, strand, step, path, states, dataset[[i]]$positions))
  }
  
  # drop small predictions
  if (!is.na(min.len)) {
    sizes = bed[,3] - bed[,2]
    bed = bed[sizes >= min.len,]
  }
  
  # name predictions
  pnames = vector(mode="character", length=dim(bed)[1])
  idx.plus = bed[,6] == '+'
  pnames[idx.plus] = paste("preds_plus_", 1:sum(idx.plus), sep='')
  idx.minus = !idx.plus
  pnames[idx.minus] = paste("preds_minus_", 1:sum(idx.minus), sep='')

  bed[,4] = pnames
  
  return(bed)
}

#' Write BED data.frame to disk
#'
#' @param bed bed data.frame
#' @param file filename or connection
#' @export
write.bed <- function(bed, file) {
  # make sure coordinates are integer values
  starts = as.integer(bed[,2])
  ends = as.integer(bed[,3])
  bed[,2] = starts
  bed[,3] = ends

  # write to disk
  write.table(bed, file=file, col.names=F, row.names=F, quote=F, sep='\t')
}

#' Write BED6 data.frame to disk as a UCSC browser track
#'
#' @param bed6 BED6 data.frame
#' @param name track name
#' @param description track description
#' @param filename track filename
#' @param extra extra track header arguments
#' @param with.color if TRUE, columns for BED interval color are added
#' (+ as red, - as blue) and itemRgb=On is added to track flags.
#' @export
write.track <- function(bed6, name, description, filename, extra="", with.color = TRUE) {
  # make sure coordinates are integer values
  starts = as.integer(bed6[,2])
  ends = as.integer(bed6[,3])
  bed6[,2] = starts
  bed6[,3] = ends

  # handle colors
  if (with.color) {
    extra = paste(extra, "itemRgb=On")
    # extend bed6 to include colors
    colors = rep("197,0,11", dim(bed6)[1])
    colors[bed6[,6] == '-'] = "0,132,209"
    
    bed6 = cbind(bed6, bed6[,2:3], colors)
  }

  # write to file
  fout = file(filename, "w")
  cat("track name=", name, " description=\"", description, "\" ",extra, "\n", sep='', file=fout)
  write.table(bed6, file=fout, col.names=F, row.names=F, quote=F, sep='\t')
  close(fout)
}

#' Write BED6 data.frame to disk as a UCSC browser track
#'
#' The two input files 'bed6.long' and 'bed6.short' should correspond to the same set of predictions, but with different extensions. These are then used to produced a "thickStart" and "thickEnd" columns in accordance with strand information.
#'
#' @param bed6.long BED6 data.frame (extended version)
#' @param bed6.short BED6 data.frame (compact version)
#' @param name track name
#' @param description track description
#' @param filename track filename
#' @param extra extra track header arguments
#' @param with.color if TRUE, columns for BED interval color are added
#' (+ as red, - as blue) and itemRgb=On is added to track flags.
#' @export
write.extended.track <- function(bed6.long, bed6.short, name, description, filename, extra="", with.color = TRUE) {
  bed6 = bed6.long
  
  # make sure coordinates are integer values
  starts = as.integer(bed6[,2])
  ends = as.integer(bed6[,3])
  bed6[,2] = starts
  bed6[,3] = ends
  
  # extend BED6 with thickStart,thickEnd
  tstart = bed6[,2]
  tend = bed6[,3]
  
  chrom.bed6.short = split.by.chrom.bed(bed6.short, cols = 1:6)
  
  foreach.bed(bed6, function(i, chrom, start, end, strand) {
    bl.i = chrom.bed6.short[[chrom]]
    
    if (strand == '+') {
      idx = which(bl.i[,2] == start & bl.i[,6] == strand)
      tend[i] <<- bl.i[idx, 3]
    } else {
      idx = which(bl.i[,3] == end & bl.i[,6] == strand)
      tstart[i] <<- bl.i[idx, 2]
    }
  })
  
  bed6 = data.frame(bed6, as.integer(tstart), as.integer(tend))

  # handle colors
  if (with.color) {
    extra = paste(extra, "itemRgb=On")
    # extend bed6 to include colors
    colors = rep("197,0,11", dim(bed6)[1])
    colors[bed6[,6] == '-'] = "0,132,209"
    
    bed6 = data.frame(bed6, colors)
  }

  # write to file
  fout = file(filename, "w")
  cat("track name=", name, " description=\"", description, "\" ",extra, "\n", sep='', file=fout)
  write.table(bed6, file=fout, col.names=F, row.names=F, quote=F, sep='\t')
  close(fout)
}

#' Write BED6 data.frame to disk as a UCSC browser track
#'
#' The two input files 'bed6.long' and 'bed6.short' should correspond to the same set of predictions, but with different extensions. These are then used to produced a "thickStart" and "thickEnd" columns in accordance with strand information.
#'
#' @param bed6.long BED6 data.frame (extended version)
#' @param bed6.inner BED6 data.frame (compact version)
#' @param name track name
#' @param description track description
#' @param filename track filename
#' @param extra extra track header arguments
#' @param with.color if TRUE, columns for BED interval color are added
#' (+ as red, - as blue) and itemRgb=On is added to track flags.
#' @export
write.extended2.track <- function(bed6.long, bed6.inner, name, description, filename, extra="", with.color = TRUE) {
  bed6 = bed6.long
  
  # make sure coordinates are integer values
  starts = as.integer(bed6[,2])
  ends = as.integer(bed6[,3])
  bed6[,2] = starts
  bed6[,3] = ends
  
  # extend BED6 with thickStart,thickEnd
  tstart = bed6[,2]
  tend = bed6[,3]
  
  chrom.bed6.inner = split.by.chrom.bed(bed6.inner, cols = 1:6)
  
  foreach.bed(bed6, function(i, chrom, start, end, strand) {
    bl.i = chrom.bed6.inner[[chrom]]
    
    idx = which(bl.i[,2] >= start & bl.i[,3] <= end & bl.i[,6] == strand)
    tstart[i] <<- max(start, bl.i[idx, 2])
    tend[i] <<- min(end, bl.i[idx, 3])
  })
  
  bed6 = data.frame(bed6, as.integer(tstart), as.integer(tend))

  # handle colors
  if (with.color) {
    extra = paste(extra, "itemRgb=On")
    # extend bed6 to include colors
    colors = rep("197,0,11", dim(bed6)[1])
    colors[bed6[,6] == '-'] = "0,132,209"
    
    bed6 = data.frame(bed6, colors)
  }

  # write to file
  fout = file(filename, "w")
  cat("track name=", name, " description=\"", description, "\" ",extra, "\n", sep='', file=fout)
  write.table(bed6, file=fout, col.names=F, row.names=F, quote=F, sep='\t')
  close(fout)
}

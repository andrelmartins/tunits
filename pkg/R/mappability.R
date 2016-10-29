#
# collect mappability info in a given chrom or range
#

process.map.values <- function(map.values, strand, step, read.len, read.left.edge, threshold.fraction) {
  do.shift = (strand == '-' && read.left.edge) || (strand == '+' && !read.left.edge)
  
  if (do.shift)
    map.values = c(rep(0, read.len - 1), map.values[1:(length(map.values) - read.len + 1)])

  # steppify
  if (step > 1) {
    N = length(map.values) %/% step
    thresh = threshold.fraction * step
    sapply(1:N, function(idx) {
      sstart = ((idx - 1) * step + 1)
      send = sstart + step - 1
      
      count = sum(map.values[sstart:send], na.rm=T)
      if (count >= thresh)
        1
      else
        0
    })
  } else
  map.values
}

#' Compute an vector of unmappable positions (=1)
#'
#' @param bwMap mappability bigWig object (0 if mappable, 1 if unmappable)
#' @param chrom chromosome name
#' @param strand target strand
#' @param step step size in bp
#' @param read.len length of reads used in mapping (must match the read length used to compute the mappability bigWig object)
#' @param read.left.edge logical, use TRUE for GRO-seq style, FALSE for PRO-seq style
#' @param threshold.fraction fraction of reads that must be unmappable for a step to be considered unmappable
#' @return integer vector with 0 on mappable positions and 1 in unmappable positions
#' @export
chrom.unmappable <- function(bwMap, chrom, strand, step, read.len, read.left.edge, threshold.fraction) {
  map.values = chromStepSum.bigWig(bwMap, chrom, defaultValue=0, step=1)

  return(process.map.values(map.values, strand, step, read.len, read.left.edge, threshold.fraction))
}

#' Compute an vector of unmappable positions (=1)
#'
#' @param bwMap mappability bigWig object (0 if mappable, 1 if unmappable)
#' @param chrom chromosome name
#' @param start start coordiante of region (0-based)
#' @param end end coordinate of region (1-based)
#' @param strand target strand
#' @param step step size in bp
#' @param read.len length of reads used in mapping (must match the read length used to compute the mappability bigWig object)
#' @param read.left.edge logical, use TRUE for GRO-seq style, FALSE for PRO-seq style
#' @param threshold.fraction fraction of reads that must be unmappable for a step to be considered unmappable
#' @return integer vector with 0 on mappable positions and 1 in unmappable positions
#' @export
range.unmappable <- function(bwMap, chrom, start, end, strand, step, read.len, read.left.edge, threshold.fraction) {
  map.values = queryByStep.bigWig(bwMap, chrom, start, end, 1, do.sum = F, default.null = F, defaultValue = 0)

  return(process.map.values(map.values, strand, step, read.len, read.left.edge, threshold.fraction))
}

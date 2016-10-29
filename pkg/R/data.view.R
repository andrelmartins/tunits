#
# Data view : functions to visualize data transformation
#

#' Apply the data transformation function to a given region
#' @export
region.transform <- function(bw, chrom, start, end, strand, step = 50, bwmap = NULL) {
  # collect reads
  reads = step.bpQuery.bigWig(bw, chrom, start, end, step, abs.value = TRUE)
  
  tdata = NULL
  if (is.null(bwmap)) {
    tdata = transform_data(reads, reverse = (strand == '-'))
  } else {
    unmap = step.bpQuery.bwMap(bwmap, chrom, start, end, step, strand)
    tdata = transform_data(reads, reverse = (strand == '-'), unmap)
  }
  
  return(c(
      list(chrom = chrom, strand = strand, step = step),
      tdata,
      list(missing = NULL)))
}

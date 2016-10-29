#
# Functions that prepare data for TU HMMs
#

# Details about data transformation:
#
# - per chromosome structure
# - step size
# - auxiliary position vector (in case of transformative HMM)
# - strand
# - chrom
# - optional mappability track

transform_data <- function(reads, reverse = FALSE, unmap = NULL) {
  positions = which(reads > 0)
  n = length(positions)
  
  if (n == 0) return(NULL)
  if (!reverse) {
    distances = positions[2:n] - positions[1:(n-1)]

    final.positions = positions[1:(n-1)]

    if (!is.null(unmap)) {
      unmap.dist = sapply(2:n, function(idx) {
      x1 = positions[idx - 1] + 1
      x2 = positions[idx] - 1
      
      if (x2 >= x1) {
        sum(unmap[x1:x2])
      } else
        0
      })
      
      X = rbind(distances, reads[final.positions])
      row.names(X) <- c("distance", "reads")
      return(list(positions = final.positions, unmap.distance = unmap.dist, data = X))
    } else {
      X = rbind(distances, reads[final.positions])
      row.names(X) <- c("distance", "reads")
      return(list(positions = final.positions, data = X))
    }
  } else {
    distances = positions[2:n] - positions[1:(n-1)]

    final.positions = positions[2:n]
    
    if (!is.null(unmap)) {
      unmap.dist = sapply(2:n, function(idx) {
      x1 = positions[idx - 1] + 1
      x2 = positions[idx] - 1
      
      if (x2 >= x1) {
        sum(unmap[x1:x2])
      } else
        0
      })

      X = rbind(rev(distances), rev(reads[final.positions]))
      row.names(X) <- c("distance", "reads")
      return(list(positions = rev(final.positions), unmap.distance = rev(unmap.dist), data = X))
      
    } else {
      X = rbind(rev(distances), rev(reads[final.positions]))
      row.names(X) <- c("distance", "reads")
      return(list(positions = rev(final.positions), data = X))
    }
  }
}

load_data_sequence <- function(chrom, bw, step, transform, bwMap, map.left.edge, map.read.len, strand, offset = 0) {
  threshold.fraction = 0.5 # for mappability

  # collect reads
  reads = step.bpQuery.bigWig(bw, chrom, NULL, NULL, step, abs.value = TRUE) + offset

  # process reads / mappability
  if (transform) {
    map.info = NULL
    
    # for transformed data, mappability comes in the form of an additional
    # dataset item, 'unmap.distances', produced by the 'transform' function
    # that accounts for how much of the distance between consecutive reads
    # is _potentially_ due to unmappable positions
    #
    unmap = NULL
    if (!is.null(bwMap)) {
      unmap = step.bpQuery.bwMap(bwMap, chrom, NULL, NULL, step, strand)
    }
    
    return(c(
      list(chrom = chrom, strand = strand, step = step),
      transform_data(reads, reverse = (strand == '-'), unmap),
      list(missing = NULL)))
  } else {
    map.info = NULL
    if (!is.null(bwMap))
      #map.info = chrom.unmappable(bwMap, chrom, strand, step, map.read.len, map.left.edge, threshold.fraction)
      map.info = step.bpQuery.bwMap(bwMap, chrom, NULL, NULL, step, strand)
      
    if (strand == '+')
      return(list(chrom = chrom, strand = strand, step = step, positions = NULL, data = reads, missing = map.info))
    else
      return(list(chrom = chrom, strand = strand, step = step, positions = NULL, data = rev(reads), missing = rev(map.info)))
  }
}

#' Load data for TUnit HMMs
#'
#' @param chroms character vector with chromosome names
#' @param bwPlus bigWig object for plus strand
#' @param bwMinus bigWig object for minus strand
#' @param step step size in base pairs
#' @param transform apply non-zero reads / distance transformation
#' @param bwMap bigWig object with mappability information (0 for mappable, 1 for unmappable)
#' @param map.left.edge TRUE for GRO-seq files, FALSE for PRO-seq files
#' @param map.read.len length of reads in source bigWig files (and bwMap input)
#' @param offset integer value to sum to read counts
#' @return list with a per chrom/strand pair list containing the fields: chrom, step, positions, data and map.info
#' @export
load.data <- function(chroms, bwPlus, bwMinus, step, transform = FALSE, bwMap = NULL, map.left.edge = TRUE, map.read.len = NA, offset = 0) {
  nChroms = length(chroms)
  result = vector(mode="list", length=nChroms*2)
  
  for (i in 1:nChroms) {
    j1 = (i - 1) * 2 + 1
    j2 = j1 + 1
    
    chrom = chroms[i]
    cat(" *", chrom, "+\n")
    result[[j1]] = load_data_sequence(chrom, bwPlus, step, transform, bwMap, map.left.edge, map.read.len, "+", offset)
    
    cat(" *", chrom, "-\n")
    result[[j2]] = load_data_sequence(chrom, bwMinus, step, transform, bwMap, map.left.edge, map.read.len, "-", offset)
  }
  
  return(result)
}

#' Load data for TUnit HMMs
#'
#' @param bwDataSet a data frame with the following components:
#' \describe{
#' \item{read.len}{integer holding the read length of the PRO/GRO-seq dataset}
#' \item{left.edge}{boolean indicating if the left edge of the read is the one being mapped}
#' \item{bw.plus}{path to the bigWig file containing the pile-up for the plus strand}
#' \item{bw.minus}{path to the bigWig file containing the pile-up for the minus strand}
#' \item{bw.map}{path to the bigWig file containing the mappability information}
#' }
#' @param offset integer value to sum to read counts
#' @export
load.dataset <- function(chroms, bwDataSet, step, transform = FALSE, with.map = FALSE, offset = 0) {
  # load bigWig files
  bwPlus = load.bigWig(bwDataSet$bw.plus)
  bwMinus = load.bigWig(bwDataSet$bw.minus)

  # get data
  result = NULL
  if (with.map) {
    threshold.fraction = 0.5 # for mappability
    
    #bwMap = load.bigWig(bwDataSet$bw.map)
    bwMap = load.bwMap(bwDataSet$bw.map, bwDataSet$read.len, bwDataSet$left.edge, threshold.fraction)
    result = load.data(chroms, bwPlus, bwMinus, step, transform, bwMap, bwDataSet$left.edge, bwDataSet$read.len, offset = offset)
    unload.bwMap(bwMap)
  } else
    result = load.data(chroms, bwPlus, bwMinus, step, transform, offset = offset)
  
  # clean up
  unload.bigWig(bwPlus)
  unload.bigWig(bwMinus)
  
  #
  return(result)
}

#' View of a dataset suitable for em.qhmm use
#' @export
train.dataset <- function(dataset) {
  lapply(dataset, function(d) d$data)
}

#' Scan mappability bigWig to collect totals
#'
#' @return pair: total mappable bases, total bases
#' @export
mapSummary <- function(bwPath) {
  bwMap = load.bigWig(bwPath)
  
  totalMap = 0
  totalSize = 0
  for (chrom in bwMap$chroms) {
    tmp = chromStepSum.bigWig(bwMap, chrom, 1, 0)
    totalSize = totalSize + length(tmp)
    totalMap = totalMap + (length(tmp) - sum(tmp))
  }
  
  return(c(totalMap, totalSize))
}


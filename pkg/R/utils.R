#
# utils
#

#' Split BED set by chromosome
#'
#' @param bed data.frame with BED content
#' @param chroms character vector with optional set of chroms to consider (otherwise all chroms in the BED are used)
#' @param cols integer vector with optional set of columns to return (otherwise all columns are returned)
#' @param with.backidx logical value, if TRUE, an extra column with the indexes of each entry is added to the result
#' @return list of BED data.frames, one per chrom
#' @export split.by.chrom.bed
split.by.chrom.bed <- function(bed, chroms = levels(bed[,1]), cols = NULL, with.backidx = FALSE) {
  if (is.null(cols))
    cols = 1:(dim(bed)[2])

  result = NULL
  if (with.backidx)
      result = lapply(chroms, function(chrom) {
        idxs = which(bed[,1] == chrom)
        cbind(bed[idxs, cols], idxs)
      })
  else
      result = lapply(chroms, function(chrom) bed[bed[,1] == chrom, cols])
  names(result) <- chroms

  return(result)
}


#' Quick plot of HMM using Rgraphviz
#'
#' See http://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html
#'
#' @param hmm QHMM object instance
#' @param state.names character vector of state names
#' @param with.probs logical, if TRUE, add transition probabilities to edges
#' @export
plot.qhmm <- function(hmm, state.names = as.character(seq(1:hmm$n.states)), with.probs = FALSE) {
  library(Rgraphviz)
  
  edge.to.prob <- function(edge) { # This does not seem to be working all that well
    params = collect.params.qhmm(hmm)
    parts = unlist(strsplit(edge, "~"))
    
    src = which(state.names == parts[1])
    tgt = which(state.names == parts[2])

    src.outdegree = sum(hmm$valid.transitions[src,] > 0)
    src.pars = params$transitions[[src]]
    
    if (length(src.pars) == 1) { # assume auto-corr
      if (src == tgt)
        as.character(src.pars)
      else
        as.character(round((1 - src.pars) / (src.outdegree - 1), 2))
    } else if (length(src.pars) == src.outdegree) {
      idx = hmm$valid.transitions[src, tgt]
      as.character(src.pars[idx])
    } else ""
  }
  
  aMat = hmm$valid.transitions
  colnames(aMat) <- state.names
  rownames(aMat) <- state.names
  gr = new("graphAM", adjMat = aMat, edgemode="directed")
  
  if (with.probs) {
    pars = collect.params.qhmm(hmm)
    
    labels = edgeNames(gr)
    names(labels) <- sapply(labels, edge.to.prob)
    gr <- layoutGraph(gr, edgeAttrs=list(label=labels))
    print(names(labels))
  }
  plot(gr, attrs = list(edge = list(arrowsize = 0.5)))
  #renderGraph(gr)
}

#' Produce a flattened vector of non-fixed HMM parameters
#'
#' @param params parameter object from 'collect.params.qhmm'
#' @param name parameter name prefix
#' @export params.flatten
params.flatten <- function(params, name = "") {
    if (is.list(params)) {
       res = NULL
       
       # build up names
       n = names(params)
       if (is.null(n))
         n = as.character(1:length(params))
       else
         n = substr(n, 1, 1)
        
       for (i in 1:length(params)) {
         part = params[[i]]
         ni = paste(name, n[i], sep='.')
         res = c(res, params.flatten(part, ni))
       }
       return(res)
    } else if (is.null(params)) {
       return(NULL)
    } else {
       aux = attr(params, "fixed")
       if (!is.null(aux)) {
         if (sum(!aux) > 0) {
            res = params[!aux]
            ni = name
            if (sum(!aux) > 1) {
                ni = paste(name, which(!aux), sep='.')
            }
            names(res) <- ni
            return(res)
         }
       }
       return(NULL)
    }
}

#' Collect parameters from all chromosomes
#' 
#' Assuming that parameters are in files '<chrom>.params.Rdata' for
#' each chromosome, collect all of them into a table
#' @export
chromparams.table <- function() {
  files = list.files(pattern=glob2rx("*.params.Rdata"))
  chroms = sapply(strsplit(files, "[.]"), function(parts) parts[1])

  res = sapply(files, function(filename) {
    load(filename) # defines hmm.params
    params.flatten(hmm.params)
  })
  
  res = t(res)
  rownames(res) <- chroms
  
  res
}

#' Regex expressions to filter 'chromparams.table' column names in transformative HMMs
#' @param type character string with type of regex to get
#' @export
chromparams.rx <- function(type = "ecounts") {
  if (!(type %in% c("ecounts", "edist", "trans")))
    stop("invalid type: ", type, " allowed types: ecounts, edist, trans")
  
  if (type == "ecounts")
    "[.]e[.].[.]2"
  else if (type == "edist")
    "[.]e[.].[.]1"
  else if (type == "trans")
    "[.]t"
}

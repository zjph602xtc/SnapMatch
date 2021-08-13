#' @import optmatch
#' @import rcbalance
#' @import lpSolve
#' @importFrom stats glm.fit rbinom rnorm

#' @export
match_snap <-
  function (x, ...)
    UseMethod("match_snap")

# optimized version
match_internal <- function(logitpsn, id, nt=nt, n = n, controls=1){
  dn <- logitpsn
  matchedn <- rep(NA, nt*controls)
  for (j in 1:controls){
    ind <- seq(j, nt*controls, controls)
    for(i in 1:nt){
      tmp <- which(dn[i,] == min(dn[i,]))
      matchedn[ind[i]] <- tmp
      idn <- id[tmp]
      dn[, id==idn] <- Inf
    }
  }
  matchedn
}


#' @export
match_snap.matrix <- function(x, controls = 1, data, id){
  dat_t <- data[rownames(data) %in% rownames(x),]
  dat_c <- data[rownames(data) %in% colnames(x),]
  # if (dat_c)
  # id is not aggregated, throw error here.

  dat <- rbind(dat_t, dat_c)
  dat_c$myidxxx <- as.numeric(factor(dat_c[,id],unique(dat_c[,id])))
  matched <- match_internal(as.matrix(x), id = dat_c$myidxxx, nt = NROW(dat_t), n = length(unique(dat_c[,'id'])), controls=controls)
  dat$pmxxx <- NA
  dat$pmxxx[matched+NROW(dat_t)] <- as.character(dat_t[,id])
  dat$pmxxx[1:NROW(dat_t)] <- as.character(dat[1:NROW(dat_t),id])
  dat$pmxxx <- factor(dat$pmxxx)
  data <- merge(data[,c('id','trt')],dat[,'pmxxx',drop=F],by=0,sort=FALSE)
  return(data$pmxxx)
}

#' @export
match_snap.formula <- function(x, controls = 1, data, id){
  x <- match_on(x, data=data, method="mahalanobis")
  match_snap.matrix(x, controls = controls, data, id)
}

#' @export
match_snap.glm <- match_snap.formula
#' @export
match_snap.bigglm <- match_snap.formula
#' @export
match_snap.optmatch.dlist <- match_snap.matrix
#' @export
match_snap.InfinitySparseMatrix <- match_snap.matrix
#' @export
match_snap.BlockedInfinitySparseMatrix <- match_snap.matrix

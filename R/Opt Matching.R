#' @import optmatch
#' @import rcbalance
#' @import lpSolve
#' @importFrom stats glm.fit rbinom rnorm

#' @export
optmatch_snap <-
  function (x, ...)
    UseMethod("optmatch_snap")

# optimized version
opt_match_internal <- function(logitpsn, id, nt, n, tol = 1e-6, controls = 1){
  ng <- n # number of subjects in control
  n <- length(id) # number control measure
  if (nt*controls > ng)
    stop("number of treatment * controls > number of control subjects!")
  dn <- logitpsn

  from <- c(rep(1,ng), id+n+nt+1, rep(2:(n+1),each=nt), (n+2):(n+nt+1))
  to <- c(seq(n+nt+2,n+nt+ng+1), (1:n)+1, rep((n+2):(n+nt+1),n), rep(n+nt+ng+2, nt))
  # create edgelist
  edgelist <- list(
    startn = from,
    endn = to,
    ucap = c(rep(1, length(from)-nt), rep(controls, nt)),
    cost = c(rep(0, ng+n), as.numeric(dn), rep(0, nt)),
    b = c(controls*nt,rep(0,n+nt+ng),-controls*nt)
  )
  cost <- edgelist$cost
  if (any(cost != round(cost))) {
    intcost <- round(cost/tol)
    searchtol <- 10^(-c(1:floor(log10(.Machine$integer.max))))
    searchtol <- searchtol[searchtol > tol]
    for (newtol in searchtol) {
      new.intcost <- round(intcost * tol/newtol)
      if (any(new.intcost != intcost * tol/newtol))
        break
      tol <- newtol
      intcost <- new.intcost
    }
    cost <- intcost
  }
  if (any(is.na(as.integer(cost)))) {
    stop("Integer overflow in cost!")
  }
  edgelist$cost <- cost
  solution_r <- callrelax(edgelist)

  # if (solution_r$crash==1 | solution_r$feasible==0){
  #   warnings('RELAX 4 solver fails. We will try the LP solver.')
  #   edgelist <- data.frame(
  #     from = from,
  #     to = to,
  #     ID = seq(1, length(from), 1),
  #     capacity = c(rep(1, length(from)-nt), rep(controls, nt)),
  #     cost = c(rep(0, ng+n), as.numeric(dn), rep(0, nt))
  #   )
  #   t_n <- length(dn)+n+ng+nt
  #   rit <- matrix(c(rep(1:t_n, 2), rep(1, t_n)), ncol = 3)
  #   rit <- rbind(rit, matrix(c(id+t_n, (ng+1):(ng+length(id)), rep(-1, length(id))), ncol = 3))
  #   rit <- rbind(rit, matrix(c(unique(id)+t_n, 1:ng, rep(1, ng)), ncol = 3))
  #   rit <- rbind(rit, matrix(c((1:n)+t_n+ng, (ng+1):(ng+n), rep(1, length(id))), ncol = 3))
  #   rit <- rbind(rit, matrix(c(rep((1:n)+t_n+ng, each = nt), (ng+n+1):(ng+n+n*nt), rep(-1, n*nt)), ncol = 3))
  #   rit <- rbind(rit, matrix(c(rep((t_n+ng+n+1):(t_n+ng+n+nt), n+1), (ng+n+1):(ng+n+nt*n+nt), rep(1, nt*n), rep(-1, nt)), ncol = 3))
  #   rit <- rbind(rit, matrix(c(rep((t_n+ng+n+nt+1), ng), 1:ng, rep(-1, ng)), ncol = 3))
  #   rit <- rbind(rit, matrix(c(rep((t_n+ng+n+nt+2), nt), (ng+n+nt*n+1):(ng+n+nt*n+nt), rep(1, nt)), ncol = 3))
  #   # Run lpSolve to find best solution
  #   solution <- lp(
  #     direction = 'min',
  #     objective.in = edgelist$cost,
  #     const.dir = c(rep('<=', nrow(edgelist)),rep('==', n+nt+ng+2)),
  #     const.rhs = c(rep(1, nrow(edgelist)),rep(0, n+nt+ng),-nt,nt),
  #     dense.const = rit, scale = 0)
  # }
  edgetb <- data.frame(
    from = from,
    to = to
  )
  # if (solution_r$crash==0 & solution_r$feasible==1){
  #   edgetb$flow <- solution_r$x
  # }else{
  #   edgetb$flow <- solution$solution
  # }
  edgetb$flow <- solution_r$x
  edgetb <- edgetb[(ng+n+1):(ng+n+nt*n),]
  edgetb <- edgetb[edgetb$flow==1,]
  edgetb <- edgetb[order(edgetb$to),]
  if (NROW(edgetb)==nt*controls){
    return(edgetb$from-1)
  }else{
    warning('inaccurate result')
    # missing <- setdiff((n+2):(n+nt+1),edgetb$to)
    # matchid <- edgetb$from-1
    # matchid <- unique(id[matchid])
    # ctln[(id %in% matchid)] <- Inf
    # missedid <- missing - n - 1
    # trtn <- trtn[missedid]
    # dn <- outer(trtn,ctln, function(x,y) abs(x - y))
    # matchedn <- rep(NA, length(missedid))
    # for(i in 1:length(missedid)){
    #   matchedn[i] <- which(dn[i,] == min(dn[i,]))[1]
    #   idn <- id[matchedn[i]]
    #   dn[, id==idn] <- Inf
    # }
    # edgetb <- rbind(edgetb[,1:2], data.frame(from=matchedn+1, to=missing))
    # edgetb <- edgetb[order(edgetb$to),]
    return(edgetb$from-1)
  }
}

demo <- function(){
  set.seed(123)
  dat_t <- data.frame(id=paste0('t',1:20), trt=1, X1=rnorm(20,-1), X2=rbinom(20,5,0.2))
  dat_c <- data.frame(id=rep(paste0('c',1:100),each=5), trt=0, X1=rnorm(500,0), X2=rbinom(500,5,0.4))
  dat <- rbind(dat_t,dat_c)
}

#' @export
optmatch_snap.matrix <- function(x, controls = 1, data, id, tol=1e-6){
  dat_t <- data[rownames(data) %in% rownames(x),]
  dat_c <- data[rownames(data) %in% colnames(x),]
  # if (dat_c)
  # id is not aggregated, throw error here.

  dat <- rbind(dat_t, dat_c)
  dat_c$myidxxx <- as.numeric(factor(dat_c[,id],unique(dat_c[,id])))
  x <- as.matrix(x)
  x[is.infinite(x)] <- 20*max(max(x[!is.infinite(x)]))
  matched <- opt_match_internal(x, id = dat_c$myidxxx, nt = NROW(dat_t), n = length(unique(dat_c[,'id'])), controls = controls)
  dat$pmxxx <- NA
  dat$pmxxx[matched+NROW(dat_t)] <- as.character(dat_t[,id])
  dat$pmxxx[1:NROW(dat_t)] <- as.character(dat[1:NROW(dat_t),id])
  dat$pmxxx <- factor(dat$pmxxx)
  data <- merge(data[,c('id','trt')],dat[,'pmxxx',drop=F],by=0,sort=FALSE)
  return(data$pmxxx)
}

#' @export
optmatch_snap.formula <- function(x, controls = 1, data, id, tol = 1e-6){
  x <- match_on(x, data = data, method = "mahalanobis")
  optmatch_snap.matrix(x, controls, data, id, tol)
}

#' @export
optmatch_snap.glm <- optmatch_snap.formula
#' @export
optmatch_snap.bigglm <- optmatch_snap.formula
#' @export
optmatch_snap.optmatch.dlist <- optmatch_snap.matrix
#' @export
optmatch_snap.InfinitySparseMatrix <- optmatch_snap.matrix
#' @export
optmatch_snap.BlockedInfinitySparseMatrix <- optmatch_snap.matrix

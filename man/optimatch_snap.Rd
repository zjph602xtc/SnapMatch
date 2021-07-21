\name{optimatch_snap}
\alias{optimatch_snap}

\title{Snapshot Matching}
\description{
Longitudinal snapshot matching.
}
\usage{
optimatch_snap(
  x,
  data,
  id,
  tol=1e-6)
}

\arguments{
  \item{x}{Any valid input to \code{match_on}. \code{optimatch_snap} will use
\code{x} and any optional arguments to generate a distance before performing
the matching.

Alternatively, a precomputed distance may be entered. A matrix of
non-negative discrepancies, each indicating the permissibility and
desirability of matching the unit corresponding to its row (a 'treatment') to
the unit corresponding to its column (a 'control'); or, better, a distance
specification as produced by \code{\link{match_on}}.}
  \item{data}{\code{data.frame} to use to get order
of the final matching factor. If a \code{data.frame}, the \code{rownames}
are used. Useful to pass if you want to
combine a match (using, e.g., \code{cbind}) with the data that were used to
generate it (for example, in a propensity score matching).}
  \item{id}{A string to indicate which column of the \code{data} is the subject ID.}
  \item{tol}{The minimum difference in distance that \code{optimatch_snap} will differentiate.}
}

\details{The usage of this function is very similar to the \code{fullmatch} in the \pkg{optmatch} package.}

\value{
A vector indicating matched groups.
}

\author{
Tianchen Xu, Yiyue Lou
}

\examples{
dat <- SnapMatch:::demo()

### Match With Covariate Distance
  # calculate the distance matrix
dis <- match_on(trt~X1+X2, data = dat, method="mahalanobis")
  # conduct optimal matching. 'Id' is the id column for each subject
pm <- optimatch_snap(dis, data = dat, id = 'id')
  # plot result
# plot(xBalance(trt ~ X1+X2+strata(pm), data=dat), xlim=c(-2,1), ylim=c(0,3))

### Match With Propensity Score
  # calculate propensity score
prop <- glm(trt ~ X1+X2, data = dat, family = binomial())
  # conduct optimal matching by propensity score. 'Id' is the id column for each subject
pm <- optimatch_snap(prop, data = dat, id = 'id')
  # plot result
# plot(xBalance(trt ~ X1+X2+strata(pm), data=dat), xlim=c(-2,1), ylim=c(0,3))

### Complex Match With ‘match_on’ Function
  # calculate distance by propensity score
dis <- match_on(glm(trt ~ X1+X2, data = dat, family = binomial()))
  # add caliper
dis <- caliper(dis,2) + dis
  # conduct optimal snapshot matching
dis <- as.matrix(dis)
dis[which(is.infinite(dis))] <- 100
pm <- optimatch_snap(dis, data = dat, id = 'id')
  # plot result
# plot(xBalance(trt ~ X1+X2+strata(pm), data=dat), xlim=c(-2,1), ylim=c(0,3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matching}% use one of  RShowDoc("KEYWORDS")
\keyword{propensity score}

#' @name NNT
#' @title Nearest Neighbor Test
#' @description
#' \code{NNT} A multivariate test for equal distributions is based on
#' nearest neighbors.The nearest neighbor (NN) tests are a type of test
#'  based on ordered distances between sample elements, which can be
#'  applied when the distributions are continuous.
#'
#' @param z  a matrix(n*d) of data points from the first and second group.
#'           each row represent a data point;
#' @param ix  a permutation of row indices of z.
#' @param size  a vector of sample sizes;
#' @param custom  the parameter which is customized. For this test, custom
#' means a integer used for the number of neighbor in nearest neighbor test
#'
#' @export
#' @return  Statistics obtained after calculation
#' \item{T_nnt}{The nearest neighbor (NN) test statistics}
NNT <- function(z, ix, size,custom) {
  n1 = size[1]
  n2 = size[2]
  n = n1 + n2
  z = z[ix,]
  o = rep(0, nrow(z))
  z = as.data.frame(cbind(z, o))
  NN = nn2(z, z, custom)
  block1 = NN$nn.idx[1:n1,-1]
  block2 = NN$nn.idx[(n1+1):n,-1]
  i1 = sum(block1 < n1 + 0.5)
  i2 = sum(block2 > n1 + 0.5)
  T_nnt = (i1 + i2) / (custom * n)
  return(T_nnt)
}

#' @name edist
#' @title Energy Distance Test
#' @description
#' \code{edist} A Energy Distance Test for equal distributions
#'
#' @param x  a matrix(n*d) of data points from the first and second group.
#'           each row represent a data point;
#' @param ix  a permutation of row indices of x.
#' @param size  a vector of sample sizes;
#' @param custom  the parameter which is customized. For this test, custom
#' means a integer used for the number of Lp distance
#'
#' @export
#' @return  Statistics obtained after calculation
#' \item{T_edist}{The Energy Distance Test statistics}
edist <- function(x,ix,size,custom) {
  n1 = size[1]
  n2 = size[2]
  ii = ix[1:n1]
  jj = ix[(n1+1):(n1+n2)]
  w = n1 * n2 / (n1 + n2)
  if(is.null(custom)) {
    custom = 2
  }
  dst = as.matrix(dist(x,method='minkowski',p=custom))
  m11 = sum(dst[ii, ii]) / (n1 * n1)
  m22 = sum(dst[jj, jj]) / (n2 * n2)
  m12 = sum(dst[ii, jj]) / (n1 * n2)
  T_edist = w * ((m12 + m12) - (m11 + m22))
  return (T_edist)
}

#' @name HT_sqT
#' @title Hotelling’s T-square test
#' @description
#' \code{HT_sqT} A Hotelling’s T-square test for equal distributions
#'
#' @param x  a matrix(n*d) of data points from the first and second group.
#'          each row represent a data point;
#' @param ix  a permutation of row indices of x.
#' @param size  a vector of sample sizes;
#' @param custom  the parameter which is customized. For this test, custom
#'                means nothing
#'
#' @export
#' @return  Statistics obtained after calculation
#' \item{T_sq}{The EHotelling’s T-square test statistics}
HT_sqT <- function(z,ix,size,custom = NULL) {
  n1 = size[1]
  n2 = size[2]
  n = n1 + n2
  z = z[ix,]
  X = z[1:n1,]
  Y = z[n1+1:n2,]
  X_ = colMeans(X)
  Y_ = colMeans(Y)
  sigma_X = cov(X)
  sigma_Y = cov(Y)
  sigma_hat = ((n1-1)*sigma_X+(n2-1)*sigma_Y)/(n-2)
  T_sq = t(X_-Y_)%*%solve(sigma_hat)%*%(X_-Y_)
  T_sq = T_sq*n1*n2/n
  T_sq = as.double(T_sq)
  return(T_sq)
}

#' @name GraphT
#' @title Graph-based Two Sample Test
#' @description
#' \code{GraphT} A graph-based two sample test for equal distributions
#'
#'
#' @param x  a matrix(n*d) of data points from the first and second group.
#'           each row represent a data point;
#' @param ix  a permutation of row indices of x.
#' @param size  a vector of sample sizes;
#' @param custom  the parameter which is customized. For this test,
#'                custom means a vector used for the number of distance
#'                and the threshold Q in graph-based test
#'
#' @export
#' @return  Statistics obtained after calculation
#' \item{Tr}{The graph-based two sample test statistics}
GraphT <- function(x,ix,size,custom) {
  n1 = size[1]
  n2 = size[2]
  n = n1 + n2
  x = x[ix,]
  ii = ix[1:n1]
  jj = ix[(n1+1):(n1+n2)]
  if(is.null(custom)) {
    custom = c(2,3)
  }
  if(length(custom)<2) {
    custom = c(custom,3)
  }
  dst = as.matrix(dist(x,method='minkowski',p=custom[1]))
  Edge = ifelse((dst < custom[2]) == "TRUE", 1, 0)
  num_edge = n+(sum(Edge)-n)/2
  sum1 = n1+(sum(Edge[ii,ii])-n1)/2
  sum2 = n2+(sum(Edge[jj,jj])-n2)/2
  Tr = (sum1+sum2)/num_edge
  return (Tr)
}

#' @name permT
#' @title Permutation Test
#' @description
#' \code{permT} A permutation test for calculate P-value
#'
#'
#' @param x  a matrix(n*d) of data points from the first and
#'           second group. each row represent a data point;
#' @param B: number of repetitions;
#' @param ix  a permutation of row indices of x.
#' @param size  a vector of sample sizes;
#' @param custom  the parameter which is customized. For this test,
#'                custom means parameter used for fun's custom
#' @param Fun  the function used for this test, including
#'             NNT, edist, HT_sqT, GraphT
#'
#'
#' @export
#' @return  P-value after calculation
#' \item{P}{The P-value}
permT <- function(x,B,ix,size,custom=NULL,Fun) {
  T0 = Fun(x,ix,size,custom)
  Tperm=c(rep(0,B))
  for(b in c(1:B)) {
    per_index = sample(1:nrow(x))
    Tperm[b] = Fun(x,per_index,size,custom)
  }
  P = mean(c(T0,Tperm)>=T0)
  return(P)
}

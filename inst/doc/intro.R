## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(SA23204166)
library(igraph)
library(microbenchmark)
library(Rcpp)

## ----eval=FALSE---------------------------------------------------------------
#  RS_PA <- function(X,T,a) {
#    n=dim(X)[1]
#    sig=svd(X)$d
#    sig0=matrix(rep(0,n*T),nrow = T)
#    for (t in 1:T) {
#      mat=matrix(2*rbinom(n^2,1,0.5)-1,nrow = n)
#      upper_triangle <- upper.tri(mat)
#      s_m <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#      s_m[upper_triangle] <- mat[upper_triangle]
#      s_m[lower.tri(s_m)] <- t(s_m)[lower.tri(s_m)]
#      sig0[t,]=svd(X*s_m)$d
#    }
#    sig1=apply(sig0, 2, quantile, probs = a)
#    k=which(sig < sig1)[1]-1
#    embedding_mat=svd(X)$u[,c(1:k)]
#    return(list(k=k,embedding_mat=embedding_mat))
#  }

## -----------------------------------------------------------------------------
set.seed(0)
mat=matrix(rbinom(50^2,1,0.2),ncol=50)
upper_triangle <- upper.tri(mat)
sim_data <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
sim_data[upper_triangle] <- mat[upper_triangle]
sim_data[lower.tri(sim_data)] <- t(sim_data)[lower.tri(sim_data)]

## -----------------------------------------------------------------------------
print(RS_PA(sim_data,300,0.95))

## -----------------------------------------------------------------------------
data(dolphins)
dolphins_mat=get.adjacency(dolphins,sparse = FALSE)
RS_PA(dolphins_mat,300,0.95)

## -----------------------------------------------------------------------------
data(football)
football_mat=get.adjacency(football,sparse = FALSE)
RS_PA(football_mat,300,0.95)

## -----------------------------------------------------------------------------
data(karate)
karate_mat=get.adjacency(karate,sparse = FALSE)
RS_PA(karate_mat,300,0.95)

## -----------------------------------------------------------------------------
data(polbooks)
polbooks_mat=get.adjacency(polbooks,sparse = FALSE)
RS_PA(polbooks_mat,300,0.95)

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
#  List CommunityC(NumericMatrix embedding_mat, int k, int iter, CharacterVector s) {
#    int n = embedding_mat.nrow();
#    int p = embedding_mat.ncol();
#    IntegerVector seq = seq_len(n);
#    NumericVector x=as<NumericVector>(seq);
#    NumericVector ind=sample(x,k,false,R_NilValue);
#    NumericMatrix centroids(k, p);
#    for (int i = 0; i < k; i++) {
#      centroids(i, _) = embedding_mat(ind[i], _);
#    }
#    IntegerVector assignments(n);
#    for (int iteration = 0; iteration < iter; iteration++) {
#      NumericMatrix distances(n, k);
#      for (int i = 0; i < k; i++) {
#        for (int j = 0; j < n; j++) {
#          distances(j, i) = euclidean_distance(embedding_mat(j, _), centroids(i, _));
#        }
#      }
#      for (int j = 0; j < n; ++j) {
#        assignments[j] = which_min(distances(j, _));
#      }
#      for (int i = 0; i < k; ++i) {
#        NumericVector centroid_update(p, 0.0);
#        int count = 0;
#        for (int j = 0; j < n; ++j) {
#          if (assignments[j] == i) {
#            for (int l = 0; l < p; ++l) {
#              centroid_update[l] += embedding_mat(j, l);
#            }
#            count++;
#          }
#        }
#        if (count > 0) {
#          for (int l = 0; l < p; ++l) {
#            centroids(i, l) = centroid_update[l] / count;
#          }
#        }
#      }
#    }
#    List result(k);
#    for (int i = 0; i < k; ++i) {
#      CharacterVector group_members = s[assignments == i];
#      result[i] = group_members;
#    }
#    return result;
#  }

## -----------------------------------------------------------------------------
set.seed(0)
n <- 100
p <- 3
embedding_mat <- matrix(rnorm(n * p), nrow = n, ncol = p)
s <- paste0("node", 1:n)
k <- 3
iter <- 50
CommunityC(embedding_mat, k, iter, s)

## -----------------------------------------------------------------------------
data(dolphins)
dolphins_mat=get.adjacency(dolphins,sparse = FALSE)
s=get.vertex.attribute(dolphins)$label
k=RS_PA(dolphins_mat,300,0.95)$k
embedding_mat=RS_PA(dolphins_mat,300,0.95)$embedding_mat
iter <- 50
CommunityC(embedding_mat, k, iter, s)


#include <Rcpp.h>
using namespace Rcpp;
double euclidean_distance(NumericVector point1, NumericVector point2) {
  return sqrt(sum(pow(point1 - point2, 2)));
}
//' @title CommunityC
//' @name CommunityC
//' @description Community Detection by k-means clustering using k and embedding matrix obtained in RS_PA
//' @param embedding_mat the data embedding matrix
//' @param k the number of communities
//' @param iter the number of iterations
//' @param s the list of node name
//' @return result the list of node name which in the same community
//' @examples
//' \dontrun{
//' set.seed(0)
//' n <- 100
//' p <- 3
//' embedding_mat <- matrix(rnorm(n * p), nrow = n, ncol = p)
//' s <- paste0("node", 1:n)
//' k <- 3
//' iter <- 50
//' CommunityC(embedding_mat, k, iter, s)
//' data(dolphins)
//' dolphins_mat=get.adjacency(dolphins,sparse = FALSE)
//' s=get.vertex.attribute(dolphins)$label
//' k=RS_PA(dolphins_mat,300,0.95)$k
//' embedding_mat=RS_PA(dolphins_mat,300,0.95)$embedding_mat
//' iter <- 50
//' CommunityC(embedding_mat, k, iter, s)
//' }
//' @export
// [[Rcpp::export]]
List CommunityC(NumericMatrix embedding_mat, int k, int iter, CharacterVector s) {
  int n = embedding_mat.nrow();
  int p = embedding_mat.ncol();
  IntegerVector seq = seq_len(n);
  NumericVector x=as<NumericVector>(seq);
  NumericVector ind=sample(x,k,false,R_NilValue);
  NumericMatrix centroids(k, p);
  for (int i = 0; i < k; i++) {
    centroids(i, _) = embedding_mat(ind[i], _);
  }
  IntegerVector assignments(n);
  for (int iteration = 0; iteration < iter; iteration++) {
    NumericMatrix distances(n, k);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < n; j++) {
        distances(j, i) = euclidean_distance(embedding_mat(j, _), centroids(i, _));
      }
    }
    for (int j = 0; j < n; ++j) {
      assignments[j] = which_min(distances(j, _));
    }
    for (int i = 0; i < k; ++i) {
      NumericVector centroid_update(p, 0.0);
      int count = 0;
      for (int j = 0; j < n; ++j) {
        if (assignments[j] == i) {
          for (int l = 0; l < p; ++l) {
            centroid_update[l] += embedding_mat(j, l);
          }
          count++;
        }
      }
      if (count > 0) {
        for (int l = 0; l < p; ++l) {
          centroids(i, l) = centroid_update[l] / count;
        }
      }
    }
  }
  List result(k);
  for (int i = 0; i < k; ++i) {
    CharacterVector group_members = s[assignments == i];
    result[i] = group_members;
  }
  return result; 
}
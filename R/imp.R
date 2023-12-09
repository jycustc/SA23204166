#' @import boot
#' @import bootstrap
#' @import DAAG
#' @import coda
#' @import stats
#' @import Rcpp
#' @useDynLib SA23204166
#' @importFrom igraph get.adjacency
#' @importFrom igraph get.vertex.attribute
#' @import microbenchmark
NULL


#' @title Dolphins dataset
#' @name dolphins
#' @description This dataset contains an undirected social network of frequent associations between 62 dolphins in a community living off Doubtful Sound, New Zealand, as compiled by Lusseau et al. (2003).
#' @examples
#' \dontrun{
#' data(dolphins)
#' RS_PA(dolphins)
#' }
NULL

#' @title Football dataset
#' @name football
#' @description This dataset contains the network of American football games between Division IA colleges during regular season Fall 2000, as compiled by M. Girvan and M. Newman.
#' @examples
#' \dontrun{
#' data(football)
#' RS_PA(football)
#' }
NULL

#' @title Karate dataset
#' @name karate
#' @description This dataset contains the network of friendships between the 34 members of a karate club at a US university, as described by Wayne Zachary in 1977. 
#' @examples
#' \dontrun{
#' data(karate)
#' RS_PA(karate)
#' }
NULL


#' @title Books about US politics dataset
#' @name polbooks
#' @description Nodes represent books about US politics sold by the online bookseller Amazon.com.  Edges represent frequent co-purchasing of books by the same buyers, as indicated by the "customers who bought this book also bought these other books" feature on Amazon.
#' @examples
#' \dontrun{
#' data(polbooks)
#' RS_PA(polbooks)
#' }
NULL
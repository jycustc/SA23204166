// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CommunityC
List CommunityC(NumericMatrix embedding_mat, int k, int iter, CharacterVector s);
RcppExport SEXP _SA23204166_CommunityC(SEXP embedding_matSEXP, SEXP kSEXP, SEXP iterSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type embedding_mat(embedding_matSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(CommunityC(embedding_mat, k, iter, s));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA23204166_CommunityC", (DL_FUNC) &_SA23204166_CommunityC, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA23204166(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
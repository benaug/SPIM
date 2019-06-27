// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// inoutCpp
bool inoutCpp(NumericVector sx, NumericVector sy, NumericMatrix vertices);
RcppExport SEXP _SPIM_inoutCpp(SEXP sxSEXP, SEXP sySEXP, SEXP verticesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sx(sxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sy(sySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vertices(verticesSEXP);
    rcpp_result_gen = Rcpp::wrap(inoutCpp(sx, sy, vertices));
    return rcpp_result_gen;
END_RCPP
}
// MCMC1
List MCMC1(double lam0, double sigma, NumericMatrix y, NumericMatrix lamd, IntegerVector z, NumericMatrix X, int K, NumericMatrix D, IntegerVector knownvector, NumericMatrix s, NumericVector psi, NumericVector xlim, NumericVector ylim, bool useverts, NumericMatrix vertices, double proplam0, double propsigma, double propsx, double propsy, int niter, int nburn, int nthin, int obstype, IntegerVector tf, bool storeLatent);
RcppExport SEXP _SPIM_MCMC1(SEXP lam0SEXP, SEXP sigmaSEXP, SEXP ySEXP, SEXP lamdSEXP, SEXP zSEXP, SEXP XSEXP, SEXP KSEXP, SEXP DSEXP, SEXP knownvectorSEXP, SEXP sSEXP, SEXP psiSEXP, SEXP xlimSEXP, SEXP ylimSEXP, SEXP usevertsSEXP, SEXP verticesSEXP, SEXP proplam0SEXP, SEXP propsigmaSEXP, SEXP propsxSEXP, SEXP propsySEXP, SEXP niterSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP obstypeSEXP, SEXP tfSEXP, SEXP storeLatentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lamd(lamdSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type knownvector(knownvectorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xlim(xlimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ylim(ylimSEXP);
    Rcpp::traits::input_parameter< bool >::type useverts(usevertsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vertices(verticesSEXP);
    Rcpp::traits::input_parameter< double >::type proplam0(proplam0SEXP);
    Rcpp::traits::input_parameter< double >::type propsigma(propsigmaSEXP);
    Rcpp::traits::input_parameter< double >::type propsx(propsxSEXP);
    Rcpp::traits::input_parameter< double >::type propsy(propsySEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type obstype(obstypeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tf(tfSEXP);
    Rcpp::traits::input_parameter< bool >::type storeLatent(storeLatentSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC1(lam0, sigma, y, lamd, z, X, K, D, knownvector, s, psi, xlim, ylim, useverts, vertices, proplam0, propsigma, propsx, propsy, niter, nburn, nthin, obstype, tf, storeLatent));
    return rcpp_result_gen;
END_RCPP
}
// MCMC2side
List MCMC2side(double lam01, double lam02, double sigma, NumericMatrix lamd1, NumericMatrix lamd2, NumericMatrix y_both, NumericMatrix y_left_true, NumericMatrix y_right_true, NumericMatrix y_left_obs, NumericMatrix y_right_obs, IntegerVector z, NumericMatrix X, IntegerVector tf, NumericMatrix D, int Nfixed, IntegerVector knownvector, IntegerVector ID_L, IntegerVector ID_R, int swap, double swaptol, NumericMatrix s, NumericVector psi, NumericVector xlim, NumericVector ylim, bool useverts, NumericMatrix vertices, double proplam01, double proplam02, double propsigma, double propsx, double propsy, int niter, int nburn, int nthin, LogicalVector updates, bool storeLatent);
RcppExport SEXP _SPIM_MCMC2side(SEXP lam01SEXP, SEXP lam02SEXP, SEXP sigmaSEXP, SEXP lamd1SEXP, SEXP lamd2SEXP, SEXP y_bothSEXP, SEXP y_left_trueSEXP, SEXP y_right_trueSEXP, SEXP y_left_obsSEXP, SEXP y_right_obsSEXP, SEXP zSEXP, SEXP XSEXP, SEXP tfSEXP, SEXP DSEXP, SEXP NfixedSEXP, SEXP knownvectorSEXP, SEXP ID_LSEXP, SEXP ID_RSEXP, SEXP swapSEXP, SEXP swaptolSEXP, SEXP sSEXP, SEXP psiSEXP, SEXP xlimSEXP, SEXP ylimSEXP, SEXP usevertsSEXP, SEXP verticesSEXP, SEXP proplam01SEXP, SEXP proplam02SEXP, SEXP propsigmaSEXP, SEXP propsxSEXP, SEXP propsySEXP, SEXP niterSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP updatesSEXP, SEXP storeLatentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lam01(lam01SEXP);
    Rcpp::traits::input_parameter< double >::type lam02(lam02SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lamd1(lamd1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lamd2(lamd2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y_both(y_bothSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y_left_true(y_left_trueSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y_right_true(y_right_trueSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y_left_obs(y_left_obsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y_right_obs(y_right_obsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tf(tfSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type Nfixed(NfixedSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type knownvector(knownvectorSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ID_L(ID_LSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ID_R(ID_RSEXP);
    Rcpp::traits::input_parameter< int >::type swap(swapSEXP);
    Rcpp::traits::input_parameter< double >::type swaptol(swaptolSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xlim(xlimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ylim(ylimSEXP);
    Rcpp::traits::input_parameter< bool >::type useverts(usevertsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vertices(verticesSEXP);
    Rcpp::traits::input_parameter< double >::type proplam01(proplam01SEXP);
    Rcpp::traits::input_parameter< double >::type proplam02(proplam02SEXP);
    Rcpp::traits::input_parameter< double >::type propsigma(propsigmaSEXP);
    Rcpp::traits::input_parameter< double >::type propsx(propsxSEXP);
    Rcpp::traits::input_parameter< double >::type propsy(propsySEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type updates(updatesSEXP);
    Rcpp::traits::input_parameter< bool >::type storeLatent(storeLatentSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC2side(lam01, lam02, sigma, lamd1, lamd2, y_both, y_left_true, y_right_true, y_left_obs, y_right_obs, z, X, tf, D, Nfixed, knownvector, ID_L, ID_R, swap, swaptol, s, psi, xlim, ylim, useverts, vertices, proplam01, proplam02, propsigma, propsx, propsy, niter, nburn, nthin, updates, storeLatent));
    return rcpp_result_gen;
END_RCPP
}
// MCMCtf2
List MCMCtf2(double lam01, double lam02, double sigma, NumericMatrix lamd1, NumericMatrix lamd2, NumericMatrix y_both, arma::cube y_left_true, arma::cube y_right_true, arma::cube y_left_obs, arma::cube y_right_obs, IntegerVector z, NumericMatrix X, IntegerMatrix tf1, IntegerVector tf2, NumericMatrix D, int Nfixed, IntegerVector knownvector, IntegerVector ID_L, IntegerVector ID_R, int swap, double swaptol, NumericMatrix s, NumericVector psi, NumericVector xlim, NumericVector ylim, bool useverts, NumericMatrix vertices, double proplam01, double proplam02, double propsigma, double propsx, double propsy, int niter, int nburn, int nthin, LogicalVector updates, bool storeLatent);
RcppExport SEXP _SPIM_MCMCtf2(SEXP lam01SEXP, SEXP lam02SEXP, SEXP sigmaSEXP, SEXP lamd1SEXP, SEXP lamd2SEXP, SEXP y_bothSEXP, SEXP y_left_trueSEXP, SEXP y_right_trueSEXP, SEXP y_left_obsSEXP, SEXP y_right_obsSEXP, SEXP zSEXP, SEXP XSEXP, SEXP tf1SEXP, SEXP tf2SEXP, SEXP DSEXP, SEXP NfixedSEXP, SEXP knownvectorSEXP, SEXP ID_LSEXP, SEXP ID_RSEXP, SEXP swapSEXP, SEXP swaptolSEXP, SEXP sSEXP, SEXP psiSEXP, SEXP xlimSEXP, SEXP ylimSEXP, SEXP usevertsSEXP, SEXP verticesSEXP, SEXP proplam01SEXP, SEXP proplam02SEXP, SEXP propsigmaSEXP, SEXP propsxSEXP, SEXP propsySEXP, SEXP niterSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP updatesSEXP, SEXP storeLatentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lam01(lam01SEXP);
    Rcpp::traits::input_parameter< double >::type lam02(lam02SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lamd1(lamd1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lamd2(lamd2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y_both(y_bothSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_left_true(y_left_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_right_true(y_right_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_left_obs(y_left_obsSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_right_obs(y_right_obsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type tf1(tf1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tf2(tf2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type Nfixed(NfixedSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type knownvector(knownvectorSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ID_L(ID_LSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ID_R(ID_RSEXP);
    Rcpp::traits::input_parameter< int >::type swap(swapSEXP);
    Rcpp::traits::input_parameter< double >::type swaptol(swaptolSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xlim(xlimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ylim(ylimSEXP);
    Rcpp::traits::input_parameter< bool >::type useverts(usevertsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vertices(verticesSEXP);
    Rcpp::traits::input_parameter< double >::type proplam01(proplam01SEXP);
    Rcpp::traits::input_parameter< double >::type proplam02(proplam02SEXP);
    Rcpp::traits::input_parameter< double >::type propsigma(propsigmaSEXP);
    Rcpp::traits::input_parameter< double >::type propsx(propsxSEXP);
    Rcpp::traits::input_parameter< double >::type propsy(propsySEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type updates(updatesSEXP);
    Rcpp::traits::input_parameter< bool >::type storeLatent(storeLatentSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMCtf2(lam01, lam02, sigma, lamd1, lamd2, y_both, y_left_true, y_right_true, y_left_obs, y_right_obs, z, X, tf1, tf2, D, Nfixed, knownvector, ID_L, ID_R, swap, swaptol, s, psi, xlim, ylim, useverts, vertices, proplam01, proplam02, propsigma, propsx, propsy, niter, nburn, nthin, updates, storeLatent));
    return rcpp_result_gen;
END_RCPP
}
// HazACup
List HazACup(NumericMatrix y, NumericMatrix pd, NumericMatrix X, NumericMatrix D, NumericMatrix s, IntegerVector z, NumericMatrix ll_y_curr, double beta0, double beta1, double sigma, int M, int J, IntegerVector tf, double props, NumericVector xlim, NumericVector ylim);
RcppExport SEXP _SPIM_HazACup(SEXP ySEXP, SEXP pdSEXP, SEXP XSEXP, SEXP DSEXP, SEXP sSEXP, SEXP zSEXP, SEXP ll_y_currSEXP, SEXP beta0SEXP, SEXP beta1SEXP, SEXP sigmaSEXP, SEXP MSEXP, SEXP JSEXP, SEXP tfSEXP, SEXP propsSEXP, SEXP xlimSEXP, SEXP ylimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pd(pdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ll_y_curr(ll_y_currSEXP);
    Rcpp::traits::input_parameter< double >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tf(tfSEXP);
    Rcpp::traits::input_parameter< double >::type props(propsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xlim(xlimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ylim(ylimSEXP);
    rcpp_result_gen = Rcpp::wrap(HazACup(y, pd, X, D, s, z, ll_y_curr, beta0, beta1, sigma, M, J, tf, props, xlim, ylim));
    return rcpp_result_gen;
END_RCPP
}
// intlikRcpp
double intlikRcpp(NumericVector parm, NumericMatrix ymat, IntegerMatrix X, int K, NumericMatrix G, NumericMatrix D, int n);
RcppExport SEXP _SPIM_intlikRcpp(SEXP parmSEXP, SEXP ymatSEXP, SEXP XSEXP, SEXP KSEXP, SEXP GSEXP, SEXP DSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parm(parmSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(intlikRcpp(parm, ymat, X, K, G, D, n));
    return rcpp_result_gen;
END_RCPP
}
// findPossible2D
LogicalVector findPossible2D(IntegerVector z, IntegerMatrix G_true, IntegerVector G_obs_true, int M, int ncat);
RcppExport SEXP _SPIM_findPossible2D(SEXP zSEXP, SEXP G_trueSEXP, SEXP G_obs_trueSEXP, SEXP MSEXP, SEXP ncatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type G_true(G_trueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type G_obs_true(G_obs_trueSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type ncat(ncatSEXP);
    rcpp_result_gen = Rcpp::wrap(findPossible2D(z, G_true, G_obs_true, M, ncat));
    return rcpp_result_gen;
END_RCPP
}
// arma_setdiff
Rcpp::NumericVector arma_setdiff(arma::uvec& x, arma::uvec& y);
RcppExport SEXP _SPIM_arma_setdiff(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(arma_setdiff(x, y));
    return rcpp_result_gen;
END_RCPP
}
// indLL
List indLL(IntegerVector z, arma::cube yI, arma::cube yP, double lambdaI, IntegerVector ID);
RcppExport SEXP _SPIM_indLL(SEXP zSEXP, SEXP yISEXP, SEXP yPSEXP, SEXP lambdaISEXP, SEXP IDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type yI(yISEXP);
    Rcpp::traits::input_parameter< arma::cube >::type yP(yPSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaI(lambdaISEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ID(IDSEXP);
    rcpp_result_gen = Rcpp::wrap(indLL(z, yI, yP, lambdaI, ID));
    return rcpp_result_gen;
END_RCPP
}
// siteup
List siteup(IntegerVector z1, IntegerVector z2, arma::cube yI, arma::cube yP, arma::cube samp_mat, double lambdaI, IntegerVector packID, arma::cube ll_yi, arma::cube ll_yp, NumericMatrix lamd);
RcppExport SEXP _SPIM_siteup(SEXP z1SEXP, SEXP z2SEXP, SEXP yISEXP, SEXP yPSEXP, SEXP samp_matSEXP, SEXP lambdaISEXP, SEXP packIDSEXP, SEXP ll_yiSEXP, SEXP ll_ypSEXP, SEXP lamdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z2(z2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type yI(yISEXP);
    Rcpp::traits::input_parameter< arma::cube >::type yP(yPSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type samp_mat(samp_matSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaI(lambdaISEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type packID(packIDSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_yi(ll_yiSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_yp(ll_ypSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lamd(lamdSEXP);
    rcpp_result_gen = Rcpp::wrap(siteup(z1, z2, yI, yP, samp_mat, lambdaI, packID, ll_yi, ll_yp, lamd));
    return rcpp_result_gen;
END_RCPP
}
// calcllyevent
arma::cube calcllyevent(arma::cube y_true, arma::cube y_event, arma::cube Xcov3D, arma::cube pi1, arma::cube pi2, arma::cube pi3, IntegerVector z, int M, int J, int K);
RcppExport SEXP _SPIM_calcllyevent(SEXP y_trueSEXP, SEXP y_eventSEXP, SEXP Xcov3DSEXP, SEXP pi1SEXP, SEXP pi2SEXP, SEXP pi3SEXP, SEXP zSEXP, SEXP MSEXP, SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type y_true(y_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_event(y_eventSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Xcov3D(Xcov3DSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi3(pi3SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(calcllyevent(y_true, y_event, Xcov3D, pi1, pi2, pi3, z, M, J, K));
    return rcpp_result_gen;
END_RCPP
}
// calcllyevent2D
NumericMatrix calcllyevent2D(IntegerMatrix y_true, IntegerMatrix y_event, IntegerMatrix Xcov3D, NumericMatrix pi1, NumericMatrix pi2, NumericMatrix pi3, int J, int K);
RcppExport SEXP _SPIM_calcllyevent2D(SEXP y_trueSEXP, SEXP y_eventSEXP, SEXP Xcov3DSEXP, SEXP pi1SEXP, SEXP pi2SEXP, SEXP pi3SEXP, SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type y_true(y_trueSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y_event(y_eventSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Xcov3D(Xcov3DSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pi3(pi3SEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(calcllyevent2D(y_true, y_event, Xcov3D, pi1, pi2, pi3, J, K));
    return rcpp_result_gen;
END_RCPP
}
// occDepUp
List occDepUp(IntegerMatrix occidx, arma::cube pd, int upOcc, arma::cube y_true, arma::cube y_Occ_true, arma::cube y_SCR, arma::cube y_event, arma::cube pi1, arma::cube pi2, arma::cube pi3, arma::cube ll_y, arma::cube ll_y_event, IntegerVector z, int M, int J, int K);
RcppExport SEXP _SPIM_occDepUp(SEXP occidxSEXP, SEXP pdSEXP, SEXP upOccSEXP, SEXP y_trueSEXP, SEXP y_Occ_trueSEXP, SEXP y_SCRSEXP, SEXP y_eventSEXP, SEXP pi1SEXP, SEXP pi2SEXP, SEXP pi3SEXP, SEXP ll_ySEXP, SEXP ll_y_eventSEXP, SEXP zSEXP, SEXP MSEXP, SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type occidx(occidxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pd(pdSEXP);
    Rcpp::traits::input_parameter< int >::type upOcc(upOccSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_true(y_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_Occ_true(y_Occ_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_SCR(y_SCRSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_event(y_eventSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi3(pi3SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y(ll_ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y_event(ll_y_eventSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(occDepUp(occidx, pd, upOcc, y_true, y_Occ_true, y_SCR, y_event, pi1, pi2, pi3, ll_y, ll_y_event, z, M, J, K));
    return rcpp_result_gen;
END_RCPP
}
// occDepUpMbradius
List occDepUpMbradius(IntegerMatrix occidx, NumericMatrix s, NumericMatrix X, arma::cube pd, double upDist, arma::cube y_true, arma::cube y_Occ_true, arma::cube y_SCR, arma::cube y_event, arma::cube pi1, arma::cube pi2, arma::cube pi3, arma::cube ll_y, arma::cube ll_y_event, arma::cube state, arma::cube linp_p0, arma::cube p0, NumericVector beta_p0, IntegerVector sex, arma::cube Xcov3D1, arma::cube Xcov3D2, arma::cube D, IntegerVector z, NumericVector sigma, int M, int J, int K);
RcppExport SEXP _SPIM_occDepUpMbradius(SEXP occidxSEXP, SEXP sSEXP, SEXP XSEXP, SEXP pdSEXP, SEXP upDistSEXP, SEXP y_trueSEXP, SEXP y_Occ_trueSEXP, SEXP y_SCRSEXP, SEXP y_eventSEXP, SEXP pi1SEXP, SEXP pi2SEXP, SEXP pi3SEXP, SEXP ll_ySEXP, SEXP ll_y_eventSEXP, SEXP stateSEXP, SEXP linp_p0SEXP, SEXP p0SEXP, SEXP beta_p0SEXP, SEXP sexSEXP, SEXP Xcov3D1SEXP, SEXP Xcov3D2SEXP, SEXP DSEXP, SEXP zSEXP, SEXP sigmaSEXP, SEXP MSEXP, SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type occidx(occidxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pd(pdSEXP);
    Rcpp::traits::input_parameter< double >::type upDist(upDistSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_true(y_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_Occ_true(y_Occ_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_SCR(y_SCRSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_event(y_eventSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi3(pi3SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y(ll_ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y_event(ll_y_eventSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type state(stateSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type linp_p0(linp_p0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_p0(beta_p0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sex(sexSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Xcov3D1(Xcov3D1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Xcov3D2(Xcov3D2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type D(DSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(occDepUpMbradius(occidx, s, X, pd, upDist, y_true, y_Occ_true, y_SCR, y_event, pi1, pi2, pi3, ll_y, ll_y_event, state, linp_p0, p0, beta_p0, sex, Xcov3D1, Xcov3D2, D, z, sigma, M, J, K));
    return rcpp_result_gen;
END_RCPP
}
// occDepUpMb4Dradius
List occDepUpMb4Dradius(IntegerMatrix occidx, NumericMatrix s, NumericMatrix X, arma::cube pd, double upDist, arma::cube y_true, arma::cube y_Occ_true, arma::cube y_SCR, arma::cube y_event, arma::cube pi1, arma::cube pi2, arma::cube pi3, arma::cube ll_y, arma::cube ll_y_event, arma::cube state, arma::cube linp_p0, arma::cube p0, NumericVector beta_p0, IntegerVector sex, arma::cube Xcov3D1, arma::cube Xcov3D2, arma::cube D, IntegerVector z, NumericVector sigma, int M, int J, int K);
RcppExport SEXP _SPIM_occDepUpMb4Dradius(SEXP occidxSEXP, SEXP sSEXP, SEXP XSEXP, SEXP pdSEXP, SEXP upDistSEXP, SEXP y_trueSEXP, SEXP y_Occ_trueSEXP, SEXP y_SCRSEXP, SEXP y_eventSEXP, SEXP pi1SEXP, SEXP pi2SEXP, SEXP pi3SEXP, SEXP ll_ySEXP, SEXP ll_y_eventSEXP, SEXP stateSEXP, SEXP linp_p0SEXP, SEXP p0SEXP, SEXP beta_p0SEXP, SEXP sexSEXP, SEXP Xcov3D1SEXP, SEXP Xcov3D2SEXP, SEXP DSEXP, SEXP zSEXP, SEXP sigmaSEXP, SEXP MSEXP, SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type occidx(occidxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pd(pdSEXP);
    Rcpp::traits::input_parameter< double >::type upDist(upDistSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_true(y_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_Occ_true(y_Occ_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_SCR(y_SCRSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_event(y_eventSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi3(pi3SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y(ll_ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y_event(ll_y_eventSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type state(stateSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type linp_p0(linp_p0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_p0(beta_p0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sex(sexSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Xcov3D1(Xcov3D1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Xcov3D2(Xcov3D2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type D(DSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(occDepUpMb4Dradius(occidx, s, X, pd, upDist, y_true, y_Occ_true, y_SCR, y_event, pi1, pi2, pi3, ll_y, ll_y_event, state, linp_p0, p0, beta_p0, sex, Xcov3D1, Xcov3D2, D, z, sigma, M, J, K));
    return rcpp_result_gen;
END_RCPP
}
// occDepUpRadius
List occDepUpRadius(IntegerMatrix occidx, NumericMatrix s, NumericMatrix X, arma::cube pd, double upDist, arma::cube y_true, arma::cube y_Occ_true, arma::cube y_SCR, arma::cube y_event, arma::cube pi1, arma::cube pi2, arma::cube pi3, arma::cube ll_y, arma::cube ll_y_event, IntegerVector z, int M, int J, int K);
RcppExport SEXP _SPIM_occDepUpRadius(SEXP occidxSEXP, SEXP sSEXP, SEXP XSEXP, SEXP pdSEXP, SEXP upDistSEXP, SEXP y_trueSEXP, SEXP y_Occ_trueSEXP, SEXP y_SCRSEXP, SEXP y_eventSEXP, SEXP pi1SEXP, SEXP pi2SEXP, SEXP pi3SEXP, SEXP ll_ySEXP, SEXP ll_y_eventSEXP, SEXP zSEXP, SEXP MSEXP, SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type occidx(occidxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pd(pdSEXP);
    Rcpp::traits::input_parameter< double >::type upDist(upDistSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_true(y_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_Occ_true(y_Occ_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_SCR(y_SCRSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y_event(y_eventSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type pi3(pi3SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y(ll_ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ll_y_event(ll_y_eventSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(occDepUpRadius(occidx, s, X, pd, upDist, y_true, y_Occ_true, y_SCR, y_event, pi1, pi2, pi3, ll_y, ll_y_event, z, M, J, K));
    return rcpp_result_gen;
END_RCPP
}
// typeCount
List typeCount(int ncat, int nrep, IntegerVector nlevels, List locustype, List ptype, IntegerMatrix G_obs_true, arma::cube G_obs);
RcppExport SEXP _SPIM_typeCount(SEXP ncatSEXP, SEXP nrepSEXP, SEXP nlevelsSEXP, SEXP locustypeSEXP, SEXP ptypeSEXP, SEXP G_obs_trueSEXP, SEXP G_obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nlevels(nlevelsSEXP);
    Rcpp::traits::input_parameter< List >::type locustype(locustypeSEXP);
    Rcpp::traits::input_parameter< List >::type ptype(ptypeSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type G_obs_true(G_obs_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type G_obs(G_obsSEXP);
    rcpp_result_gen = Rcpp::wrap(typeCount(ncat, nrep, nlevels, locustype, ptype, G_obs_true, G_obs));
    return rcpp_result_gen;
END_RCPP
}
// typeCountSamptype
List typeCountSamptype(int ncat, int nrep, IntegerVector nlevels, List locustype, List ptype, IntegerMatrix G_obs_true, arma::cube G_obs, IntegerVector samptype, int samplevels);
RcppExport SEXP _SPIM_typeCountSamptype(SEXP ncatSEXP, SEXP nrepSEXP, SEXP nlevelsSEXP, SEXP locustypeSEXP, SEXP ptypeSEXP, SEXP G_obs_trueSEXP, SEXP G_obsSEXP, SEXP samptypeSEXP, SEXP samplevelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nlevels(nlevelsSEXP);
    Rcpp::traits::input_parameter< List >::type locustype(locustypeSEXP);
    Rcpp::traits::input_parameter< List >::type ptype(ptypeSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type G_obs_true(G_obs_trueSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type G_obs(G_obsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type samptype(samptypeSEXP);
    Rcpp::traits::input_parameter< int >::type samplevels(samplevelsSEXP);
    rcpp_result_gen = Rcpp::wrap(typeCountSamptype(ncat, nrep, nlevels, locustype, ptype, G_obs_true, G_obs, samptype, samplevels));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SPIM_inoutCpp", (DL_FUNC) &_SPIM_inoutCpp, 3},
    {"_SPIM_MCMC1", (DL_FUNC) &_SPIM_MCMC1, 25},
    {"_SPIM_MCMC2side", (DL_FUNC) &_SPIM_MCMC2side, 36},
    {"_SPIM_MCMCtf2", (DL_FUNC) &_SPIM_MCMCtf2, 37},
    {"_SPIM_HazACup", (DL_FUNC) &_SPIM_HazACup, 16},
    {"_SPIM_intlikRcpp", (DL_FUNC) &_SPIM_intlikRcpp, 7},
    {"_SPIM_findPossible2D", (DL_FUNC) &_SPIM_findPossible2D, 5},
    {"_SPIM_arma_setdiff", (DL_FUNC) &_SPIM_arma_setdiff, 2},
    {"_SPIM_indLL", (DL_FUNC) &_SPIM_indLL, 5},
    {"_SPIM_siteup", (DL_FUNC) &_SPIM_siteup, 10},
    {"_SPIM_calcllyevent", (DL_FUNC) &_SPIM_calcllyevent, 10},
    {"_SPIM_calcllyevent2D", (DL_FUNC) &_SPIM_calcllyevent2D, 8},
    {"_SPIM_occDepUp", (DL_FUNC) &_SPIM_occDepUp, 16},
    {"_SPIM_occDepUpMbradius", (DL_FUNC) &_SPIM_occDepUpMbradius, 27},
    {"_SPIM_occDepUpMb4Dradius", (DL_FUNC) &_SPIM_occDepUpMb4Dradius, 27},
    {"_SPIM_occDepUpRadius", (DL_FUNC) &_SPIM_occDepUpRadius, 18},
    {"_SPIM_typeCount", (DL_FUNC) &_SPIM_typeCount, 7},
    {"_SPIM_typeCountSamptype", (DL_FUNC) &_SPIM_typeCountSamptype, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_SPIM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

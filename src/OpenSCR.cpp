////////////////////////////////////in polygon functions/////////////////////
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
bool intersectCppOpen(NumericVector sx, NumericVector sy,NumericVector vertex1, NumericVector vertex2) {
  NumericVector swap;
  double m_blue;
  double m_red;
  bool out;
  if((sy(0)==vertex1[1])|(sy(0)==vertex2[1])){
    sy(0)=sy(0)+0.000001;
  }
  if(vertex1[1]>vertex2[1]){
    swap=vertex1;
    vertex1=vertex2;
    vertex2=swap;
  }
  if((sy(0)<vertex1[1])|(sy(0)>vertex2[1])){
    out=FALSE;
  }else if((sx(0) > vertex1[0]) & (sx(0)> vertex2[0])){
    out=FALSE;
  }else{
    if((sx(0) < vertex1[0]) & (sx(0) < vertex2[0])){
      out=TRUE;
    }else{
      if(vertex1[0]!=vertex2[0]){
        m_red=(vertex2(1)-vertex1(1))/(vertex2(0)-vertex1(0));
      }else{
        m_red=1000000000000000000;
      }
      if(vertex1[0]!=sx(0)){
        m_blue=(sy(0)-vertex1(1))/(sx(0)-vertex1(0));
      }else{
        m_blue=1000000000000000000;
      }
      if(m_blue>=m_red){
        out=TRUE;
      }else{
        out=FALSE;
      }
    }
  }
  return out;
}


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
bool inoutCppOpen(NumericVector sx,NumericVector sy,NumericMatrix vertices) {
  int count=0;
  int I=vertices.nrow();
  for(int i=0; i<(I-1); i++) {
    if(intersectCppOpen(sx,sy,vertices(i,_),vertices(i+1,_))){
      count=count+1;
    }
  }
  bool to_return;
  if(count % 2 != 0){
    to_return = true;
  }else{
    to_return = false;
  }
  return to_return;
}




//calculate lamd all years at once
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::cube calclamd3D(NumericVector lam0, NumericVector sigma, arma::cube D) {
  //Preallocate
  int M = size(D)[0];
  int J = size(D)[1];
  int t = size(D)[2];
  arma::cube lamd(M,J,t);
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      for(int k=0; k<t; k++){
        lamd(i,j,k)=lam0[k]*exp(-D(i,j,k)*D(i,j,k)/(2*sigma[k]*sigma[k]));
      }
    }
  }
  return lamd;
}
//calculate lamd one at a time
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat calclamd2D(NumericVector lam0, NumericVector sigma, arma::cube D,int idx) {
  //Preallocate
  int M = size(D)[0];
  int J = size(D)[1];
  arma::mat lamd(M,J);
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      lamd(i,j)=lam0(0)*exp(-D(i,j,idx)*D(i,j,idx)/(2*sigma(0)*sigma(0)));
    }
  }
  return lamd;
}

//df likelihood calculation
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::cube OpenSCRlik(IntegerMatrix z,arma::cube lamd,arma::cube y,NumericVector K,NumericVector Xidx) {
  //Preallocate
  int M = size(lamd)[0];
  int J = size(lamd)[1];
  int t = size(lamd)[2];
  arma::cube pd(M,J,t);
  arma::cube ll_y(M,J,t);
  //  Calculate likelihood
  for(int l=0; l<t; l++){
    for(int i=0; i<M; i++) {
      if(z(i,l)==1){ //if in pop
        for(int j=0; j<Xidx(l); j++){
          pd(i,j,l)=1-exp(-lamd(i,j,l));
          ll_y(i,j,l)=y(i,j,l)*log(pd(i,j,l))+(K[l]-y(i,j,l))*log(1-pd(i,j,l));
        }
      }
    }
  }
  return ll_y;
}

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double sum_ll_y(IntegerMatrix z,arma::cube ll_y,NumericVector Xidx) {
  //Preallocate
  int M = size(ll_y)[0];
  int t = size(ll_y)[2];
  double ll_sum=0;
  //  Calculate likelihood
  for(int l=0; l<t; l++){
    for(int i=0; i<M; i++) {
      if(z(i,l)==1){ //if in pop
        for(int j=0; j<Xidx(l); j++){
          if(ll_y(i,j,l)==ll_y(i,j,l)){
            ll_sum+=ll_y(i,j,l);
          }
        }
      }
    }
  }
  return ll_sum;
}

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
NumericVector sum_ll_y_t(IntegerMatrix z,arma::cube ll_y,NumericVector Xidx) {
  //Preallocate
  int M = size(ll_y)[0];
  int t = size(ll_y)[2];
  NumericVector ll_sum(3);
  //  Calculate likelihood
  for(int l=0; l<t; l++){
    for(int i=0; i<M; i++) {
      if(z(i,l)==1){ //if in pop
        for(int j=0; j<Xidx(l); j++){
          if(ll_y(i,j,l)==ll_y(i,j,l)){
            ll_sum(l)+=ll_y(i,j,l);
          }
        }
      }
    }
  }
  return ll_sum;
}

// update lam01, lam02, and sigma

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateDfun3D(NumericVector lam0, NumericVector sigma, arma::cube lamd, arma::cube y,
                  IntegerMatrix z,NumericVector Xidx, NumericVector K,arma::cube D) {
  RNGScope scope;
  //Preallocate
  int M = size(lamd)[0];
  int J = size(lamd)[1];
  int t = size(lamd)[2];
  NumericVector lam0cand(t);
  NumericVector sigmacand(t);
  NumericVector rand;
  NumericVector rand2;
  arma::cube pd=zeros<cube>(M,J,t);
  arma::cube ll_y_curr=zeros<cube>(M,J,t);
  arma::cube lamdcand=zeros<cube>(M,J,t);
  arma::cube pdcand=zeros<cube>(M,J,t);
  arma::cube ll_y_cand=zeros<cube>(M,J,t);
  double likcurr=0;
  double likcand=0;
  double explldiff=0;
  //Structures for dealing with year specific parms
  NumericVector likcurr2D(t);
  NumericVector likcand2D(t);
  NumericVector sigmause(t);
  NumericVector lam0use(t);
  for(int i=0; i<t; i++) {
    if(sigma.size()==1){
      sigmause(i)=sigma(0);
    }else{
      sigmause(i)=sigma(i);
    }
  }
  for(int i=0; i<t; i++) {
    if(lam0.size()==1){
      lam0use(i)=lam0(0);
    }else{
      lam0use(i)=lam0(i);
    }
  }
  //  Calculate and sum current likelihood
  for(int l=0; l<t; l++){
    likcurr2D(l)=0;
    for(int i=0; i<M; i++) {
      for(int j=0; j<Xidx(l); j++){
        pd(i,j,l)=1-exp(-lamd(i,j,l));
        ll_y_curr(i,j,l)=z(i,l)*(y(i,j,l)*log(pd(i,j,l))+(K[l]-y(i,j,l))*log(1-pd(i,j,l)));
        if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
          likcurr2D(l)+=ll_y_curr(i,j,l); //year-specific components
        }
      }
    }
    likcurr+=likcurr2D(l); //full likelihood sum
  }

  //Update lam0
  if(lam0.size()==1){//fixed lambda
    rand=Rcpp::rnorm(1,lam0(0),0.025);
    if(rand(0) > 0){
      likcand=0;
      lam0cand(0)=rand(0);
      //  Update lamd and calculate cand likelihood
      for(int l=0; l<t; l++){
        likcand2D(l)=0;
        for(int i=0; i<M; i++) {
          for(int j=0; j<Xidx(l); j++){
            lamdcand(i,j,l)=lam0cand(0)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmause(l)*sigmause(l)));
            pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
            ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
            if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
              likcand2D(l)+=ll_y_cand(i,j,l);
            }
          }
        }
        likcand+=likcand2D(l);
      }
      rand2=Rcpp::runif(1);
      explldiff = exp(likcand-likcurr);
      if(rand2(0)<explldiff){
        lam0(0)=lam0cand(0);
        lamd=lamdcand;
        pd=pdcand;
        ll_y_curr=ll_y_cand;
        likcurr=likcand;
        likcurr2D=likcand2D;
      }
    }
  }else{//year-specific lambdas
    likcurr=0;
    for(int l=0; l<t; l++){
      likcand2D(l)=0;
      rand=Rcpp::rnorm(1,lam0(l),0.025);
      if(rand(0) > 0){
        lam0cand(l)=rand(0);
        //  Calculate likelihood
        for(int i=0; i<M; i++) {
          for(int j=0; j<Xidx(l); j++){
            lamdcand(i,j,l)=lam0cand(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmause(l)*sigmause(l)));
            pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
            ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
            if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
              likcand2D(l)+=ll_y_cand(i,j,l);
            }
          }
        }
        rand2=Rcpp::runif(1);
        explldiff = exp(likcand2D(l)-likcurr2D(l));
        if(rand2(0)<explldiff){
          lam0(l)=lam0cand(l);
          ll_y_curr.slice(l)=ll_y_cand.slice(l);
          pd.slice(l)=pdcand.slice(l);
          lamd.slice(l)=lamdcand.slice(l);
          likcurr2D(l)=likcand2D(l);
        }
      }
      likcurr+=likcurr2D(l);
    }
  }
  //fill lam0use
  for(int i=0; i<t; i++) {
    if(lam0.size()==1){
      lam0use(i)=lam0(0);
    }else{
      lam0use(i)=lam0(i);
    }
  }
  // Update sigma
  if(sigma.size()==1){//fixed sigma
    rand=Rcpp::rnorm(1,sigma(0),0.025);
    if(rand(0) > 0){
      sigmacand(0)=rand(0);
      likcand=0;
      //  Update lamd and calculate cand likelihood
      for(int l=0; l<t; l++){
        likcand2D(l)=0;
        for(int i=0; i<M; i++) {
          for(int j=0; j<Xidx(l); j++){
            lamdcand(i,j,l)=lam0use(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmacand(0)*sigmacand(0)));
            pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
            ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
            if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
              likcand2D(l)+=ll_y_cand(i,j,l);
            }
          }
        }
        likcand+=likcand2D(l);
      }
      rand2=Rcpp::runif(1);
      explldiff = exp(likcand-likcurr);
      if(rand2(0)<explldiff){
        sigma(0)=sigmacand(0);
        lamd=lamdcand;
        pd=pdcand;
        ll_y_curr=ll_y_cand;
        likcurr=likcand;
        likcurr2D=likcand2D;
      }
    }
  }else{//year-specific sigmas
    likcurr=0;
    for(int l=0; l<t; l++){
      rand=Rcpp::rnorm(1,sigma(l),0.025);
      likcand2D(l)=0;
      if(rand(0) > 0){
        sigmacand(l)=rand(0);
        //  Calculate likelihood
        for(int i=0; i<M; i++) {
          for(int j=0; j<Xidx(l); j++){
            lamdcand(i,j,l)=lam0use(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmacand(l)*sigmacand(l)));
            pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
            ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
            if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
              likcand2D(l)+=ll_y_cand(i,j,l);
            }
          }
        }
        rand2=Rcpp::runif(1);
        explldiff = exp(likcand2D(l)-likcurr2D(l));
        if(rand2(0)<explldiff){
          sigma(l)=sigmacand(l);
          ll_y_curr.slice(l)=ll_y_cand.slice(l);
          pd.slice(l)=pdcand.slice(l);
          lamd.slice(l)=lamdcand.slice(l);
          likcurr2D(l)=likcand2D(l);
        }
      }
      likcurr+=likcurr2D(l);
    }
  }
  List to_return(4);
  to_return[0] = lam0;
  to_return[1] = sigma;
  to_return[2] = lamd;
  to_return[3] = pd;
  return to_return;
}


//Update psi,z, gamma and phi
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::export]]
List upZ(IntegerMatrix z, IntegerMatrix a,IntegerMatrix knownmatrix,IntegerVector Xidx,
         arma::cube y, arma::cube pd, IntegerVector K,NumericMatrix Ez,
         NumericVector gamma,NumericVector gammaprime, NumericVector phi, double psi,IntegerVector N,int maxJ,
         NumericVector propz, double propgamma) {
  RNGScope scope;
  //Preallocate
  int M = z.nrow();
  int t = z.ncol();
  NumericVector z1curr(M);
  NumericVector z1cand(M);
  NumericVector z1tmp(M);
  NumericVector a1tmp(M);
  NumericVector gammaprimecand(t);
  arma::cube ll_y_curr(M,maxJ,t);
  NumericMatrix ll_z(M,t);
  NumericVector rand;
  NumericMatrix Ezcand(M,t);
  NumericMatrix ll_z_cand(M,t);
  arma::cube ll_y_cand(M,maxJ,t);
  double llysum=0;
  double llycandsum=0;
  double llzsum=0;
  double llzcandsum=0;
  int sumz;
  int sumz1tmp=0;
  int suma1tmp=0;
  int suma=0;
  LogicalVector warn(M,FALSE);
  //Preallocate z[,2+]
  NumericVector pr_zt(M);
  IntegerVector at_cand(M);
  NumericVector zt_cand(M);
  double prop_probs;
  double back_probs;
  int navail=0;
  int idx=0;
  //Preallocate phi and gamma
  int survive=0;
  int dead=0;
  bool gamma_cand_ok=TRUE;
  double gamma_cand=0;
  NumericVector gamma_prime_cand(t-1);
  //Fixed or year-specific stuff
  NumericVector gammacand(t);
  NumericVector phicand(t);
  NumericVector gammause(t);
  NumericVector phiuse(t);
  for(int i=0; i<t; i++) {
    if(gamma.size()==1){
      gammause(i)=gamma(0);
    }else{
      gammause(i)=gamma(i);
    }
  }
  for(int i=0; i<t; i++) {
    if(phi.size()==1){
      phiuse(i)=phi(0);
    }else{
      phiuse(i)=phi(i);
    }
  }
  ////////////////////Z1 stuff////////////////////
  //Figure out who can be updated
  LogicalVector upz(M);
  IntegerVector latecaps(M,0);
  if(t==2){
    for(int i=0; i<M; i++){
      upz(i)=(knownmatrix(i,0)==0);
    }
  }else{
    for(int i=0; i<M; i++){
      for(int j=2; j<t; j++){
        latecaps(i)+=z(i,j);
      }
      upz(i)=(!((z(i,0)==0)&(z(i,1)==0)&(latecaps(i)>0)))&(knownmatrix(i,0)==0);
    }
  }
  IntegerVector allon(t);
  IntegerVector alloff(t);
  for(int i=0; i<t; i++){
    allon(i)=1;
    alloff(i)=0;
  }
  //Calculate ll.z. Can move outside MCMC later.
  //z1
  for(int i=0; i<M; i++){
    ll_z(i,0)= z(i,0)*log(psi)+(1-z(i,0))*log(1-psi);
  }
  //z2+
  for(int l=1; l<t; l++){
    for(int i=0; i<M; i++){
      ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
      if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
        ll_z(i,l)=0;
      }
    }
  }
  //same for ll.y
  for(int i=0; i<M; i++){
    for(int j=0; j<Xidx(0); j++){
      for(int l=0; l<t; l++){
        ll_y_curr(i,j,l) =z(i,l)*(y(i,j,l)*log(pd(i,j,l))+(K(l)-y(i,j,l))*log(1-pd(i,j,l)));
        ll_y_cand(i,j,l) =ll_y_curr(i,j,l);
      }
    }
  }

  //update z[,1]
  N(0)=0;
  for(int i=0; i<M; i++){
    if(upz(i)){
      z1curr = z(_,0);
      z1cand = 1-z1curr;
      gammaprimecand = gammaprime;
      z1tmp = z1curr;
      z1tmp(i) = z1cand(i);
      a1tmp = a(_,0);
      a1tmp(i) = 1-z1cand(i);
      sumz1tmp=0;
      suma1tmp=0;
      for(int i2=0; i2<M; i2++){
        sumz1tmp+=z1tmp(i2);
        suma1tmp+=a1tmp(i2);
      }
      gammaprimecand(0)=sumz1tmp*gammause(0)/suma1tmp;
      if(gammaprimecand(1) > 1) { // E(Recruits) must be < nAvailable
        warn(i)=TRUE;
      }
      if(warn(i)==FALSE){
        Ezcand=clone(Ez);
        ll_z_cand=clone(ll_z);
        ll_z_cand(i,0) = z1cand(i)*log(psi)+(1-z1cand(i))*log(1-psi);
        for(int i2=0; i2<M; i2++){//update Ezcand and ll_z_cand for everyone due to gammaprimecand update
          Ezcand(i2,0)=z(i2,0)*phiuse(0) + a(i2,0)*gammaprimecand(0);
          ll_z_cand(i2,1) = z(i2,1)*log(Ezcand(i2,0))+(1-z(i2,1))*log(1-Ezcand(i2,0));
        }
        //Now update for focal individual
        Ezcand(i,0) = z1cand(i)*phiuse(0) + (1-z1cand(i))*gammaprimecand(0);
        ll_z_cand(i,1) = z(i,1)*log(Ezcand(i,0))+(1-z(i,1))*log(1-Ezcand(i,0));
        //sum z ll
        llzcandsum=ll_z_cand(i,0);//add ll.z[i,1] for focal guy
        llzsum=ll_z(i,0);
        //then add ll.z[,2] for all guys
        for(int i2=0; i2<M; i2++){
          if(ll_z_cand(i2,1)==ll_z_cand(i2,1)){
            llzcandsum+=ll_z_cand(i2,1);
          }
          if(ll_z(i2,1)==ll_z(i2,1)){
            llzsum+=ll_z(i2,1);
          }
        }
        //sum y ll
        llycandsum=0;
        llysum=0;
        for(int j=0; j<Xidx(0); j++){
          ll_y_cand(i,j,0) =z1cand(i)*(y(i,j,0)*log(pd(i,j,0))+(K(0)-y(i,j,0))*log(1-pd(i,j,0)));
          if(ll_y_cand(i,j,0)==ll_y_cand(i,j,0)){
            llycandsum+=ll_y_cand(i,j,0);
          }
          ll_y_curr(i,j,0) =z(i,0)*(y(i,j,0)*log(pd(i,j,0))+(K(0)-y(i,j,0))*log(1-pd(i,j,0)));
          if(ll_y_curr(i,j,0)==ll_y_curr(i,j,0)){
            llysum+=ll_y_curr(i,j,0);
          }
        }
        rand=Rcpp::runif(1);
        if(rand(0) < exp((llycandsum+ llzcandsum)-(llysum+llzsum ))) {
          ll_y_curr.subcube(i,0,0,i,maxJ-1,0) = ll_y_cand.subcube(i,0,0,i,maxJ-1,0);
          Ez =clone(Ezcand);//cloning whole thing prob not most efficient. need first 2 dim
          ll_z = clone(ll_z_cand);//same here
          z(i,0) = z1cand(i);
          sumz=0;
          for(int j=0; j<t; j++){
            sumz+= z(i,j);
          }
          gammaprime(0) = gammaprimecand(0);
          if(z1cand(i)==1){
            a(i,_)=alloff;//if caught occ 1, never available to recruit
          }else if (sumz==0){  //if never caught, turn availability all on
            a(i,_)=allon;
          }else {  //if not caught on occ1, but caught later, available to recruit in occ2
            a(i,0)=1;
          }
        }
      }
    }
    N(0)+=z(i,0);
  }
  // Finally, need to update gamma.prime[,2:t], Ez[,2:t-1], and ll.z[,3:t]
  // for z1 guys now caught in first year, previously not captured (all a turned off) and guys now never caught
  // but previously only captured on occasion 1 (all a turned on)
  // But only if t>2
if(t>2){
  for(int l=2; l<t; l++){
    suma=0;
    for(int i=0; i<M; i++){
      suma+=a(i,l-1);
    }
    gammaprime(l-1)=(N(l-1)*gammause(l-1)) / suma;
    for(int i=0; i<M; i++){
      Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
      ll_z(i,l)=z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
      if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
        ll_z(i,l)=0;
      }
    }
  }
}
//update z[,2+]
for(int l=1; l<t; l++){
  //figure out who can be updated with upz
  //always remove dead guys
  if(t==2){
    for(int i=0; i<M; i++){
      upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));
    }
  }else if(t==3){
    if(l==1){
      for(int i=0; i<M; i++){
        upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));  //remove guys that are on before and after l
      }
    }else{
      for(int i=0; i<M; i++){
        upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));
      }
    }
  }else{//t>3
    if(l==(t-1)){ //can update anyone on last occasion unless they're dead
      for(int i=0; i<M; i++){
        upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));
      }
    }else if(l==(t-2)){ //second to last occasion
      for(int i=0; i<M; i++){
        upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));  //remove guys that are on before and after l
      }
    }else{//l between 2 and t-2
      for(int i=0; i<M; i++){
        latecaps(i)=0;
        for(int j=(l+2); j<t; j++){
          latecaps(i)+=z(i,j);
        }
        upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!((z(i,l)==0)&(z(i,l+1)==0)&(latecaps(i)>0)))&(!( (z(i,l-1)==0) & (a(i,l-1)==0))); //remove guys that are on before and after l or 0,0,1
      }
    }
  }
  //Can we swap?
  navail=0;
  for(int i=0; i<M; i++){
    navail+=upz(i);
  }
  int propzuse=0;
  //How many to swap?
  if(navail < propz(l-1)) {
    propzuse=navail;
  }else{
    propzuse=propz(l-1);
  }
  //get upz2
  IntegerVector upz2(navail,0);
  idx=0;
  for(int i=0; i<M; i++){
    if(upz(i)){
      upz2(idx)=i;
      idx+=1;
    }
  }
  if(navail>0){//if so, proceed with swapping
    //Who to swap? structure size varies
    IntegerVector swapzidx(propzuse,0);
    IntegerVector choose=Rcpp::seq(0,(navail-1));
    swapzidx=Rcpp::RcppArmadillo::sample(choose,propzuse,FALSE);
    IntegerVector swapz(propzuse,0);
    for(int i=0; i<propzuse; i++){
      swapz(i)=upz2(swapzidx(i));
    }
    for(int i=0; i<M; i++){
      pr_zt(i)=z(i,l);
      zt_cand(i)=z(i,l);
    }
    for(int i=0; i<propzuse; i++){
      pr_zt(swapz(i))=Ez(swapz(i),l-1);
      rand=Rcpp::rbinom(1,1,pr_zt(swapz(i)));
      zt_cand(swapz(i))=rand(0);
    }
    at_cand=1*((a(_,l-1)==1)&(zt_cand==0)); //who was available on last occasion and not proposed to be captured?
    prop_probs=0;
    back_probs=0;
    ll_z_cand=clone(ll_z);
    //These are different ll sums. Only for updated individuals
    llycandsum=0;
    llysum=0;
    llzcandsum=0;
    llzsum=0;
    for(int i=0; i<propzuse; i++){
      prop_probs+=zt_cand(swapz(i))*log(pr_zt(swapz(i)))+(1-zt_cand(swapz(i)))*log(1-pr_zt(swapz(i)));
      back_probs+=z(swapz(i),l)*log(pr_zt(swapz(i)))+(1-z(swapz(i),l))*log(1-pr_zt(swapz(i)));
      for(int j=0; j<Xidx(l); j++){
        ll_y_cand(swapz(i),j,l)=zt_cand(swapz(i))*(y(swapz(i),j,l)*log(pd(swapz(i),j,l))+(K(0)-y(swapz(i),j,l))*log(1-pd(swapz(i),j,l)));
        //Add up y likelihood contributions curr and cand for swapped guys
        if(ll_y_cand(swapz(i),j,l)==ll_y_cand(swapz(i),j,l)){
          llycandsum+=ll_y_cand(swapz(i),j,l);
        }
        if(ll_y_curr(swapz(i),j,l)==ll_y_curr(swapz(i),j,l)){
          llysum+=ll_y_curr(swapz(i),j,l);
        }
      }
      //Add up the z likelihood contributions curr and cand for swapped guys
      ll_z_cand(swapz(i),l) = zt_cand(swapz(i))*log(Ezcand(swapz(i),l-1))+(1-zt_cand(swapz(i)))*log(1-Ezcand(swapz(i),l-1));
      if(ll_z(swapz(i),l)!=ll_z(swapz(i),l)){//Turn NaNs to 0s
        ll_z(swapz(i),l)=0;
      }
      llzcandsum+=ll_z_cand(swapz(i),l);//prior.z.cand
      llzsum+=ll_z(swapz(i),l);//prior.z
    }
    //If we make this change to z[,l] and a[,l], how does it change ll.z[,l+1]?
    if(l<(t-1)){ //extra stuff if this is true
      sumz1tmp=0;
      suma1tmp=0;
      for(int i=0; i<M; i++){
        sumz1tmp+=zt_cand(i);
        suma1tmp+=at_cand(i);
      }
      gammaprimecand(l) = sumz1tmp*gammause(l) / suma1tmp;
      if(gammaprimecand(l) <= 1){ //is gamma prime permissible?
        for(int i=0; i<M; i++){
          //Add on contributions to z ll from l+1 for cand and curr
          Ezcand(i,l) = zt_cand(i)*phiuse(l) + at_cand(i)*gammaprimecand(l);
          ll_z_cand(i,l+1) = z(i,l+1)*log(Ezcand(i,l))+(1-z(i,l+1))*log(1-Ezcand(i,l));
          if(ll_z(i,l+1)==ll_z(i,l+1)){
            llzsum+=ll_z(i,l+1);//prior.z
          }
          if(ll_z_cand(i,l+1)==ll_z_cand(i,l+1)){
            llzcandsum+=ll_z_cand(i,l+1);//prior.z.cand
          }
        }
        rand=Rcpp::runif(1);
        if(rand(0) < exp(( llycandsum+llzcandsum+back_probs)-( llysum+llzsum+prop_probs) )) {
          for(int i=0; i<propzuse; i++){
            for(int j=0; j<Xidx(l); j++){
              ll_y_curr(swapz(i),j,l) = ll_y_cand(swapz(i),j,l);
            }
          }
          ll_z(_,l)=ll_z_cand(_,l);
          z(_,l)=zt_cand;
          a(_,l)=at_cand;
          //Now update t<l stuff
          ll_z(_,l+1)=ll_z_cand(_,l+1);
          Ez(_,l)=Ezcand(_,l);
          gammaprime(l)=gammaprimecand(l);
          if(t>3){//more houskeeping if t>3. if t=3 you can update a[,3] but it doesn't matter. can't recruit after survey is over
            // for the same reason, you don't need to update if (l+1)==t and t>3. I exclude those below so a is updated
            // for the last occasion and my sanity check works.
            // turn off availability for (l+1):t if you enter population but weren't there before
            LogicalVector fix1(swapz.size());
            bool fix3=FALSE;
            for(int i=0; i<propzuse; i++){
              suma=0;
              for(int l2=l; l2<t; l2++){
                suma+=a(swapz(i),l2);
              }
              if((zt_cand(swapz(i))==1)&(suma>0)){
                fix1(i)=TRUE;
                for(int l2=l; l2<t; l2++){
                  a(swapz(i),l2)=0;
                }
                fix3=TRUE;
              }else{
                fix1(i)=FALSE;
              }
            }
            // turn on availability for (l+1):t if you leave the population
            LogicalVector fix2(M,swapz.size());
            for(int i=0; i<propzuse; i++){
              suma=0;
              for(int l2=0; l2<t; l2++){
                suma+=a(swapz(i),l2);
              }
              sumz=0;
              for(int l2=0; l2<t; l2++){
                sumz+=z(swapz(i),l2);
              }
              if((sumz==0)&(suma!=t)){
                fix2(i)=TRUE;
                for(int l2=l; l2<t; l2++){
                  a(swapz(i),l2)=1;
                }
                fix3=TRUE;
              }else{
                fix2(i)=FALSE;
              }
            }
            // If we updated anything here,a[,(l+1):t] has changed and
            // we have to recalculate gamma.prime[,(l+1):(t-1)], Ez.cand[,(l+1):(t-1)], and ll.z.cand[,(l+2):t]
            if(((l+1)!=(t-1))){//t-1 instead of t for Rcpp
              if(fix3){
                for(int l2=l+1; l2<(t-1); l2++){
                  sumz1tmp=0;
                  suma1tmp=0;
                  for(int i=0; i<M; i++){
                    sumz1tmp+=z(i,l2);
                    suma1tmp+=a(i,l2);
                  }
                  gammaprime(l2) = sumz1tmp*gammause(l2) / suma1tmp;
                  Ez(_,l2)=z(_,l2)*phiuse(l2) + a(_,l2)*gammaprime(l2);
                  for(int i=0; i<M; i++){
                    //Add on contributions to z ll from l+1 for cand and curr
                    Ez(i,l2) = z(i,l2)*phiuse(l2) + a(i,l2)*gammaprime(l2);
                    ll_z(i,l2+1) = z(i,l2+1)*log(Ez(i,l2))+(1-z(i,l2+1))*log(1-Ez(i,l2));
                  }
                }
              }
            }
          }
        }
      }
    }else{//if l==t don't need to consider next time step
      rand=Rcpp::runif(1);
      if(rand(0) < exp(( llycandsum+llzcandsum+back_probs)-( llysum+llzsum+prop_probs) )) {
        for(int i=0; i<propzuse; i++){
          for(int j=0; j<Xidx(l); j++){
            ll_y_curr(swapz(i),j,l) = ll_y_cand(swapz(i),j,l);
          }
        }
        ll_z(_,l)=ll_z_cand(_,l);
        z(_,l)=zt_cand;
        a(_,l)=at_cand;
      }
    }
    N(l)=0;
    for(int i=0; i<M; i++){
      N(l)+=z(i,l);
    }
    //Update psi
    rand=Rcpp::rbeta(1, 1+N(0), 1+M-N(1));
    psi=rand(0);
    for(int i=0; i<M; i++){
      ll_z(i,0)= z(i,0)*log(psi)+(1-z(i,0))*log(1-psi);
    }

  }
}
//Update phi
if(phi.size()==(t-1)){//if time-specific survival
  for(int l=1; l<t; l++){
    survive=0;
    dead=0;
    for(int i=0; i<M; i++){
      survive+=(z(i,l-1)==1)&(z(i,l)==1);
      dead+=(z(i,l-1)==1)&(z(i,l)==0);
    }
    rand=Rcpp::rbeta(1, 1+survive, 1+dead);
    phi(l-1)=rand(0);
    phiuse(l-1)=rand(0);
  }
}else{
  survive=0;
  dead=0;
  for(int l=1; l<t; l++){
    for(int i=0; i<M; i++){
      survive+=(z(i,l-1)==1)&(z(i,l)==1);
      dead+=(z(i,l-1)==1)&(z(i,l)==0);
    }
  }
  rand=Rcpp::rbeta(1, 1+survive, 1+dead);
  phi(0)=rand(0);
  for(int l=1; l<(t+1); l++){
    phiuse(l-1)=rand(0);
  }
}
//  Update gamma
// NOTE: Must update ll.z, Ez, etc...
if(gamma.size()==1){
  rand=Rcpp::rnorm(1, gamma(0), propgamma);
  gamma_cand=rand(0);
  gamma_cand_ok=TRUE;
  for(int l=1; l<t; l++){
    suma=0;
    for(int i=0; i<M; i++){
      suma+=a(i,l-1);
    }
    gamma_prime_cand(l-1)=(N(l-1)*gamma_cand) / suma;
    if(gamma_prime_cand(l-1) > 1){  //Note don't break loop b/c ll.z needs updating because phi changed
      gamma_cand_ok=FALSE;
    }
    llzsum=0;
    for(int i=0; i<M; i++){
      Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
      ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
      if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
        ll_z(i,l)=0;
      }
      llzsum+=ll_z(i,l);
    }
  }
  if((gamma_cand>0)&(gamma_cand_ok)){
    llzcandsum=0;
    for(int l=1; l<t; l++){
      for(int i=0; i<M; i++){
        Ezcand(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gamma_prime_cand(l-1);
        ll_z_cand(i,l)= z(i,l)*log(Ezcand(i,l-1))+(1-z(i,l))*log(1-Ezcand(i,l-1));
        if(ll_z_cand(i,l)!=ll_z_cand(i,l)){//Turn NaNs to 0s
          ll_z_cand(i,l)=0;
        }
        llzcandsum+=ll_z_cand(i,l);
      }
    }
    rand=Rcpp::runif(1);
    if(rand(0) < exp(llzcandsum - llzsum)){
      gamma(0)=gamma_cand;
      for(int l=1; l<t; l++){
        gammaprime(l-1)=gamma_prime_cand(l-1);
        gammause(l-1)=gamma_prime_cand(l-1);
      }
      Ez=clone(Ezcand);
      ll_z=clone(ll_z_cand);
    }
  }
}else{
  for(int l=1; l<t; l++){
    rand=Rcpp::rnorm(1, gamma(l-1), propgamma);
    gamma_cand=rand(0);
    gamma_cand_ok=TRUE;
    suma=0;
    for(int i=0; i<M; i++){
      suma+=a(i,l-1);
    }
    gamma_prime_cand(l-1)=(N(l-1)*gamma_cand) / suma;
    if(gamma_prime_cand(l-1) > 1){  //Note don't break loop b/c ll.z needs updating because phi changed
      gamma_cand_ok=FALSE;
    }
    llzsum=0;
    for(int i=0; i<M; i++){
      Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
      ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
      if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
        ll_z(i,l)=0;
      }
      llzsum+=ll_z(i,l);
    }
    if((gamma_cand>0)&(gamma_cand_ok)){
      llzcandsum=0;
      for(int i=0; i<M; i++){
        Ezcand(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gamma_prime_cand(l-1);
        ll_z_cand(i,l)= z(i,l)*log(Ezcand(i,l-1))+(1-z(i,l))*log(1-Ezcand(i,l-1));
        if(ll_z_cand(i,l)!=ll_z_cand(i,l)){//Turn NaNs to 0s
          ll_z_cand(i,l)=0;
        }
        llzcandsum+=ll_z_cand(i,l);
      }
      rand=Rcpp::runif(1);
      if(rand(0) < exp(llzcandsum - llzsum)){
        gamma(l-1)=gamma_cand;
        gammaprime(l-1)=gamma_prime_cand(l-1);
        gammause(l-1)=gamma_prime_cand(l-1);
        for(int i=0; i<M; i++){
          Ez(i,l-1)=Ezcand(i,l-1);
          ll_z(i,l-1)=ll_z_cand(i,l-1);
        }
      }
    }
  }
}


List to_return(12);
to_return[0] = ll_z;
to_return[1] = Ez;
to_return[2] = gammaprime;
to_return[3] = z;
to_return[4] = a;
to_return[5] = N;
to_return[6] = psi;
to_return[7] = phi;
to_return[8] = gamma;
to_return[9] = gammause;
to_return[10] = gammaprime;
to_return[11] = phiuse;
return to_return;
}

//Full open pop MCMC sampler
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::export]]
List mcmc_Open(NumericVector lam0, NumericVector sigma, NumericVector gamma,NumericVector gammaprime, NumericVector phi,
               arma::cube D,arma::cube lamd, arma::cube y,IntegerMatrix z,IntegerMatrix a, NumericMatrix s1,arma::cube s2,
               bool metamu, bool useverts,NumericMatrix vertices,NumericVector xlim,NumericVector ylim,
               IntegerMatrix knownmatrix,NumericVector Xidx, arma::cube Xcpp,NumericVector K,NumericMatrix Ez, double psi,
               IntegerVector N,NumericVector proplam0, NumericVector propsig,NumericVector propz, NumericVector propgamma,double props1x,
               double props1y,double props2x,double props2y, double propsigma_t,NumericVector sigma_t,
               int niter, int nburn, int nthin,int npar,NumericVector each) {
  RNGScope scope;
  //Preallocate detection function
  int M = size(lamd)[0];
  int J = size(lamd)[1];
  int t = size(lamd)[2];
  NumericVector lam0cand(t);
  NumericVector sigmacand(t);
  NumericVector rand;
  NumericVector rand2;
  arma::cube pd=zeros<cube>(M,J,t);
  arma::cube ll_y_curr=zeros<cube>(M,J,t);
  arma::cube lamdcand=zeros<cube>(M,J,t);
  arma::cube pdcand=zeros<cube>(M,J,t);
  arma::cube ll_y_cand=zeros<cube>(M,J,t);
  //Structures for dealing with year specific parms
  NumericVector likcurr2D(t);
  NumericVector likcand2D(t);
  NumericVector sigmause(t);
  NumericVector lam0use(t);
  for(int i=0; i<t; i++) {
    if(sigma.size()==1){
      sigmause(i)=sigma(0);
    }else{
      sigmause(i)=sigma(i);
    }
  }
  for(int i=0; i<t; i++) {
    if(lam0.size()==1){
      lam0use(i)=lam0(0);
    }else{
      lam0use(i)=lam0(i);
    }
  }
  //Preallocate Z update
  NumericVector z1curr(M);
  NumericVector z1cand(M);
  NumericVector z1tmp(M);
  NumericVector a1tmp(M);
  NumericVector gammaprimecand(t);
  NumericMatrix ll_z(M,t);
  NumericMatrix Ezcand(M,t);
  NumericMatrix ll_z_cand(M,t);
  double llysum=0;
  double llycandsum=0;
  double llzsum=0;
  double llzcandsum=0;
  int sumz;
  int sumz1tmp=0;
  int suma1tmp=0;
  int suma=0;
  LogicalVector warn(M,FALSE);  ////fix this!!
  LogicalVector upz(M);
  IntegerVector latecaps(M,0);
  IntegerVector allon(t);
  IntegerVector alloff(t);
  for(int i=0; i<t; i++){
    allon(i)=1;
    alloff(i)=0;
  }
  //Preallocate z[,2+]
  NumericVector pr_zt(M);
  IntegerVector at_cand(M);
  NumericVector zt_cand(M);
  double prop_probs;
  double back_probs;
  int navail=0;
  int idx=0;
  //Preallocate phi and gamma
  int survive=0;
  int dead=0;
  bool gamma_cand_ok=TRUE;
  double gamma_cand=0;
  NumericVector gamma_prime_cand(t-1);
  //Fixed or year-specific stuff
  NumericVector gammacand(t-1);
  NumericVector phicand(t-1);
  NumericVector gammause(t-1);
  NumericVector phiuse(t-1);
  for(int l=0; l<(t-1); l++) {
    if(gamma.size()==1){
      gammause(l)=gamma(0);
    }else{
      gammause(l)=gamma(l);
    }
  }
  for(int l=0; l<(t-1); l++) {
    if(phi.size()==1){
      phiuse(l)=phi(0);
    }else{
      phiuse(l)=phi(l);
    }
  }
  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericMatrix dtmp(J,t);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  int zany=FALSE;
  //metamu stuff
  NumericMatrix ll_s2(M,t);
  NumericMatrix ll_s2_cand(M,t);
  double lls2sum=0;
  double lls2candsum=0;
  NumericVector sigma_t_cand(1);
  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,npar);
  NumericMatrix s1xout(nstore,M);
  NumericMatrix s1yout(nstore,M);
  arma::cube s2xout(nstore,M,t);
  arma::cube s2yout(nstore,M,t);
  arma::cube zout(nstore,M,t);
  int iteridx=0;
  //////Calculate starting log likelihoods///////
  //ll.s2
  if(metamu){
    for(int l=0; l<t; l++){
      for(int i=0; i<M; i++) { //X and Y normal log-likelihood simplified
        ll_s2(i,l)=pow(sigma_t(0),2.0)-(1/(2*pow(sigma_t(0),2.0)))*(pow(s2(i,l,0)-s1(i,0),2.0)+pow(s2(i,l,1)-s1(i,1),2.0));
      }
    }
  }
  //  Detection function
  for(int l=0; l<t; l++){
    likcurr2D(l)=0;
    for(int i=0; i<M; i++) {
      for(int j=0; j<Xidx(l); j++){
        pd(i,j,l)=1-exp(-lamd(i,j,l));
        ll_y_curr(i,j,l)=z(i,l)*(y(i,j,l)*log(pd(i,j,l))+(K[l]-y(i,j,l))*log(1-pd(i,j,l)));
        if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
          likcurr2D(l)+=ll_y_curr(i,j,l); //year-specific components
        }
      }
    }
    llysum+=likcurr2D(l); //full likelihood sum
  }
  //ll.z.
  for(int i=0; i<M; i++){//z1
    ll_z(i,0)= z(i,0)*log(psi)+(1-z(i,0))*log(1-psi);
  }
  for(int l=1; l<t; l++){//z2+
    for(int i=0; i<M; i++){
      ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
      if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
        ll_z(i,l)=0;
      }
    }
  }

//Here we go!
  for(int iter=0; iter<niter; iter++){
    /////////////////Detection function update///////////////////////
    //Need to resum the ll_y on each iter
    llysum=0;
    for(int l=0; l<t; l++){
      likcurr2D(l)=0;
      for(int i=0; i<M; i++) {
        for(int j=0; j<Xidx(l); j++){
          if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
            likcurr2D(l)+=ll_y_curr(i,j,l);
          }
        }
      }
      llysum+=likcurr2D(l); //full likelihood sum
    }
    //Update lam0
    if(lam0.size()==1){//fixed lambda
      rand=Rcpp::rnorm(1,lam0(0),proplam0(0));
      if(rand(0) > 0){
        llycandsum=0;
        lam0cand(0)=rand(0);
        //  Update lamd and calculate cand likelihood
        for(int l=0; l<t; l++){
          likcand2D(l)=0;
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0cand(0)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmause(l)*sigmause(l)));
              pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
              ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          llycandsum+=likcand2D(l);
        }
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(llycandsum-llysum)){
          lam0(0)=lam0cand(0);
          lamd=lamdcand;
          pd=pdcand;
          ll_y_curr=ll_y_cand;
          llysum=llycandsum;
          likcurr2D=likcand2D;
        }
      }
    }else{//year-specific lambdas
      llysum=0;
      for(int l=0; l<t; l++){
        likcand2D(l)=0;
        rand=Rcpp::rnorm(1,lam0(l),proplam0(l));
        if(rand(0) > 0){
          lam0cand(l)=rand(0);
          //  Calculate likelihood
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0cand(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmause(l)*sigmause(l)));
              pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
              ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          rand2=Rcpp::runif(1);
          if(rand2(0)<exp(likcand2D(l)-likcurr2D(l))){
            lam0(l)=lam0cand(l);
            ll_y_curr.slice(l)=ll_y_cand.slice(l);
            pd.slice(l)=pdcand.slice(l);
            lamd.slice(l)=lamdcand.slice(l);
            likcurr2D(l)=likcand2D(l);
          }
        }
        llysum+=likcurr2D(l);
      }
    }
    //fill lam0use
    for(int i=0; i<t; i++) {
      if(lam0.size()==1){
        lam0use(i)=lam0(0);
      }else{
        lam0use(i)=lam0(i);
      }
    }
    // Update sigma
    if(sigma.size()==1){//fixed sigma
      rand=Rcpp::rnorm(1,sigma(0),propsig(0));
      if(rand(0) > 0){
        sigmacand(0)=rand(0);
        llycandsum=0;
        //  Update lamd and calculate cand likelihood
        for(int l=0; l<t; l++){
          likcand2D(l)=0;
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0use(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmacand(0)*sigmacand(0)));
              pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
              ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          llycandsum+=likcand2D(l);
        }
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(llycandsum-llysum)){
          sigma(0)=sigmacand(0);
          lamd=lamdcand;
          pd=pdcand;
          ll_y_curr=ll_y_cand;
          llysum=llycandsum;
          likcurr2D=likcand2D;
        }
      }
    }else{//year-specific sigmas
      llysum=0;
      for(int l=0; l<t; l++){
        rand=Rcpp::rnorm(1,sigma(l),propsig(l));
        likcand2D(l)=0;
        if(rand(0) > 0){
          sigmacand(l)=rand(0);
          //  Calculate likelihood
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0use(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmacand(l)*sigmacand(l)));
              pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
              ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          rand2=Rcpp::runif(1);
          if(rand2(0)<exp(likcand2D(l)-likcurr2D(l))){
            sigma(l)=sigmacand(l);
            ll_y_curr.slice(l)=ll_y_cand.slice(l);
            pd.slice(l)=pdcand.slice(l);
            lamd.slice(l)=lamdcand.slice(l);
            likcurr2D(l)=likcand2D(l);
          }
        }
        llysum+=likcurr2D(l);
      }
    }
    //fill sigmause
    for(int l=0; l<t; l++) {
      if(sigma.size()==1){
        sigmause(l)=sigma(0);
      }else{
        sigmause(l)=sigma(l);
      }
    }
    ////////////////////Z1 stuff////////////////////
    // Figure out who can be updated
    if(t==2){
      for(int i=0; i<M; i++){
        upz(i)=(knownmatrix(i,0)==0);
      }
    }else{
      for(int i=0; i<M; i++){
        for(int j=2; j<t; j++){
          latecaps(i)+=z(i,j);
        }
        upz(i)=(!((z(i,0)==0)&(z(i,1)==0)&(latecaps(i)>0)))&(knownmatrix(i,0)==0);
      }
    }
    //update z[,1]
    N(0)=0;
    for(int i=0; i<M; i++){
      ll_z(i,0)= z(i,0)*log(psi)+(1-z(i,0))*log(1-psi);
    }
    //z2+
    for(int l=1; l<t; l++){
      for(int i=0; i<M; i++){
        ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
        if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
          ll_z(i,l)=0;
        }
      }
    }
    //same for ll.y
    for(int i=0; i<M; i++){
      for(int j=0; j<Xidx(0); j++){
        for(int l=0; l<t; l++){
          ll_y_curr(i,j,l) =z(i,l)*(y(i,j,l)*log(pd(i,j,l))+(K(l)-y(i,j,l))*log(1-pd(i,j,l)));
          ll_y_cand(i,j,l) =ll_y_curr(i,j,l);
        }
      }
    }
    for(int i=0; i<M; i++){
      if(upz(i)){
        z1curr = z(_,0);
        z1cand = 1-z1curr;
        gammaprimecand = clone(gammaprime);
        z1tmp = z1curr;
        z1tmp(i) = z1cand(i);
        a1tmp = a(_,0);
        a1tmp(i) = 1-z1cand(i);
        sumz1tmp=0;
        suma1tmp=0;
        for(int i2=0; i2<M; i2++){
          sumz1tmp+=z1tmp(i2);
          suma1tmp+=a1tmp(i2);
        }
        gammaprimecand(0)=sumz1tmp*gammause(0)/suma1tmp;
        if(gammaprimecand(0) > 1) { // E(Recruits) must be < nAvailable
          warn(i)=TRUE;
        }
        if(warn(i)==FALSE){
          Ezcand=clone(Ez);
          ll_z_cand=clone(ll_z);
          ll_z_cand(i,0) = z1cand(i)*log(psi)+(1-z1cand(i))*log(1-psi);
          for(int i2=0; i2<M; i2++){//update Ezcand and ll_z_cand for everyone due to gammaprimecand update
            Ezcand(i2,0)=z(i2,0)*phiuse(0) + a(i2,0)*gammaprimecand(0);
            ll_z_cand(i2,1) = z(i2,1)*log(Ezcand(i2,0))+(1-z(i2,1))*log(1-Ezcand(i2,0));
          }
          //Now update for focal individual
          Ezcand(i,0) = z1cand(i)*phiuse(0) + (1-z1cand(i))*gammaprimecand(0);
          ll_z_cand(i,1) = z(i,1)*log(Ezcand(i,0))+(1-z(i,1))*log(1-Ezcand(i,0));
          //sum z ll
          llzcandsum=ll_z_cand(i,0);//add ll.z[i,1] for focal guy
          llzsum=ll_z(i,0);
          //then add ll.z[,2] for all guys
          for(int i2=0; i2<M; i2++){
            if(ll_z_cand(i2,1)==ll_z_cand(i2,1)){
              llzcandsum+=ll_z_cand(i2,1);
            }
            if(ll_z(i2,1)==ll_z(i2,1)){
              llzsum+=ll_z(i2,1);
            }
          }
          //sum y ll across j dimension for each i and l=0
          llysum=0;
          llycandsum=0;
          for(int j=0; j<Xidx(0); j++){
            ll_y_cand(i,j,0) =z1cand(i)*(y(i,j,0)*log(pd(i,j,0))+(K(0)-y(i,j,0))*log(1-pd(i,j,0)));
            if(ll_y_cand(i,j,0)==ll_y_cand(i,j,0)){
              llycandsum+=ll_y_cand(i,j,0);
            }
            if(ll_y_curr(i,j,0)==ll_y_curr(i,j,0)){
              llysum+=ll_y_curr(i,j,0);
            }
          }
          rand=Rcpp::runif(1);
          if(rand(0) < exp((llycandsum+ llzcandsum)-(llysum+llzsum ))) {
            ll_y_curr.subcube(i,0,0,i,J-1,0) = ll_y_cand.subcube(i,0,0,i,J-1,0);
            Ez =clone(Ezcand);//cloning whole thing prob not most efficient. need first 2 dim
            ll_z = clone(ll_z_cand);//same here
            z(i,0) = z1cand(i);
            sumz=0;
            for(int j=0; j<t; j++){
              sumz+= z(i,j);
            }
            gammaprime(0) = gammaprimecand(0);
            if(z1cand(i)==1){
              a(i,_)=alloff;//if caught occ 1, never available to recruit
            }else if (sumz==0){  //if never caught, turn availability all on
              a(i,_)=allon;
            }else {  //if not caught on occ1, but caught later, available to recruit in occ2
              a(i,0)=1;
            }
          }
        }
      }
      N(0)+=z(i,0);
    }
    // Finally, need to update gamma.prime[,2:t], Ez[,2:t-1], and ll.z[,3:t]
    // for z1 guys now caught in first year, previously not captured (all a turned off) and guys now never caught
    // but previously only captured on occasion 1 (all a turned on)
    // But only if t>2
    if(t>2){
      for(int l=2; l<t; l++){
        suma=0;
        for(int i=0; i<M; i++){
          suma+=a(i,l-1);
        }
        gammaprime(l-1)=(N(l-1)*gammause(l-1)) / suma;
        for(int i=0; i<M; i++){
          Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
          ll_z(i,l)=z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
          if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
            ll_z(i,l)=0;
          }
        }
      }
    }
    //update z[,2+]
    for(int l=1; l<t; l++){
      //figure out who can be updated with upz
      //always remove dead guys
      if(t==2){
        for(int i=0; i<M; i++){
          upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));
        }
      }else if(t==3){
        if(l==1){
          for(int i=0; i<M; i++){
            upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));  //remove guys that are on before and after l
          }
        }else{
          for(int i=0; i<M; i++){
            upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));
          }
        }
      }else{//t>3
        if(l==(t-1)){ //can update anyone on last occasion unless they're dead
          for(int i=0; i<M; i++){
            upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));
          }
        }else if(l==(t-2)){ //second to last occasion
          for(int i=0; i<M; i++){
            upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));  //remove guys that are on before and after l
          }
        }else{//l between 2 and t-2
          for(int i=0; i<M; i++){
            latecaps(i)=0;
            for(int j=(l+2); j<t; j++){
              latecaps(i)+=z(i,j);
            }
            upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!((z(i,l)==0)&(z(i,l+1)==0)&(latecaps(i)>0)))&(!( (z(i,l-1)==0) & (a(i,l-1)==0))); //remove guys that are on before and after l or 0,0,1
          }
        }
      }
      //Can we swap?
      navail=0;
      for(int i=0; i<M; i++){
        navail+=upz(i);
      }
      int propzuse=0;
      //How many to swap?
      if(navail < propz(l-1)) {
        propzuse=navail;
      }else{
        propzuse=propz(l-1);
      }
      //get upz2
      IntegerVector upz2(navail,0);
      idx=0;
      for(int i=0; i<M; i++){
        if(upz(i)){
          upz2(idx)=i;
          idx+=1;
        }
      }
      if(navail>0){//if so, proceed with swapping
        //Who to swap?
        IntegerVector swapzidx(propzuse,0);
        IntegerVector choose=Rcpp::seq(0,(navail-1));
        swapzidx=Rcpp::RcppArmadillo::sample(choose,propzuse,FALSE);
        IntegerVector swapz(propzuse,0);
        for(int i=0; i<propzuse; i++){
          swapz(i)=upz2(swapzidx(i));
        }
        for(int i=0; i<M; i++){
          pr_zt(i)=z(i,l);
          zt_cand(i)=z(i,l);
        }
        for(int i=0; i<propzuse; i++){
          pr_zt(swapz(i))=Ez(swapz(i),l-1);
          rand=Rcpp::rbinom(1,1,pr_zt(swapz(i)));
          zt_cand(swapz(i))=rand(0);
        }
        at_cand=1*((a(_,l-1)==1)&(zt_cand==0)); //who was available on last occasion and not proposed to be captured?
        prop_probs=0;
        back_probs=0;
        ll_z_cand=clone(ll_z);
        //sum ll across j for each l and chosen i
        llycandsum=0;
        llysum=0;
        llzcandsum=0;
        llzsum=0;
        for(int i=0; i<propzuse; i++){
          prop_probs+=zt_cand(swapz(i))*log(pr_zt(swapz(i)))+(1-zt_cand(swapz(i)))*log(1-pr_zt(swapz(i)));
          back_probs+=z(swapz(i),l)*log(pr_zt(swapz(i)))+(1-z(swapz(i),l))*log(1-pr_zt(swapz(i)));
          for(int j=0; j<Xidx(l); j++){
            ll_y_cand(swapz(i),j,l)=zt_cand(swapz(i))*(y(swapz(i),j,l)*log(pd(swapz(i),j,l))+(K(0)-y(swapz(i),j,l))*log(1-pd(swapz(i),j,l)));
            //Add up y likelihood contributions curr and cand for swapped guys
            if(ll_y_cand(swapz(i),j,l)==ll_y_cand(swapz(i),j,l)){
              llycandsum+=ll_y_cand(swapz(i),j,l);
            }
            if(ll_y_curr(swapz(i),j,l)==ll_y_curr(swapz(i),j,l)){
              llysum+=ll_y_curr(swapz(i),j,l);
            }
          }
          //Add up the z likelihood contributions curr and cand for swapped guys
          ll_z_cand(swapz(i),l) = zt_cand(swapz(i))*log(Ezcand(swapz(i),l-1))+(1-zt_cand(swapz(i)))*log(1-Ezcand(swapz(i),l-1));
          if(ll_z(swapz(i),l)!=ll_z(swapz(i),l)){//Turn NaNs to 0s
            ll_z(swapz(i),l)=0;
          }
          llzcandsum+=ll_z_cand(swapz(i),l);//prior.z.cand
          llzsum+=ll_z(swapz(i),l);//prior.z
        }
        //If we make this change to z[,l] and a[,l], how does it change ll.z[,l+1]?
        if(l<(t-1)){ //extra stuff if this is true
          sumz1tmp=0;
          suma1tmp=0;
          for(int i=0; i<M; i++){
            sumz1tmp+=zt_cand(i);
            suma1tmp+=at_cand(i);
          }
          gammaprimecand(l) = sumz1tmp*gammause(l) / suma1tmp;
          if(gammaprimecand(l) <= 1){ //is gamma prime permissible?
            for(int i=0; i<M; i++){
              //Add on contributions to z ll from l+1 for cand and curr
              Ezcand(i,l) = zt_cand(i)*phiuse(l) + at_cand(i)*gammaprimecand(l);
              ll_z_cand(i,l+1) = z(i,l+1)*log(Ezcand(i,l))+(1-z(i,l+1))*log(1-Ezcand(i,l));
              if(ll_z(i,l+1)==ll_z(i,l+1)){
                llzsum+=ll_z(i,l+1);//prior.z
              }
              if(ll_z_cand(i,l+1)==ll_z_cand(i,l+1)){
                llzcandsum+=ll_z_cand(i,l+1);//prior.z.cand
              }else{
                ll_z_cand(i,l+1)=0;
              }
            }
            rand=Rcpp::runif(1);
            if(rand(0) < exp(( llycandsum+llzcandsum+back_probs)-( llysum+llzsum+prop_probs) )) {
              for(int i=0; i<propzuse; i++){
                for(int j=0; j<Xidx(l); j++){
                  ll_y_curr(swapz(i),j,l) = ll_y_cand(swapz(i),j,l);
                }
              }
              ll_z(_,l)=ll_z_cand(_,l);
              z(_,l)=zt_cand;
              a(_,l)=at_cand;
              //Now update t<l stuff
              ll_z(_,l+1)=ll_z_cand(_,l+1);
              Ez(_,l)=Ezcand(_,l);
              gammaprime(l)=gammaprimecand(l);
              if(t>3){//more houskeeping if t>3. if t=3 you can update a[,3] but it doesn't matter. can't recruit after survey is over
                // for the same reason, you don't need to update if (l+1)==t and t>3. I exclude those below so a is updated
                // for the last occasion and my sanity check works.
                // turn off availability for (l+1):t if you enter population but weren't there before
                LogicalVector fix1(swapz.size());
                bool fix3=FALSE;
                for(int i=0; i<propzuse; i++){
                  suma=0;
                  for(int l2=l; l2<t; l2++){
                    suma+=a(swapz(i),l2);
                  }
                  if((zt_cand(swapz(i))==1)&(suma>0)){
                    fix1(i)=TRUE;
                    for(int l2=l; l2<t; l2++){
                      a(swapz(i),l2)=0;
                    }
                    fix3=TRUE;
                  }else{
                    fix1(i)=FALSE;
                  }
                }
                // turn on availability for (l+1):t if you leave the population
                LogicalVector fix2(M,swapz.size());
                for(int i=0; i<propzuse; i++){
                  suma=0;
                  for(int l2=0; l2<t; l2++){
                    suma+=a(swapz(i),l2);
                  }
                  sumz=0;
                  for(int l2=0; l2<t; l2++){
                    sumz+=z(swapz(i),l2);
                  }
                  if((sumz==0)&(suma!=t)){
                    fix2(i)=TRUE;
                    for(int l2=l; l2<t; l2++){
                      a(swapz(i),l2)=1;
                    }
                    fix3=TRUE;
                  }else{
                    fix2(i)=FALSE;
                  }
                }
                // If we updated anything here,a[,(l+1):t] has changed and
                // we have to recalculate gamma.prime[,(l+1):(t-1)], Ez.cand[,(l+1):(t-1)], and ll.z.cand[,(l+2):t]
                if(((l+1)!=(t-1))){//t-1 instead of t for Rcpp
                  if(fix3){
                    for(int l2=l+1; l2<(t-1); l2++){
                      sumz1tmp=0;
                      suma1tmp=0;
                      for(int i=0; i<M; i++){
                        sumz1tmp+=z(i,l2);
                        suma1tmp+=a(i,l2);
                      }
                      gammaprime(l2) = sumz1tmp*gammause(l2) / suma1tmp;
                      Ez(_,l2)=z(_,l2)*phiuse(l2) + a(_,l2)*gammaprime(l2);
                      for(int i=0; i<M; i++){
                        //Add on contributions to z ll from l+1 for cand and curr
                        Ez(i,l2) = z(i,l2)*phiuse(l2) + a(i,l2)*gammaprime(l2);
                        ll_z(i,l2+1) = z(i,l2+1)*log(Ez(i,l2))+(1-z(i,l2+1))*log(1-Ez(i,l2));
                      }
                    }
                  }
                }
              }
            }
          }
        }else{//if l==t don't need to consider next time step
          rand=Rcpp::runif(1);
          if(rand(0) < exp(( llycandsum+llzcandsum+back_probs)-( llysum+llzsum+prop_probs) )) {
            for(int i=0; i<propzuse; i++){
              for(int j=0; j<Xidx(l); j++){
                ll_y_curr(swapz(i),j,l) = ll_y_cand(swapz(i),j,l);
              }
            }
            ll_z(_,l)=ll_z_cand(_,l);
            z(_,l)=zt_cand;
            a(_,l)=at_cand;
          }
        }
      }
    }
    // Update N
    for(int l=1; l<t; l++){
      N(l)=0;
      for(int i=0; i<M; i++){
        N(l)+=z(i,l);
      }
    }
    //Update psi
    rand=Rcpp::rbeta(1, 1+N(0), 1+M-N(1));
    psi=rand(0);
    for(int i=0; i<M; i++){
      ll_z(i,0)= z(i,0)*log(psi)+(1-z(i,0))*log(1-psi);
    }
    //Update phi
    if(phi.size()==(t-1)){//if time-specific survival
      for(int l=1; l<t; l++){
        survive=0;
        dead=0;
        for(int i=0; i<M; i++){
          survive+=(z(i,l-1)==1)&(z(i,l)==1);
          dead+=(z(i,l-1)==1)&(z(i,l)==0);
        }
        rand=Rcpp::rbeta(1, 1+survive, 1+dead);
        phi(l-1)=rand(0);
        phiuse(l-1)=rand(0);
      }
    }else{
      survive=0;
      dead=0;
      for(int l=1; l<t; l++){
        for(int i=0; i<M; i++){
          survive+=(z(i,l-1)==1)&(z(i,l)==1);
          dead+=(z(i,l-1)==1)&(z(i,l)==0);
        }
      }
      rand=Rcpp::rbeta(1, 1+survive, 1+dead);
      phi(0)=rand(0);
      for(int l=1; l<t; l++){
        phiuse(l-1)=rand(0);
      }
    }
    //  Update gamma
    // NOTE: Must update ll.z, Ez, etc...
    if(gamma.size()==1){
      rand=Rcpp::rnorm(1, gamma(0), propgamma(0));
      gamma_cand=rand(0);
      gamma_cand_ok=TRUE;
      llzsum=0;
      for(int l=1; l<t; l++){
        suma=0;
        for(int i=0; i<M; i++){
          suma+=a(i,l-1);
        }
        gamma_prime_cand(l-1)=(N(l-1)*gamma_cand) / suma;
        if(gamma_prime_cand(l-1) > 1){  //Note don't break loop b/c ll.z needs updating because phi changed
          gamma_cand_ok=FALSE;
        }
        for(int i=0; i<M; i++){
          Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
          ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
          if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
            ll_z(i,l)=0;
          }
          llzsum+=ll_z(i,l);
        }
      }
      if((gamma_cand>0)&(gamma_cand_ok)){
        llzcandsum=0;
        for(int l=1; l<t; l++){
          for(int i=0; i<M; i++){
            Ezcand(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gamma_prime_cand(l-1);
            ll_z_cand(i,l)= z(i,l)*log(Ezcand(i,l-1))+(1-z(i,l))*log(1-Ezcand(i,l-1));
            if(ll_z_cand(i,l)!=ll_z_cand(i,l)){//Turn NaNs to 0s
              ll_z_cand(i,l)=0;
            }
            llzcandsum+=ll_z_cand(i,l);
          }
        }
        rand=Rcpp::runif(1);
        if(rand(0) < exp(llzcandsum - llzsum)){
          gamma(0)=gamma_cand;
          for(int l=1; l<t; l++){
            gammaprime(l-1)=gamma_prime_cand(l-1);
            gammause(l-1)=gamma_cand;
          }
          Ez=Ezcand;
          ll_z=ll_z_cand;
        }
      }
    }else{
      for(int l=1; l<t; l++){
        rand=Rcpp::rnorm(1, gamma(l-1), propgamma(l-1));
        gamma_cand=rand(0);
        gamma_cand_ok=TRUE;
        suma=0;
        for(int i=0; i<M; i++){
          suma+=a(i,l-1);
        }
        gamma_prime_cand(l-1)=(N(l-1)*gamma_cand) / suma;
        if(gamma_prime_cand(l-1) > 1){  //Note don't break loop b/c ll.z needs updating because phi changed
          gamma_cand_ok=FALSE;
        }
        llzsum=0;
        for(int i=0; i<M; i++){
          Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
          ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
          if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
            ll_z(i,l)=0;
          }
          llzsum+=ll_z(i,l);
        }
        if((gamma_cand>0)&(gamma_cand_ok)){
          llzcandsum=0;
          for(int i=0; i<M; i++){
            Ezcand(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gamma_prime_cand(l-1);
            ll_z_cand(i,l)= z(i,l)*log(Ezcand(i,l-1))+(1-z(i,l))*log(1-Ezcand(i,l-1));
            if(ll_z_cand(i,l)!=ll_z_cand(i,l)){//Turn NaNs to 0s
              ll_z_cand(i,l)=0;
            }
            llzcandsum+=ll_z_cand(i,l);
          }
          rand=Rcpp::runif(1);
          if(rand(0) < exp(llzcandsum - llzsum)){
            gamma(l-1)=gamma_cand;
            gammaprime(l-1)=gamma_prime_cand(l-1);
            gammause(l-1)=gamma_cand;
            for(int i=0; i<M; i++){
              Ez(i,l-1)=Ezcand(i,l-1);
              ll_z(i,l-1)=ll_z_cand(i,l-1);
            }
          }
        }
      }
    }
    ////////// Now we have to update the activity centers//////////////////
    if(metamu){
      // // Update within year ACs
      for(int i=0; i<M; i++) {
        for(int l=0; l<t; l++) {
          //Don't update if not alive on this occasion
          if(z(i,l)==1){
            ScandX=Rcpp::rnorm(1,s2(i,l,0),props2x);
            ScandY=Rcpp::rnorm(1,s2(i,l,1),props2y);
            if(useverts==FALSE){
              inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
            }else{
              inbox=inoutCppOpen(ScandX,ScandY,vertices);
            }
            if(inbox(0)){
              //sum ll across j for each i and l
              llysum=0;
              llycandsum=0;
              lls2sum=0;
              lls2candsum=0;
              for(int j=0; j<Xidx(l); j++){
                dtmp(j,l)=pow( pow(ScandX(0) - Xcpp(l,j,0), 2.0) + pow(ScandY(0)-Xcpp(l,j,1), 2.0), 0.5 );
                lamdcand(i,j,l)=lam0use(l)*exp(-dtmp(j,l)*dtmp(j,l)/(2*sigmause(l)*sigmause(l)));
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
                if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                  llycandsum+=ll_y_cand(i,j,l);
                }
                if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
                  llysum+=ll_y_curr(i,j,l);
                }
              }
              ll_s2_cand(i,l)=pow(sigma_t(0),2.0)-(1/(2*pow(sigma_t(0),2.0)))*(pow(ScandX(0)-s1(i,0),2.0)+pow(ScandY(0)-s1(i,1),2.0));
              rand=Rcpp::runif(1);
              if(rand(0)<exp((llycandsum+ll_s2_cand(i,l))-(llysum+ll_s2(i,l)))){
                s2(i,l,0)=ScandX(0);
                s2(i,l,1)=ScandY(0);
                ll_s2(i,l)=ll_s2_cand(i,l);
                for(int j=0; j<J; j++){
                  D(i,j,l) = dtmp(j,l);
                  lamd(i,j,l) = lamdcand(i,j,l);
                  pd(i,j,l) = pdcand(i,j,l);
                  ll_y_curr(i,j,l) = ll_y_cand(i,j,l);
                }
              }
            }
          }
        }
      }
      // // Update meta mus
      for(int i=0; i<M; i++) {
        zany=FALSE;
        for(int l=0; l<t; l++) {
          if(z(i,l)==1){
            zany=TRUE;
          }
        }
        if(zany){
          ScandX=Rcpp::rnorm(1,s1(i,0),props1x);
          ScandY=Rcpp::rnorm(1,s1(i,1),props1y);
          if(useverts==FALSE){
            inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
          }else{
            inbox=inoutCppOpen(ScandX,ScandY,vertices);
          }
          if(inbox(0)){
            lls2sum=0;
            lls2candsum=0;
            for(int l=0; l<t; l++) {//cancel out z(i,l)=0
              ll_s2_cand(i,l)=z(i,l)*(pow(sigma_t(0),2.0)-(1/(2*pow(sigma_t(0),2.0)))*(pow(s2(i,l,0)-ScandX(0),2.0)+pow(s2(i,l,1)-ScandY(0),2.0)));
              lls2sum+=ll_s2(i,l);
              lls2candsum+=ll_s2_cand(i,l);
            }
            rand=Rcpp::runif(1);
            if (rand(0) < exp(lls2candsum - lls2sum)) {
              s1(i,0)=ScandX(0);
              s1(i,1)=ScandY(0);
              for(int l=0; l<t; l++) {
                ll_s2(i,l)=ll_s2_cand(i,l);
              }
            }
          }
        }
      }
      // Update sigma_t
      NumericVector sigma_t_cand(1);
      sigma_t_cand=rnorm(1,sigma_t(0),propsigma_t);
      if(sigma_t_cand(0) > 0){
        lls2sum=0;
        lls2candsum=0;
        for(int i=0; i<M; i++) {
          for(int l=0; l<t; l++) {
            ll_s2_cand(i,l)=pow(sigma_t_cand(0),2.0)-(1/(2*pow(sigma_t_cand(0),2.0)))*(pow(s2(i,l,0)-s1(i,0),2.0)+pow(s2(i,l,1)-s1(i,1),2.0));
            if(z(i,l)==1){//probably can not calc line above, too, if z=0
              lls2sum+=ll_s2(i,l);
              lls2candsum+=ll_s2_cand(i,l);
            }
          }
        }
        rand=Rcpp::runif(1);
        if (rand(0) < exp(lls2candsum - lls2sum)) {
          sigma_t(0)=sigma_t_cand(0);
          for(int i=0; i<M; i++) {
            for(int l=0; l<t; l++) {
              ll_s2(i,l)=ll_s2_cand(i,l);
            }
          }
        }
      }
    }else{//Stationary ACs
      //Update Activity Centers
      for(int i=0; i<M; i++) {
        //Don't update if never alive
        zany=FALSE;
        for(int l=0; l<t; l++) {
          if(z(i,l)==1){
            zany=TRUE;
          }
        }
        if(zany){
          ScandX=Rcpp::rnorm(1,s1(i,0),props2x);
          ScandY=Rcpp::rnorm(1,s1(i,1),props2y);
          if(useverts==FALSE){
            inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
          }else{
            inbox=inoutCppOpen(ScandX,ScandY,vertices);
          }
          if(inbox(0)){
            //sum ll across j for each i and l
            llysum=0;
            llycandsum=0;
            for(int l=0; l<t; l++){
              for(int j=0; j<Xidx(l); j++){
                dtmp(j,l)=pow( pow(ScandX(0) - Xcpp(l,j,0), 2.0) + pow(ScandY(0)-Xcpp(l,j,1), 2.0), 0.5 );
                lamdcand(i,j,l)=lam0use(l)*exp(-dtmp(j,l)*dtmp(j,l)/(2*sigmause(l)*sigmause(l)));
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(K[l]-y(i,j,l))*log(1-pdcand(i,j,l)));
                if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                  llycandsum+=ll_y_cand(i,j,l);
                }
                if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
                  llysum+=ll_y_curr(i,j,l);
                }
              }
            }
            rand=Rcpp::runif(1);
            if((rand(0)<exp(llycandsum-llysum))){
              s1(i,0)=ScandX(0);
              s1(i,1)=ScandY(0);
              for(int l=0; l<t; l++){
                s2(i,l,0)=ScandX(0);
                s2(i,l,1)=ScandY(0);
                for(int j=0; j<J; j++){
                  D(i,j,l) = dtmp(j,l);
                  lamd(i,j,l) = lamdcand(i,j,l);
                  pd(i,j,l) = pdcand(i,j,l);
                  ll_y_curr(i,j,l) = ll_y_cand(i,j,l);
                }
              }
            }
          }
        }
      }
    }
    //Record output ll_y_curr.subcube(i,0,0,i,maxJ-1,0) s2yout(nstore,M,t)
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      s1xout(iteridx,_)= s1(_,0);
      s1yout(iteridx,_)= s1(_,1);
      for(int i=0; i<M; i++){
        for(int l=0; l<t; l++){
          s2xout(iteridx,i,l)=s2(i,l,0);
          s2yout(iteridx,i,l)=s2(i,l,1);
          zout(iteridx,i,l)= z(i,l);
        }
      }
      idx=0;
      //fill in lam0
      for(int l=0; l<each(0); l++){
        out(iteridx,idx)=lam0(l);
        idx=idx+1;
      }
      //fill in sigma
      for(int l=0; l<each(1); l++){
        out(iteridx,idx)=sigma(l);
        idx=idx+1;
      }
      //fill in lam0
      for(int l=0; l<each(2); l++){
        out(iteridx,idx)=gamma(l);
        idx=idx+1;
      }
      //fill in lam0
      for(int l=0; l<each(3); l++){
        out(iteridx,idx)=phi(l);
        idx=idx+1;
      }
      //fill in N
      for(int l=0; l<t; l++){
        out(iteridx,idx)=N(l);
        idx=idx+1;
      }
      iteridx=iteridx+1;
    }
  }
  List to_return(11);
  to_return[0] = out;
  to_return[1] = s1xout;
  to_return[2] = s1yout;
  to_return[3] = s1xout;
  to_return[4] = s1yout;
  to_return[5] = zout;
  to_return[6] = z;
  to_return[7] = a;
  to_return[8] = Ez;
  to_return[9] = gammaprime;
  to_return[10] = N;


  return to_return;
}

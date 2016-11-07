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


//To do.  fix Ez[,1], recalculate gammaprime after z1 update

//Update psi,z
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::export]]
List upZ(IntegerMatrix z, IntegerMatrix a,IntegerMatrix knownmatrix,IntegerVector Xidx,
         arma::cube y, arma::cube pd, IntegerVector K,NumericMatrix Ez,
         NumericVector gamma,NumericVector gammaprime, NumericVector phi, double psi,IntegerVector N,int maxJ,
         NumericVector propz) {
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
  IntegerVector pr_zt(M);
  IntegerVector at_cand(M);
  NumericVector zt_cand(M);
  bool skip;
  double prop_probs;
  double back_probs;

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
      gammaprimecand(0)=sumz1tmp*gamma(0)/suma1tmp;
      if(gammaprimecand(1) > 1) { // E(Recruits) must be < nAvailable
        warn(i)=TRUE;
      }
      if(warn(i)==FALSE){
        Ezcand=clone(Ez);
        ll_z_cand=clone(ll_z);
        ll_z_cand(i,0) = z1cand(i)*log(psi)+(1-z1cand(i))*log(1-psi);
        for(int i2=0; i2<M; i2++){//update Ezcand and ll_z_cand for everyone due to gammaprimecand update
          Ezcand(i2,0)=z(i2,0)*phi(0) + a(i2,0)*gammaprimecand(0);
          ll_z_cand(i2,1) = z(i2,1)*log(Ezcand(i2,0))+(1-z(i2,1))*log(1-Ezcand(i2,0));
        }
        //Now update for focal individual
        Ezcand(i,0) = z1cand(i)*phi(0) + (1-z1cand(i))*gammaprimecand(0);
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
          ll_y_cand(i,j,1) =z1cand(i)*(y(i,j,0)*log(pd(i,j,0))+(K(0)-y(i,j,0))*log(1-pd(i,j,0)));
          if(ll_y_cand(i,j,1)==ll_y_cand(i,j,1)){
            llycandsum+=ll_y_cand(i,j,1);
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
      gammaprime(l-1)=(N(l-1)*gamma(l-1)) / suma;
      for(int i=0; i<M; i++){
        Ez(i,l-1)=z(i,l-1)*phi(l-1) + a(i,l-1)*gammaprime(l-1);
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
    if(t==2){
      for(int i=0; i<M; i++){
        upz(i)=(knownmatrix(i,0)==0);
      }
    }else if(t==3){
      if(l==1){
        for(int i=0; i<M; i++){
          upz(i)=(!((z(i,0)==1)&(z(i,t)==1)))&(knownmatrix(i,0)==0);  //remove guys that are on before and after l
        }
      }else{
        for(int i=0; i<M; i++){
          upz(i)=(knownmatrix(i,0)==0);
        }
      }
    }else{//t>3
      if(l==(t-1)){ //can update anyone on last occasion
        for(int i=0; i<M; i++){
          upz(i)=(knownmatrix(i,0)==0);
        }
      }else if(l==(t-2)){ //second to last occasion
        for(int i=0; i<M; i++){
          upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,0)==0);  //remove guys that are on before and after l
        }
      }else{//l between 2 and t-2
        for(int i=0; i<M; i++){
          latecaps(i)=0;
          for(int j=(l+2); j<t; j++){
            latecaps(i)+=z(i,j);
          }
          upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,0)==0)&(!((z(i,l)==0)&(z(i,l+1)==0)&(latecaps(i)>0))); //remove guys that are on before and after l or 0,0,1
        }
      }
    }
    //Can we swap?
    int navail=upz.size();
    if(navail>0){//if so, proceed with swapping
      int propzuse=0;
      //How many to swap?
      if(navail < propz(l-1)) {
        propzuse=navail;
      }else{
        propzuse=propz(l-1);
      }
      //Who to swap?
      IntegerVector swapzidx(propzuse,0);
      IntegerVector choose=Rcpp::seq(1,navail);
      swapzidx=Rcpp::RcppArmadillo::sample(choose,propzuse,FALSE);
      IntegerVector swapz(propzuse,0);
      for(int i=0; i<propzuse; i++){
        swapz(i)=upz(swapzidx(i)-1); //-1 since counting from 0. sample() starts at 1
      }
      pr_zt =z(_,l);//clone?
      zt_cand=clone(pr_zt);
      skip=FALSE;
      for(int i=0; i<propzuse; i++){
        pr_zt(swapz(i))=Ez(swapz(i),l-1);
        rand=Rcpp::rbinom(1,1,pr_zt(swapz(i)));
        zt_cand(swapz(i))=rand(0);
        if(zt_cand(swapz(i))==z(swapz(i),l)){
          skip=TRUE;
        }

      }
      if(!skip){
        at_cand=1*((a(_,l-1)==1)&(zt_cand==0)); //who was available on last occasion and not proposed to be captured?
        prop_probs=0;
        back_probs=0;
        ll_z_cand=clone(ll_z);
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
          ll_z_cand(swapz(i),l) = zt_cand(swapz(i))*log(Ez(swapz(i),l-1))+(1-zt_cand(swapz(i)))*log(1-Ez(swapz(i),l-1));
          if(ll_z_cand(swapz(i),l)==ll_z_cand(swapz(i),l)){
            llzcandsum+=ll_z_cand(swapz(i),l);//prior.z.cand
          }
          if(ll_z(swapz(i),l)==ll_z(swapz(i),l)){
            llzsum+=ll_z(swapz(i),l);//prior.z
          }
        }
        //If we make this change to z[,l] and a[,l], how does it change ll.z[,l+1]?
        if(l<t){ //extra stuff if this is true
          sumz1tmp=0;
          suma1tmp=0;
          for(int i=0; i<M; i++){
            sumz1tmp+=zt_cand(i);
            suma1tmp+=at_cand(i);
          }
          gammaprimecand(l) <- sumz1tmp*gamma(l) / suma1tmp;
          if(gammaprimecand(l) <= 1){ //is gamma prime permissible?
            for(int i=0; i<M; i++){
              //Add on contributions to z ll from l+1 for cand and curr
              Ezcand(i,l) <- zt_cand(i)*phi(l) + at_cand(i)*gammaprimecand(l);
              ll_z_cand(i,l+1) <- z(i,l+1)*log(Ezcand(i,l))+(1-z(i,l+1))*log(1-Ezcand(i,l));
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
              // if(t>3){//#more a houskeeping if t>3//#turn off availability if you died.
              //   if(l==2){
              //     dead=rowSums(z[swapz,l:t])==0&z[swapz,1:(l-1)]>0 //#all 0 l and later but not all zero before
              //   }else if(l==t){
              //     dead=z[swapz,t]==0&rowSums(z[swapz,1:(l-1)])>0
              //   }else{
              //     dead=rowSums(z[swapz,l:t])==0&rowSums(z[swapz,1:(l-1)])>0
              //   }
              //   a[swapz[dead],l:t]=0//#turn off availability for (l+1):t if you enter population
              //   a[swapz[z[swapz,l]==1],(l+1):t]=0
              // }
              gammaprime(l)=gammaprimecand(l);
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
  }
  List to_return(8);
  to_return[0] = ll_z;
  to_return[1] = Ez;
  to_return[2] = gammaprime;
  to_return[3] = z;
  to_return[4] = a;
  to_return[5] = N;
  to_return[6] = upz;
  to_return[7] = upz;

  return to_return;
}


// //Update activity centers
// #include <Rcpp.h>
// using namespace Rcpp;
// // [[Rcpp::export]]
// List updateACs(IntegerVector z, NumericMatrix s, IntegerVector xlim,IntegerVector ylim,NumericMatrix D, NumericMatrix lamd,
//                double lam0,double sigma, NumericMatrix y, NumericMatrix X,int K,bool useverts,NumericMatrix vertices) {
//   RNGScope scope;
//   //Preallocate
//   int M = z.size();
//   int J=X.nrow();
//   LogicalVector inbox(1);
//   NumericVector pdc(J);
//   NumericVector pdccand(J);
//   NumericVector vc(J);
//   NumericVector vcCand(J);
//   NumericVector rand(1);
//   double explldiff;
//   NumericVector dtmp(J);
//   NumericVector lamdcand(J);
//   NumericVector ScandX(1);
//   NumericVector ScandY(1);
//   double v;
//   double vcand;
//   NumericMatrix scand=Rcpp::clone(s);
//   NumericMatrix Dcand=Rcpp::clone(D);
//   NumericMatrix lamdcand=Rcpp::clone(lamd);
//
//   //  Update activity centers
//   for(int i=0; i<M; i++) {
//     ScandX=Rcpp::rnorm(1,s(i,0),.2);
//     ScandY=Rcpp::rnorm(1,s(i,1),.2);
//     if(useverts==FALSE){
//       inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
//     }else{
//       inbox=inoutCpp(ScandX,ScandY,vertices);
//     }
//     if(inbox(0)){
//       dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
//       lamdcand=lam0*exp(-dtmp*dtmp/(2*sigma*sigma));
//       if(z[i]==1){ //if in pop, otherwise keep update
//         v=0;
//         vcand=0;
//         //Calculate likelihood for original and candidate
//         for(int j=0; j<J; j++){
//           pdc(j)=1-exp(-lamd(i,j));
//           pdccand(j)=1-exp(-lamdcand(j));
//           vc(j)=y(i,j)*log(pdc(j))+(K-y(i,j))*log(1-pdc(j));
//           vcCand(j)=y(i,j)*log(pdccand(j))+(K-y(i,j))*log(1-pdccand(j));
//           if(vc(j)==vc(j)){
//             v+=vc(j);
//           }
//           if(vcCand(j)==vcCand(j)){
//             vcand+=vcCand(j);
//           }
//
//         }
//         rand=Rcpp::runif(1);
//         explldiff=exp(vcand-v);
//       }
//       if((rand(0)<explldiff)|(z[i]==0)){
//         scand(i,0)=ScandX(0);
//         scand(i,1)=ScandY(0);
//         for(int j=0; j<J; j++){
//           Dcand(i,j) = dtmp(j);
//           lamdcand(i,j) = lamdcand(j);
//         }
//       }
//     }
//   }
//   List to_return(3);
//   to_return[0] = scand;
//   to_return[1] = Dcand;
//   to_return[2] = lamdcand;
//   return to_return;
// }


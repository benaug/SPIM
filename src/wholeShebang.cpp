////////////////////////////////////in polygon functions/////////////////////
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
bool intersectCpp(NumericVector sx, NumericVector sy,NumericVector vertex1, NumericVector vertex2) {
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
bool inoutCpp(NumericVector sx,NumericVector sy,NumericMatrix vertices) {
  int count=0;
  int I=vertices.nrow();
  for(int i=0; i<(I-1); i++) {
    if(intersectCpp(sx,sy,vertices(i,_),vertices(i+1,_))){
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
///////////////////Basic SCR steps to eliminate///////
//likelihood calculation
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double SCRlik(IntegerVector z,NumericMatrix lamd,NumericMatrix y,int K) {
  //Preallocate
  int M = z.size();
  int J=lamd.ncol();
  NumericMatrix pd(M,J);
  NumericMatrix v1(M,J);
  double v2=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    if(z[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pd(i,j)=1-exp(-lamd(i,j));
        v1(i,j)=y(i,j)*log(pd(i,j))+(K-y(i,j))*log(1-pd(i,j));
        if(v1(i,j)==v1(i,j)){
          v2+=v1(i,j);
        }
      }
    }
  }
  double to_return;
  to_return = v2;
  return to_return;
}

//calculate lamd
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix calclamd(double lam, double sigma, NumericMatrix D) {
  //Preallocate
  int M = D.nrow();
  int J = D.ncol();
  NumericMatrix lamd(M,J);
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      lamd(i,j)=lam*exp(-D(i,j)*D(i,j)/(2*sigma*sigma));
    }
  }
  NumericMatrix to_return;
  to_return = lamd;
  return to_return;
}


//update lam01, lam02, and sigma

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateDfun(double lam0, double sigma, NumericMatrix lamd,
                NumericMatrix y, IntegerVector z, NumericMatrix X,int K,
                NumericMatrix D) {
  RNGScope scope;
  //Preallocate
  //int M = z.size();
  //int J=X.nrow();
  double likcurr=SCRlik( z, lamd, y, K);
  double liknew=0;
  double lam0cand=0;
  double sigmacand=0;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamdcand;
  double explldiff;
  
  //Update lam0
  rand=Rcpp::rnorm(1,lam0,0.05);
  if(rand(0) > 0){
    lam0cand=rand(0);
    lamdcand=calclamd(lam0cand,sigma,D);
    liknew=SCRlik( z, lamdcand, y, K);
    rand2=Rcpp::runif(1);
    explldiff = exp(liknew-likcurr);
    if(rand2(0)<explldiff){
      lam0=lam0cand;
      lamd=lamdcand;
      likcurr=liknew;
    }
  }
  
  //Update sigma
  rand=Rcpp::rnorm(1,sigma,0.1);
  if(rand(0) > 0){
    sigmacand=rand(0);
    lamdcand=calclamd(lam0,sigmacand,D);
    liknew=SCRlik( z, lamdcand, y, K);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(liknew-likcurr)){
      sigma=sigmacand;
      lamd=lamdcand;
      likcurr=liknew;
    }
  }
  List to_return(4);
  to_return[0] = lam0;
  to_return[1] = sigma;
  to_return[2] = lamd;
  to_return[3] = likcurr;
  return to_return;
}

//Update psi,z
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updatePsi(IntegerVector z, IntegerVector knownvector,NumericMatrix lamd,
               NumericMatrix y,IntegerMatrix X, int K,double psi) {
  RNGScope scope;
  //Preallocate
  int M = z.size();
  int J=y.ncol();
  NumericMatrix pd(M,J);
  NumericMatrix pbar(M,J);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericMatrix v1(M,J);
  //  Calculate probability of no capture and update z
  for(int i=0; i<M; i++) {
    prob0(i)=1;
    for(int j=0; j<J; j++){
      pd(i,j)=1-exp(-lamd(i,j));
      pbar(i,j)=pow(1-pd(i,j),K);
      prob0(i)*=pbar(i,j);
    }
    fc(i)=prob0(i)*psi/(prob0(i)*psi + 1-psi);
    swappable(i)=(knownvector(i)==0);
    if(swappable(i)){
      NumericVector rand=Rcpp::rbinom(1,1,fc(i));
      z(i)=rand(0);
    }
  }
  //Calculate current likelihood (no constants)
  double v=0;
  for(int i=0; i<M; i++) {
    if(z[i]==1){
      for(int j=0; j<J; j++){
        v1(i,j)=y(i,j)*log(pd(i,j))+(K-y(i,j))*log(1-pd(i,j));
        if(v1(i,j)==v1(i,j)){
          v+=v1(i,j);
        }
      }
    }
  }
  //update psi
  NumericVector psinew=Rcpp::rbeta(1, 1 + sum(z), 1 + M - sum(z));
  List to_return(3);
  to_return[0] = psinew;
  to_return[1] = v;
  to_return[2] = z;
  return to_return;
}


//Update activity centers
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateACs(IntegerVector z, NumericMatrix s, IntegerVector xlim,IntegerVector ylim,NumericMatrix D, NumericMatrix lamd,
               double lam0,double sigma, NumericMatrix y, NumericMatrix X,int K,bool useverts,NumericMatrix vertices) {
  RNGScope scope;
  //Preallocate
  int M = z.size();
  int J=X.nrow();
  LogicalVector inbox(1);
  NumericVector pdc(J);
  NumericVector pdccand(J);
  NumericVector vc(J);
  NumericVector vcCand(J);
  NumericVector rand(1);
  double explldiff;
  NumericVector dtmp(J);
  NumericVector lamdcand(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double v;
  double vcand;
  NumericMatrix snew=Rcpp::clone(s);
  NumericMatrix Dnew=Rcpp::clone(D);
  NumericMatrix lamdnew=Rcpp::clone(lamd);
  
  //  Update activity centers
  for(int i=0; i<M; i++) {
    ScandX=Rcpp::rnorm(1,s(i,0),.2);
    ScandY=Rcpp::rnorm(1,s(i,1),.2);
    if(useverts==FALSE){
      inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
    }else{
      inbox=inoutCpp(ScandX,ScandY,vertices);
    }
    if(inbox(0)){
      dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
      lamdcand=lam0*exp(-dtmp*dtmp/(2*sigma*sigma));
      if(z[i]==1){ //if in pop, otherwise keep update
        v=0;
        vcand=0;
        //Calculate likelihood for original and candidate
        for(int j=0; j<J; j++){
          pdc(j)=1-exp(-lamd(i,j));
          pdccand(j)=1-exp(-lamdcand(j));
          vc(j)=y(i,j)*log(pdc(j))+(K-y(i,j))*log(1-pdc(j));
          vcCand(j)=y(i,j)*log(pdccand(j))+(K-y(i,j))*log(1-pdccand(j));
          if(vc(j)==vc(j)){
            v+=vc(j);
          }
          if(vcCand(j)==vcCand(j)){
            vcand+=vcCand(j);
          }
          
        }
        rand=Rcpp::runif(1);
        explldiff=exp(vcand-v);
      }
      if((rand(0)<explldiff)|(z[i]==0)){
        snew(i,0)=ScandX(0);
        snew(i,1)=ScandY(0);
        for(int j=0; j<J; j++){
          Dnew(i,j) = dtmp(j);
          lamdnew(i,j) = lamdcand(j);
        }
      }
    }
  }
  List to_return(3);
  to_return[0] = snew;
  to_return[1] = Dnew;
  to_return[2] = lamdnew;
  return to_return;
}

////////////////////////////////////////Basic SCR model///////////////////////////////////////////////
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMC1(double lam0, double sigma,NumericMatrix y, NumericMatrix lamd, IntegerVector z, NumericMatrix X,int K,NumericMatrix D, IntegerVector knownvector,
           NumericMatrix s,NumericVector psi, NumericVector xlim,NumericVector ylim,bool useverts, NumericMatrix vertices, double proplam0,
           double propsigma,double propsx,double propsy, int niter, int nburn, int nthin,
           int obstype,IntegerVector tf) {
  RNGScope scope;
  int M = lamd.nrow();
  int J = lamd.ncol();
  //Preallocate for update lam0  sigma
  double lam0cand;
  double sigmacand;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamdcand(M,J);
  NumericMatrix pd(M,J);
  NumericMatrix pdcand(M,J);
  NumericMatrix ll_y_curr(M,J);
  NumericMatrix ll_y_cand(M,J);
  double llysum;
  double llycandsum;
  
  //Preallocate for Psi update
  NumericMatrix pbar(M,J);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  int N;
  
  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector dtmp(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  
  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,3);
  NumericMatrix sxout(nstore,M);
  NumericMatrix syout(nstore,M);
  NumericMatrix zout(nstore,M);
  int iteridx=0;
  //Starting log likelihood obs mod
  llysum=0;
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      if(obstype==1){
        pd(i,j)=1-exp(-lamd(i,j));
        ll_y_curr(i,j)=z(i)*(y(i,j)*log(pd(i,j))+(tf(j)-y(i,j))*log(1-pd(i,j)));
      }else{
        ll_y_curr(i,j)=z(i)*(y(i,j)*log(tf(j)*lamd(i,j))-tf(j)*lamd(i,j));
      } 
      if(ll_y_curr(i,j)==ll_y_curr(i,j)){
        llysum+=ll_y_curr(i,j); 
      }
    }
  }
  
  //start MCMC here
  int iter;
  for(iter=0; iter<niter; iter++){
    // Need to resum the ll_y on each iter
    llysum=0;
    for(int i=0; i<M; i++) {
      for(int j=0; j<J; j++){
        if(ll_y_curr(i,j)==ll_y_curr(i,j)){
          llysum+=ll_y_curr(i,j);
        }
      }
    }
    //Update lam0
    rand=Rcpp::rnorm(1,lam0,proplam0);
    if(rand(0) > 0){
      llycandsum=0;
      lam0cand=rand(0);
       // Update lamd and calculate cand likelihood
      for(int i=0; i<M; i++) {
        for(int j=0; j<J; j++){
          lamdcand(i,j)=lam0cand*exp(-D(i,j)*D(i,j)/(2*sigma*sigma));
          if(obstype==1){
            pdcand(i,j)=1-exp(-lamdcand(i,j));
            ll_y_cand(i,j)=z(i)*(y(i,j)*log(pdcand(i,j))+(tf(j)-y(i,j))*log(1-pdcand(i,j)));
          }else{
            ll_y_cand(i,j)=z(i)*(y(i,j)*log(tf(j)*lamdcand(i,j))-tf(j)*lamdcand(i,j));
          }
          if(ll_y_cand(i,j)==ll_y_cand(i,j)){
            llycandsum+=ll_y_cand(i,j);
          }
        }
      }
      //MH step
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(llycandsum-llysum)){
        lam0=lam0cand;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd(i,j)=lamdcand(i,j);
            if(obstype==1){
              pd(i,j)=pdcand(i,j);
            }
            ll_y_curr(i,j)=ll_y_cand(i,j);
          }
        }
        llysum=llycandsum;
      }
    }

    // Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      llycandsum=0;
      //  Update lamd and calculate cand likelihood
      for(int i=0; i<M; i++) {
        for(int j=0; j<J; j++){
          lamdcand(i,j)=lam0*exp(-D(i,j)*D(i,j)/(2*sigmacand*sigmacand));
          if(obstype==1){
            pdcand(i,j)=1-exp(-lamdcand(i,j));
            ll_y_cand(i,j)=z(i)*(y(i,j)*log(pdcand(i,j))+(tf(j)-y(i,j))*log(1-pdcand(i,j)));
          }else{
            ll_y_cand(i,j)=z(i)*(y(i,j)*log(tf(j)*lamdcand(i,j))-tf(j)*lamdcand(i,j));
          }
          if(ll_y_cand(i,j)==ll_y_cand(i,j)){
            llycandsum+=ll_y_cand(i,j);
          }
        }
      }
      rand=Rcpp::runif(1);
      if(rand(0)<exp(llycandsum-llysum)){
        sigma=sigmacand;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd(i,j)=lamdcand(i,j);
            if(obstype==1){
              pd(i,j)=pdcand(i,j);
            }
            ll_y_curr(i,j)=ll_y_cand(i,j);
          }
        }
        llysum=llycandsum;
      }
    }
    //Update Psi
    //  Calculate probability of no capture and update z
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        if(obstype==2){
          pd(i,j)=1-exp(-lamd(i,j));
        }
        pbar(i,j)=pow(1-pd(i,j),tf(j));
        prob0(i)*=pbar(i,j);
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0);
      if(swappable(i)){
        NumericVector rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
    }
    //Calculate updated obs mod ll
    llysum=0;
    for(int i=0; i<M; i++) {
      for(int j=0; j<J; j++){
        if(obstype==1){
          pd(i,j)=1-exp(-lamd(i,j));
          ll_y_curr(i,j)=z(i)*(y(i,j)*log(pd(i,j))+(tf(j)-y(i,j))*log(1-pd(i,j)));
        }else{
          ll_y_curr(i,j)=z(i)*(y(i,j)*log(tf(j)*lamd(i,j))-tf(j)*lamd(i,j));
        }
        if(ll_y_curr(i,j)==ll_y_curr(i,j)){
          llysum+=ll_y_curr(i,j);
        }
      }
    }
    //update psi
    N=sum(z);
    psi=Rcpp::rbeta(1, 1 + N, 1 + M - N);

    //Update Activity Centers
    for(int i=0; i<M; i++) {
      ScandX=Rcpp::rnorm(1,s(i,0),propsx);
      ScandY=Rcpp::rnorm(1,s(i,1),propsy);
      if(useverts==FALSE){
        inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      }else{
        inbox=inoutCpp(ScandX,ScandY,vertices);
      }
      if(inbox(0)){
        llysum=0;
        llycandsum=0;
        for(int j=0; j<J; j++){
          dtmp(j)=pow(pow(ScandX(0)-X(j,0), 2.0)+pow(ScandY(0)-X(j,1),2.0),0.5);
          lamdcand(i,j)=lam0*exp(-dtmp(j)*dtmp(j)/(2*sigma*sigma));
          if(obstype==1){
            pdcand(i,j)=1-exp(-lamdcand(i,j));
            ll_y_cand(i,j)=z(i)*(y(i,j)*log(pdcand(i,j))+(tf(j)-y(i,j))*log(1-pdcand(i,j)));
          }else{
            ll_y_cand(i,j)=z(i)*(y(i,j)*log(tf(j)*lamdcand(i,j))-tf(j)*lamdcand(i,j));
          }
          if(ll_y_cand(i,j)==ll_y_cand(i,j)){
            llycandsum+=ll_y_cand(i,j);
          }
          if(ll_y_curr(i,j)==ll_y_curr(i,j)){
            llysum+=ll_y_curr(i,j);
          }
        }
        rand=Rcpp::runif(1);
        if((rand(0)<exp(llycandsum-llysum))){
          s(i,0)=ScandX(0);
          s(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D(i,j) = dtmp(j);
            lamd(i,j) = lamdcand(i,j);
            if(obstype==1){
              pd(i,j) = pdcand(i,j);
            }
            ll_y_curr(i,j) = ll_y_cand(i,j);
          }
        }
      }
    }
  //Record output
  if(((iter+1)>nburn)&((iter+1) % nthin==0)){
    sxout(iteridx,_)= s(_,0);
    syout(iteridx,_)= s(_,1);
    zout(iteridx,_)= z;
    out(iteridx,0)=lam0;
    out(iteridx,1)=sigma;
    out(iteridx,2)=N;
    iteridx=iteridx+1;
  }
}
List to_return(4);
to_return[0] = out;
to_return[1] = sxout;
to_return[2] = syout;
to_return[3] = zout;
return to_return;
}

///////////////////////////Add trap history functionality////////////////////////
//likelihood calculation
#include <Rcpp.h>
    using namespace Rcpp;
//[[Rcpp::export]]
double SCRliktf(IntegerVector z,NumericMatrix lamd,NumericMatrix y,NumericVector K) {
  //Preallocate
  int M = z.size();
  int J=lamd.ncol();
  NumericMatrix pd(M,J);
  NumericMatrix v1(M,J);
  double v2=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    if(z[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pd(i,j)=1-exp(-lamd(i,j));
        v1(i,j)=y(i,j)*log(pd(i,j))+(K(j)-y(i,j))*log(1-pd(i,j));
        if(v1(i,j)==v1(i,j)){
          v2+=v1(i,j);
        }
      }
    }
  }
  double to_return;
  to_return = v2;
  return to_return;
}


////////////////////////////////////////Side matching model///////////////////////////////////////////
//Individual components first. Full, combined algorithm below.

//likelihood calculation
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fulllik(IntegerVector z,NumericMatrix lamd1,NumericMatrix lamd2,NumericMatrix yboth, NumericMatrix yleft,
               NumericMatrix yright, NumericMatrix X,int K) {
  //Preallocate
  int M = z.size();
  int J=X.nrow();
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix partboth(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);
  double v=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    if(z[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(i,j));
        if(X(j,2)==2){ //if double trap
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          pd2(i,j)=1-exp(-lamd2(i,j));
          partboth(i,j)=yboth(i,j)*log(pd2(i,j))+(K-yboth(i,j))*log(1-pd2(i,j));
          partleft(i,j)=yleft(i,j)*log(pd12(i,j))+(K-yleft(i,j))*log(1-pd12(i,j));
          partright(i,j)=yright(i,j)*log(pd12(i,j))+(K-yright(i,j))*log(1-pd12(i,j));
        }else{ //if single trap
          partboth(i,j)=0;
          partleft(i,j)=yleft(i,j)*log(pd1(i,j))+(K-yleft(i,j))*log(1-pd1(i,j));
          partright(i,j)=yright(i,j)*log(pd1(i,j))+(K-yright(i,j))*log(1-pd1(i,j));
        }
        if(partboth(i,j)==partboth(i,j)){
          v+=partboth(i,j);
        }
        if(partleft(i,j)==partleft(i,j)){
          v+=partleft(i,j);
        }
        if(partright(i,j)==partright(i,j)){
          v+=partright(i,j);
        }
      }
    }
  }
  double to_return;
  to_return = v;
  return to_return;
}

//update lam01, lam02, and sigma

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateparms(double lam01,double lam02, double sigma, NumericMatrix lamd1,NumericMatrix lamd2,
                 NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright, IntegerVector z, NumericMatrix X,int K,
                 NumericMatrix D) {
  RNGScope scope;
  //Preallocate
  //int M = z.size();
  //int J=X.nrow();
  double likcurr=fulllik( z, lamd1, lamd2, yboth,  yleft, yright,  X, K);
  double liknew=0;
  double lam01cand=0;
  double lam02cand=0;
  double sigmacand=0;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd1cand;
  NumericMatrix lamd2cand;
  double explldiff;

  //Update lam01
  rand=Rcpp::rnorm(1,lam01,0.05);
  if(rand(0) > 0){
    lam01cand=rand(0);
    lamd1cand=calclamd(lam01cand,sigma,D);
    liknew=fulllik( z, lamd1cand, lamd2, yboth,  yleft, yright,  X, K);
    rand2=Rcpp::runif(1);
    explldiff = exp(liknew-likcurr);
    if(rand2(0)<explldiff){
      lam01=lam01cand;
      lamd1=lamd1cand;
      likcurr=liknew;
    }
  }
  //Update lam02
  rand=Rcpp::rnorm(1,lam02,0.05);
  if(rand(0) > 0){
    lam02cand=rand(0);
    lamd2cand=calclamd(lam02cand,sigma,D);
    liknew=fulllik( z, lamd1, lamd2cand, yboth,  yleft, yright,  X, K);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(liknew-likcurr)){
      lam02=lam02cand;
      lamd2=lamd2cand;
      likcurr=liknew;
    }
  }
  //Update sigma
  rand=Rcpp::rnorm(1,sigma,0.1);
  if(rand(0) > 0){
    sigmacand=rand(0);
    lamd1cand=calclamd(lam01,sigmacand,D);
    lamd2cand=calclamd(lam02,sigmacand,D);
    liknew=fulllik( z, lamd1cand, lamd2cand, yboth,  yleft, yright,  X, K);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(liknew-likcurr)){
      sigma=sigmacand;
      lamd1=lamd1cand;
      lamd2=lamd2cand;
      likcurr=liknew;
    }
  }
  List to_return(6);
  to_return[0] = lam01;
  to_return[1] = lam02;
  to_return[2] = sigma;
  to_return[3] = lamd1;
  to_return[4] = lamd2;
  to_return[5] = likcurr;
  return to_return;
}

//update IDs
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List swapC(IntegerVector z,int Nfixed,IntegerVector ID_L, IntegerVector ID_R,NumericMatrix s,int swap,double swaptol,
           NumericMatrix lamd1,NumericMatrix lamd2, NumericMatrix yboth,NumericMatrix yleft, NumericMatrix yright,IntegerMatrix X, int K,
           IntegerVector left, IntegerVector right,double likcurr) {
  int M = z.size();
  RNGScope scope;
  //  Build map and candidate map for L and R separately

  //Build mapL
  IntegerMatrix mapL(M,2);
  IntegerVector v=seq_len(M);
  mapL(_,1)=ID_L;
  mapL(_,0)=v;
  IntegerMatrix candmapL=Rcpp::clone(mapL);
  for (int i=0; i<M; i++) {
    if(z[candmapL(i,1)-1]==0){
      candmapL(i,1)= -1;
    }
    if(i<Nfixed){
      candmapL(i,0)= -1;
      candmapL(i,1)= -1;
    }
    if(z(i)==0){
      candmapL(i,0)= -1;
    }
  }

  //Find outcandsL and incandsL
  IntegerVector outcandsL=v[candmapL(_,1)>0];
  IntegerVector incandsL=v[candmapL(_,0)>0];
  int insizeL=incandsL.size();
  int outsizeL=outcandsL.size();

  //Build mapR
  IntegerMatrix mapR(M,2);
  mapR(_,1)=ID_R;
  mapR(_,0)=v;
  IntegerMatrix candmapR=Rcpp::clone(mapR);
  for (int i=0; i<M; i++) {
    if(z[candmapR(i,1)-1]==0){
      candmapR(i,1)= -1;
    }
    if(i<Nfixed){
      candmapR(i,0)= -1;
      candmapR(i,1)= -1;
    }
    if(z(i)==0){
      candmapR(i,0)= -1;
    }
  }

  //Find outcandsR and incandsR
  IntegerVector outcandsR=v[candmapR(_,1)>0];
  IntegerVector incandsR=v[candmapR(_,0)>0];
  int insizeR=incandsR.size();
  int outsizeR=outcandsR.size();

  //Build left swapout bins for sample function
  //swapin bins vary in size so must be calculated in for loop
  NumericVector binsL(outsizeL+1);
  double x;
  for (int i=0; i<(outsizeL+1); i++) {
    binsL(i)=i;
    x=binsL(i)/outsizeL;
    binsL(i)=x;
  }

  //Build right swapout bins for sample function
  //swapin bins vary in size so must be calculated in for loop
  NumericVector binsR(outsizeR+1);
  for (int i=0; i<(outsizeR+1); i++) {
    binsR(i)=i;
    x=binsR(i)/outsizeR;
    binsR(i)=x;
  }

  //Prespecify storage for objects that dont change size in loop
  //Swapping structures
  NumericVector rand(1);
  int guy1;
  int guy2;
  int swapin;
  int swapout;
  int idx=0;
  double ncand;
  double ncand2;
  int match;
  double jumpprob;
  double backprob;
  NumericVector Dl(incandsL.size());
  NumericVector Dr(incandsR.size());
  IntegerVector possible;
  IntegerVector unpossible;
  NumericVector trashL(insizeL);
  NumericVector trashR(insizeR);
  IntegerVector newID(M);
  int J=yleft.ncol();
  NumericMatrix ylefttmp(M,J);
  NumericMatrix yrighttmp(M,J);  //could use same structure but dont want to get confused

  //ll structures
  NumericMatrix prepart(2,J);
  NumericMatrix postpart(2,J);
  NumericMatrix pd1(2,J);
  NumericMatrix pd12(2,J);
  NumericVector swapped(2);
  double llpre;
  double llpost;

  //MH structures
  double lldiff;

  //Housekeeping structures
  NumericVector guycounts(M);
  LogicalVector zeroguys(M);

  ///////////////// //swap left sides//////////////
  for (int l=0; l<swap; l++) {
    //Find guy 1
    rand=Rcpp::runif(1);
    idx=0;
    for (int i=0; i<(outsizeL+1); i++) {
      if(binsL(i)<rand(0)){
        idx=idx+1;
      }
    }
    guy1=outcandsL(idx-1);

    //find swapout
    swapout=mapL(guy1-1,1);

    //calculate distances and find possible switches
    for (int i=0; i<insizeL; i++) {
      Dl(i) = sqrt( pow( s(swapout-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsL(i)-1,1), 2.0) );
    }
    possible=incandsL[Dl < swaptol];
    ncand=possible.size();
    jumpprob=1/ncand;

    //swap in
    if(ncand>1){
      NumericVector binsL2(ncand+1);
      rand=Rcpp::runif(1);
      idx=0;
      for (int i=0; i<(ncand+1); i++) {
        binsL2(i)=i;
        x=binsL2(i)/ncand;
        binsL2(i)=x;
        if(x<rand(0)){
          idx=idx+1;
        }
      }
      swapin=possible(idx-1);
    }else{
      swapin=possible(0);
    }
    for (int i=0; i<incandsL.size(); i++) {
      trashL(i) = sqrt( pow( s(swapin-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsL(i)-1,1), 2.0) );
    }
    unpossible=incandsL[trashL < swaptol];
    ncand2=unpossible.size();
    backprob=1/ncand2;
    //find guy2
    for (int i=0; i<M; i++) {
      match=mapL(i,1);
      if(match==swapin){
        guy2=i+1;
      }
    }
    //update ID_L
    newID=Rcpp::clone(ID_L);
    newID(guy1-1)=swapin;
    newID(guy2-1)=swapout;

    //Make new left data set
    idx=0;
    for (int i=0; i<M; i++) {
      for(int j=0; j<J; j++) {
        ylefttmp(newID(i)-1,j)=0;
        for (int k=0; k<K; k++) {
          ylefttmp(newID(i)-1,j)+=left(idx);
          idx+=1;
        }
      }
    }
    //Calculate likelihoods pre and post switch. Only need to calc for left data
    swapped(0)=swapout;
    swapped(1)=swapin;
    llpre=0;
    llpost=0;
    for(int i=0; i<2; i++) {
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
        if(X(j,2)==2){ //if double trap
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          prepart(i,j)=yleft(swapped(i)-1,j)*log(pd12(i,j))+(K-yleft(swapped(i)-1,j))*log(1-pd12(i,j));
          postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd12(i,j))+(K-ylefttmp(swapped(i)-1,j))*log(1-pd12(i,j));
        }else{ //if single trap
          prepart(i,j)=yleft(swapped(i)-1,j)*log(pd1(i,j))+(K-yleft(swapped(i)-1,j))*log(1-pd1(i,j));
          postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd1(i,j))+(K-ylefttmp(swapped(i)-1,j))*log(1-pd1(i,j));
        }
        if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
          llpre+=prepart(i,j);
        }
        if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
          llpost+=postpart(i,j);
        }
      }
    }

    //MH step
    rand=Rcpp::runif(1);
    lldiff=llpost-llpre;
    if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
      likcurr+=lldiff;
      yleft=Rcpp::clone(ylefttmp);
      ID_L=Rcpp::clone(newID);
      mapL(guy1-1,1)=swapin;
      mapL(guy2-1,1)=swapout;
    }

  }
  ///////////////// //swap right sides/////////////////
  for (int l=0; l<swap; l++) {
    //Find guy 1
    rand=Rcpp::runif(1);
    idx=0;
    for (int i=0; i<(outsizeR+1); i++) {
      if(binsR(i)<rand(0)){
        idx=idx+1;
      }
    }
    guy1=outcandsR(idx-1);

    //find swapout
    swapout=mapR(guy1-1,1);

    //calculate distances and find possible switches
    for (int i=0; i<insizeR; i++) {
      Dr(i) = sqrt( pow( s(swapout-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsR(i)-1,1), 2.0) );
    }
    possible=incandsR[Dr < swaptol];
    ncand=possible.size();
    jumpprob=1/ncand;

    //swap in
    if(ncand>1){
      NumericVector binsR2(ncand+1);
      rand=Rcpp::runif(1);
      idx=0;
      for (int i=0; i<(ncand+1); i++) {
        binsR2(i)=i;
        x=binsR2(i)/ncand;
        binsR2(i)=x;
        if(x<rand(0)){
          idx=idx+1;
        }
      }
      swapin=possible(idx-1);
    }else{
      swapin=possible(0);
    }
    for (int i=0; i<incandsR.size(); i++) {
      trashR(i) = sqrt( pow( s(swapin-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsR(i)-1,1), 2.0) );
    }
    unpossible=incandsR[trashR < swaptol];
    ncand2=unpossible.size();
    backprob=1/ncand2;
    //find guy2
    for (int i=0; i<M; i++) {
      match=mapR(i,1);
      if(match==swapin){
        guy2=i+1;
      }
    }
    //update ID_R
    newID=Rcpp::clone(ID_R);
    newID(guy1-1)=swapin;
    newID(guy2-1)=swapout;

    //Make new right data set
    idx=0;
    for (int i=0; i<M; i++) {
      for(int j=0; j<J; j++) {
        yrighttmp(newID(i)-1,j)=0;
        for (int k=0; k<K; k++) {
          yrighttmp(newID(i)-1,j)+=right(idx);
          idx+=1;
        }
      }
    }
    //Calculate likelihoods pre and post switch. Only need to calc for right side data
    swapped(0)=swapout;
    swapped(1)=swapin;
    llpre=0;
    llpost=0;
    for(int i=0; i<2; i++){
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
        if(X(j,2)==2){ //if double trap
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          prepart(i,j)=yright(swapped(i)-1,j)*log(pd12(i,j))+(K-yright(swapped(i)-1,j))*log(1-pd12(i,j));
          postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd12(i,j))+(K-yrighttmp(swapped(i)-1,j))*log(1-pd12(i,j));
        }else{ //if single trap
          prepart(i,j)=yright(swapped(i)-1,j)*log(pd1(i,j))+(K-yright(swapped(i)-1,j))*log(1-pd1(i,j));
          postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd1(i,j))+(K-yrighttmp(swapped(i)-1,j))*log(1-pd1(i,j));
        }
        if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
          llpre+=prepart(i,j);
        }
        if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
          llpost+=postpart(i,j);
        }
      }
    }

    //MH step
    rand=Rcpp::runif(1);
    lldiff=llpost-llpre;
    if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
      likcurr+=lldiff;
      yright=Rcpp::clone(yrighttmp);
      ID_R=Rcpp::clone(newID);
      mapR(guy1-1,1)=swapin;
      mapR(guy2-1,1)=swapout;
    }

  } //end swapping rights
  //}
  //Recalculate zeroguys
  for (int i=0; i<M; i++) {
    guycounts(i)=0;
    for(int j=0; j<J; j++) {
      guycounts(i)+=yboth(i,j)+yright(i,j)+yleft(i,j);
    }
    zeroguys(i)=guycounts(i)==0;
  }

  List to_return(6);
  to_return[0] = likcurr;
  to_return[1] = ID_L;
  to_return[2] = ID_R;
  to_return[3] = yleft;
  to_return[4] = yright;
  to_return[5] = zeroguys;
  return to_return;
}

//Update psi,z
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updatePsi2(IntegerVector z, LogicalVector zeroguys, IntegerVector knownvector,NumericMatrix lamd1,
               NumericMatrix lamd2, NumericMatrix yboth,NumericMatrix yleft,
               NumericMatrix yright,IntegerMatrix X, int K,double psi) {
  RNGScope scope;
  //Preallocate
  int M = z.size();
  int J=yleft.ncol();
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix pd1traps(M,J);
  NumericMatrix pd2traps(M,J);
  NumericMatrix pbar1(M,J);
  NumericMatrix pbar2(M,J);
  NumericMatrix pbar0(M,J);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericMatrix partboth(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);
  //  Calculate probability of no capture and update z
  for(int i=0; i<M; i++) {
    prob0(i)=1;
    for(int j=0; j<J; j++){
      pd1(i,j)=1-exp(-lamd1(i,j));
      if(X(j,2)==2){ //if double trap
        pd2(i,j)=1-exp(-lamd2(i,j));
        pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
        pd1traps(i,j)=pd12(i,j);
        pd2traps(i,j)=pd2(i,j);
      }else{ //if single trap
        pd1traps(i,j)=pd1(i,j);
        pd2traps(i,j)=0;
      }
      pbar1(i,j)=pow(1-pd1traps(i,j),2*K);
      pbar2(i,j)=pow(1-pd2traps(i,j),K);
      pbar0(i,j)=pbar1(i,j)*pbar2(i,j);
      prob0(i)*=pbar0(i,j);
    }
    fc(i)=prob0(i)*psi/(prob0(i)*psi + 1-psi);
    swappable(i)=(knownvector(i)==0)&zeroguys(i);
    if(swappable(i)){
      NumericVector rand=Rcpp::rbinom(1,1,fc(i));
      z(i)=rand(0);
    }
  }
  //Calculate current likelihood (no constants)
  double v=0;
  for(int i=0; i<M; i++) {
    if(z[i]==1){
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          partboth(i,j)=yboth(i,j)*log(pd2(i,j))+(K-yboth(i,j))*log(1-pd2(i,j));
          partleft(i,j)=yleft(i,j)*log(pd12(i,j))+(K-yleft(i,j))*log(1-pd12(i,j));
          partright(i,j)=yright(i,j)*log(pd12(i,j))+(K-yright(i,j))*log(1-pd12(i,j));
        }else{ //if single trap
          partboth(i,j)=0;
          partleft(i,j)=yleft(i,j)*log(pd1(i,j))+(K-yleft(i,j))*log(1-pd1(i,j));
          partright(i,j)=yright(i,j)*log(pd1(i,j))+(K-yright(i,j))*log(1-pd1(i,j));
        }
        if(partboth(i,j)==partboth(i,j)){
          v+=partboth(i,j);
        }
        if(partleft(i,j)==partleft(i,j)){
          v+=partleft(i,j);
        }
        if(partright(i,j)==partright(i,j)){
          v+=partright(i,j);
        }
      }
    }
  }
  //update psi
  NumericVector psinew=Rcpp::rbeta(1, 1 + sum(z), 1 + M - sum(z));
  List to_return(3);
  to_return[0] = psinew;
  to_return[1] = v;
  to_return[2] = z;
  return to_return;
}

//Update activity centers
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateACs2(IntegerVector z, NumericMatrix s, IntegerVector xlim,IntegerVector ylim,bool useverts, NumericMatrix vertices,
               NumericMatrix D, NumericMatrix lamd1,NumericMatrix lamd2, double lam01, double lam02,
               double sigma, NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright, NumericMatrix X,int K) {
  RNGScope scope;
  //Preallocate
  int M = z.size();
  int J=X.nrow();
  LogicalVector inbox(1);
  NumericVector pd1c(J);
  NumericVector pd12c(J);
  NumericVector pd2c(J);
  NumericVector pd1ccand(J);
  NumericVector pd12ccand(J);
  NumericVector pd2ccand(J);
  NumericVector partbothc(J);
  NumericVector partleftc(J);
  NumericVector partrightc(J);
  NumericVector partbothcCand(J);
  NumericVector partleftcCand(J);
  NumericVector partrightcCand(J);
  NumericVector rand(1);
  double explldiff;
  NumericVector dtmp(J);
  NumericVector lamd1cand(J);
  NumericVector lamd2cand(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double v;
  double vcand;
  NumericMatrix snew=Rcpp::clone(s);
  NumericMatrix Dnew=Rcpp::clone(D);
  NumericMatrix lamd1new=Rcpp::clone(lamd1);
  NumericMatrix lamd2new=Rcpp::clone(lamd2);


  //  Update activity centers
  for(int i=0; i<M; i++) {
    ScandX=Rcpp::rnorm(1,s(i,0),.2);
    ScandY=Rcpp::rnorm(1,s(i,1),.2);
    if(useverts==FALSE){
      inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
    }else{
      inbox=inoutCpp(ScandX,ScandY,vertices);
    }
    if(inbox(0)){
        dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
        lamd1cand=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
        lamd2cand=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
        if(z[i]==1){ //if in pop, otherwise keep update
          v=0;
          vcand=0;
          //Calculate likelihood for original and candidate
          for(int j=0; j<J; j++){
            pd1c(j)=1-exp(-lamd1(i,j));
            pd1ccand(j)=1-exp(-lamd1cand(j));
            if(X(j,2)==2){ //if double trap
              pd12c(j)=2*pd1c(j)-pd1c(j)*pd1c(j);
              pd12ccand(j)=2*pd1ccand(j)-pd1ccand(j)*pd1ccand(j);
              pd2c(j)=1-exp(-lamd2(i,j));
              pd2ccand(j)=1-exp(-lamd2cand(j));
              partbothc(j)=yboth(i,j)*log(pd2c(j))+(K-yboth(i,j))*log(1-pd2c(j));
              partleftc(j)=yleft(i,j)*log(pd12c(j))+(K-yleft(i,j))*log(1-pd12c(j));
              partrightc(j)=yright(i,j)*log(pd12c(j))+(K-yright(i,j))*log(1-pd12c(j));
              partbothcCand(j)=yboth(i,j)*log(pd2ccand(j))+(K-yboth(i,j))*log(1-pd2ccand(j));
              partleftcCand(j)=yleft(i,j)*log(pd12ccand(j))+(K-yleft(i,j))*log(1-pd12ccand(j));
              partrightcCand(j)=yright(i,j)*log(pd12ccand(j))+(K-yright(i,j))*log(1-pd12ccand(j));
            }else{ //if single trap
              partbothc(j)=0;
              partleftc(j)=yleft(i,j)*log(pd1c(j))+(K-yleft(i,j))*log(1-pd1c(j));
              partrightc(j)=yright(i,j)*log(pd1c(j))+(K-yright(i,j))*log(1-pd1c(j));
              partleftcCand(j)=yleft(i,j)*log(pd1ccand(j))+(K-yleft(i,j))*log(1-pd1ccand(j));
              partrightcCand(j)=yright(i,j)*log(pd1ccand(j))+(K-yright(i,j))*log(1-pd1ccand(j));
            }
            if(partbothc(j)==partbothc(j)){
              v+=partbothc(j);
            }
            if(partleftc(j)==partleftc(j)){
              v+=partleftc(j);
            }
            if(partrightc(j)==partrightc(j)){
              v+=partrightc(j);
            }
            if(partbothcCand(j)==partbothcCand(j)){
              vcand+=partbothcCand(j);
            }
            if(partleftcCand(j)==partleftcCand(j)){
              vcand+=partleftcCand(j);
            }
            if(partrightcCand(j)==partrightcCand(j)){
              vcand+=partrightcCand(j);
            }
          }
          rand=Rcpp::runif(1);
          explldiff=exp(vcand-v);
        }
        if((rand(0)<explldiff)|(z[i]==0)){
          snew(i,0)=ScandX(0);
          snew(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            Dnew(i,j) = dtmp(j);
            lamd1new(i,j) = lamd1cand(j);
            lamd2new(i,j) = lamd2cand(j);
          }
        }
      }
    }
    List to_return(4);
    to_return[0] = snew;
    to_return[1] = Dnew;
    to_return[2] = lamd1new;
    to_return[3] = lamd2new;
    return to_return;
}



////////////////////////WholeShebang////////////////////////////////
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMC(double lam01,double lam02, double sigma,
          NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright, IntegerVector z, NumericMatrix X,int K,
          NumericMatrix D,int Nfixed, IntegerVector knownvector,IntegerVector ID_L, IntegerVector ID_R,int swap,double swaptol,
          IntegerVector left,IntegerVector right,NumericMatrix s,NumericVector psi, NumericVector xlim,NumericVector ylim,
          bool useverts, NumericMatrix vertices, double proplam01, double proplam02, double propsigma, double propsx, double propsy,
          int niter, int nburn, int nthin, LogicalVector updates) {
  RNGScope scope;
  int M = z.size();
  //Preallocate for update lam01 lam02 sigma
  double likcurr;//=fulllik( z, lamd1, lamd2, yboth,  yleft, yright,  X, K);
  double liknew;
  double lam01cand;
  double lam02cand;
  double sigmacand;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd1cand;
  NumericMatrix lamd2cand;
  //Preallocate side swapping structures
  //Swapping structures
  IntegerVector IDs=seq_len(M);
  int guy1=0;
  int guy2=0;
  int swapin;
  int swapout;
  int idx;
  double ncand;
  double ncand2;
  int match;
  double jumpprob;
  double backprob;
  IntegerVector possible;
  IntegerVector unpossible;
  IntegerVector newID(M);
  int J=yleft.ncol();
  NumericMatrix ylefttmp(M,J);
  NumericMatrix yrighttmp(M,J);  //could use same structure but dont want to get confused

  //ll structures
  NumericMatrix prepart(2,J);
  NumericMatrix postpart(2,J);
  NumericMatrix pd1b(2,J);
  NumericMatrix pd12b(2,J);
  NumericVector swapped(2);
  double llpre;
  double llpost;

  //MH structures
  double lldiff=0;

  //Housekeeping structures
  NumericVector guycounts(M);
  LogicalVector zeroguys(M);

  //Preallocate for Psi update
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix pd1traps(M,J);
  NumericMatrix pd2traps(M,J);
  NumericMatrix pbar1(M,J);
  NumericMatrix pbar2(M,J);
  NumericMatrix pbar0(M,J);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericMatrix partboth(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);
  int N;

  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector pd1c(J);
  NumericVector pd12c(J);
  NumericVector pd2c(J);
  NumericVector pd1ccand(J);
  NumericVector pd12ccand(J);
  NumericVector pd2ccand(J);
  NumericVector partbothc(J);
  NumericVector partleftc(J);
  NumericVector partrightc(J);
  NumericVector partbothcCand(J);
  NumericVector partleftcCand(J);
  NumericVector partrightcCand(J);
  NumericVector dtmp(J);
  NumericVector lamd1candc(J);
  NumericVector lamd2candc(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double vcand;
  double v;
  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,4);
  NumericMatrix sxout(nstore,M);
  NumericMatrix syout(nstore,M);
  NumericMatrix zout(nstore,M);
  NumericMatrix ID_Lout(nstore,M);
  NumericMatrix ID_Rout(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamd1=calclamd(lam01,sigma,D);
  NumericMatrix lamd2=calclamd(lam02,sigma,D);
  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    //update lam01, lam02, and sigma

    likcurr=fulllik( z, lamd1, lamd2, yboth,  yleft, yright,  X, K);
    liknew=0;
    lam01cand=0;
    lam02cand=0;
    sigmacand=0;

    //Update lam01
    if(updates(0)){
      rand=Rcpp::rnorm(1,lam01,proplam01);
      if(rand(0) > 0){
        lam01cand=rand(0);
        lamd1cand=calclamd(lam01cand,sigma,D);
        liknew=fulllik( z, lamd1cand, lamd2, yboth,  yleft, yright,  X, K);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam01=lam01cand;
          lamd1=lamd1cand;
          likcurr=liknew;
        }
      }
    }
    //Update lam02
    if(updates(1)){
      rand=Rcpp::rnorm(1,lam02,proplam02);
      if(rand(0) > 0){
        lam02cand=rand(0);
        lamd2cand=calclamd(lam02cand,sigma,D);
        liknew=fulllik( z, lamd1, lamd2cand, yboth,  yleft, yright,  X, K);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam02=lam02cand;
          lamd2=lamd2cand;
          likcurr=liknew;
        }
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamd1cand=calclamd(lam01,sigmacand,D);
      lamd2cand=calclamd(lam02,sigmacand,D);
      liknew=fulllik( z, lamd1cand, lamd2cand, yboth,  yleft, yright,  X, K);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamd1=lamd1cand;
        lamd2=lamd2cand;
        likcurr=liknew;
      }
    }

    //  Update left sides
    if(updates(2)){
      //Build mapL
      IntegerMatrix mapL(M,2);
      mapL(_,1)=ID_L;
      mapL(_,0)=IDs;
      IntegerMatrix candmapL=Rcpp::clone(mapL);
      for (int i=0; i<M; i++) {
        if(z[candmapL(i,1)-1]==0){
          candmapL(i,1)= -1;
        }
        if(i<Nfixed){
          candmapL(i,0)= -1;
          candmapL(i,1)= -1;
        }
        if(z(i)==0){
          candmapL(i,0)= -1;
        }
      }

      //Find outcandsL and incandsL
      IntegerVector outcandsL=IDs[candmapL(_,1)>0];
      IntegerVector incandsL=IDs[candmapL(_,0)>0];
      int insizeL=incandsL.size();
      int outsizeL=outcandsL.size();

      //Structures to store distances
      NumericVector Dl(incandsL.size());
      NumericVector trashL(insizeL);

      //Build left swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsL(outsizeL+1);
      double x;
      for (int i=0; i<(outsizeL+1); i++) {
        binsL(i)=i;
        x=binsL(i)/outsizeL;
        binsL(i)=x;
      }

      ///////////////// //swap left sides//////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeL+1); i++) {
          if(binsL(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsL(idx-1);

        //find swapout
        swapout=mapL(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeL; i++) {
          Dl(i) = sqrt( pow( s(swapout-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        possible=incandsL[Dl < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsL2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsL2(i)=i;
            x=binsL2(i)/ncand;
            binsL2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsL.size(); i++) {
          trashL(i) = sqrt( pow( s(swapin-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        unpossible=incandsL[trashL < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapL(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_L
        newID=Rcpp::clone(ID_L);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new left data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            ylefttmp(newID(i)-1,j)=0;
            for (int k=0; k<K; k++) {
              ylefttmp(newID(i)-1,j)+=left(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for left data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              prepart(i,j)=yleft(swapped(i)-1,j)*log(pd12(i,j))+(K-yleft(swapped(i)-1,j))*log(1-pd12(i,j));
              postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd12(i,j))+(K-ylefttmp(swapped(i)-1,j))*log(1-pd12(i,j));
            }else{ //if single trap
              prepart(i,j)=yleft(swapped(i)-1,j)*log(pd1(i,j))+(K-yleft(swapped(i)-1,j))*log(1-pd1(i,j));
              postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd1(i,j))+(K-ylefttmp(swapped(i)-1,j))*log(1-pd1(i,j));
            }
            if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
              llpre+=prepart(i,j);
            }
            if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
              llpost+=postpart(i,j);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yleft=Rcpp::clone(ylefttmp);
          ID_L=Rcpp::clone(newID);
          mapL(guy1-1,1)=swapin;
          mapL(guy2-1,1)=swapout;
        }

      }
    }
    //Update Right sides
    if(updates(3)){
      //Build mapR
      IntegerMatrix mapR(M,2);
      mapR(_,1)=ID_R;
      mapR(_,0)=IDs;
      IntegerMatrix candmapR=Rcpp::clone(mapR);
      for (int i=0; i<M; i++) {
        if(z[candmapR(i,1)-1]==0){
          candmapR(i,1)= -1;
        }
        if(i<Nfixed){
          candmapR(i,0)= -1;
          candmapR(i,1)= -1;
        }
        if(z(i)==0){
          candmapR(i,0)= -1;
        }
      }

      //Find outcandsR and incandsR
      IntegerVector outcandsR=IDs[candmapR(_,1)>0];
      IntegerVector incandsR=IDs[candmapR(_,0)>0];
      int insizeR=incandsR.size();
      int outsizeR=outcandsR.size();

      //Structures to store distances
      NumericVector Dr(incandsR.size());
      NumericVector trashR(insizeR);

      //Build right swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsR(outsizeR+1);
      double x;
      for (int i=0; i<(outsizeR+1); i++) {
        binsR(i)=i;
        x=binsR(i)/outsizeR;
        binsR(i)=x;
      }

      ///////////////// //swap right sides/////////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeR+1); i++) {
          if(binsR(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsR(idx-1);

        //find swapout
        swapout=mapR(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeR; i++) {
          Dr(i) = sqrt( pow( s(swapout-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        possible=incandsR[Dr < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsR2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsR2(i)=i;
            x=binsR2(i)/ncand;
            binsR2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsR.size(); i++) {
          trashR(i) = sqrt( pow( s(swapin-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        unpossible=incandsR[trashR < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapR(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_R
        newID=Rcpp::clone(ID_R);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new right data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            yrighttmp(newID(i)-1,j)=0;
            for (int k=0; k<K; k++) {
              yrighttmp(newID(i)-1,j)+=right(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for right side data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++){
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              prepart(i,j)=yright(swapped(i)-1,j)*log(pd12(i,j))+(K-yright(swapped(i)-1,j))*log(1-pd12(i,j));
              postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd12(i,j))+(K-yrighttmp(swapped(i)-1,j))*log(1-pd12(i,j));
            }else{ //if single trap
              prepart(i,j)=yright(swapped(i)-1,j)*log(pd1(i,j))+(K-yright(swapped(i)-1,j))*log(1-pd1(i,j));
              postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd1(i,j))+(K-yrighttmp(swapped(i)-1,j))*log(1-pd1(i,j));
            }
            if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
              llpre+=prepart(i,j);
            }
            if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
              llpost+=postpart(i,j);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yright=Rcpp::clone(yrighttmp);
          ID_R=Rcpp::clone(newID);
          mapR(guy1-1,1)=swapin;
          mapR(guy2-1,1)=swapout;
        }

      } //end swapping rights
    }

    //Recalculate zeroguys
    for (int i=0; i<M; i++) {
      guycounts(i)=0;
      for(int j=0; j<J; j++) {
        guycounts(i)+=yboth(i,j)+yright(i,j)+yleft(i,j);
      }
      zeroguys(i)=guycounts(i)==0;
    }

    //Update Psi
    //  Calculate probability of no capture and update z
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(i,j));
        if(X(j,2)==2){ //if double trap
          pd2(i,j)=1-exp(-lamd2(i,j));
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          pd1traps(i,j)=pd12(i,j);
          pd2traps(i,j)=pd2(i,j);
        }else{ //if single trap
          pd1traps(i,j)=pd1(i,j);
          pd2traps(i,j)=0;
        }
        pbar1(i,j)=pow(1-pd1traps(i,j),2*K);
        pbar2(i,j)=pow(1-pd2traps(i,j),K);
        pbar0(i,j)=pbar1(i,j)*pbar2(i,j);
        prob0(i)*=pbar0(i,j);
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0)&zeroguys(i);
      if(swappable(i)){
        rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
    }
    //Calculate current likelihood
    v=0;
    for(int i=0; i<M; i++) {
      if(z[i]==1){
        for(int j=0; j<J; j++){
          if(X(j,2)==2){ //if double trap
            partboth(i,j)=yboth(i,j)*log(pd2(i,j))+(K-yboth(i,j))*log(1-pd2(i,j));
            partleft(i,j)=yleft(i,j)*log(pd12(i,j))+(K-yleft(i,j))*log(1-pd12(i,j));
            partright(i,j)=yright(i,j)*log(pd12(i,j))+(K-yright(i,j))*log(1-pd12(i,j));
          }else{ //if single trap
            partboth(i,j)=0;
            partleft(i,j)=yleft(i,j)*log(pd1(i,j))+(K-yleft(i,j))*log(1-pd1(i,j));
            partright(i,j)=yright(i,j)*log(pd1(i,j))+(K-yright(i,j))*log(1-pd1(i,j));
          }
          if(partboth(i,j)==partboth(i,j)){
            v+=partboth(i,j);
          }
          if(partleft(i,j)==partleft(i,j)){
            v+=partleft(i,j);
          }
          if(partright(i,j)==partright(i,j)){
            v+=partright(i,j);
          }
        }
      }
    }
    //update psi
    N=sum(z);
    psi=Rcpp::rbeta(1, 1 + N, 1 + M - N);
    likcurr=v;

    //Update Activity Centers
    for(int i=0; i<M; i++) {
      ScandX=Rcpp::rnorm(1,s(i,0),propsx);
      ScandY=Rcpp::rnorm(1,s(i,1),propsy);
      if(useverts==FALSE){
        inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      }else{
        inbox=inoutCpp(ScandX,ScandY,vertices);
      }
      if(inbox(0)){
          dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
          lamd1candc=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
          lamd2candc=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
          if(z[i]==1){ //if in pop, otherwise keep update
            v=0;
            vcand=0;
            //Calculate likelihood for original and candidate
            for(int j=0; j<J; j++){
              pd1c(j)=1-exp(-lamd1(i,j));
              pd1ccand(j)=1-exp(-lamd1candc(j));
              if(X(j,2)==2){ //if double trap
                pd12c(j)=2*pd1c(j)-pd1c(j)*pd1c(j);
                pd12ccand(j)=2*pd1ccand(j)-pd1ccand(j)*pd1ccand(j);
                pd2c(j)=1-exp(-lamd2(i,j));
                pd2ccand(j)=1-exp(-lamd2candc(j));
                partbothc(j)=yboth(i,j)*log(pd2c(j))+(K-yboth(i,j))*log(1-pd2c(j));
                partleftc(j)=yleft(i,j)*log(pd12c(j))+(K-yleft(i,j))*log(1-pd12c(j));
                partrightc(j)=yright(i,j)*log(pd12c(j))+(K-yright(i,j))*log(1-pd12c(j));
                partbothcCand(j)=yboth(i,j)*log(pd2ccand(j))+(K-yboth(i,j))*log(1-pd2ccand(j));
                partleftcCand(j)=yleft(i,j)*log(pd12ccand(j))+(K-yleft(i,j))*log(1-pd12ccand(j));
                partrightcCand(j)=yright(i,j)*log(pd12ccand(j))+(K-yright(i,j))*log(1-pd12ccand(j));
              }else{ //if single trap
                partbothc(j)=0;
                partleftc(j)=yleft(i,j)*log(pd1c(j))+(K-yleft(i,j))*log(1-pd1c(j));
                partrightc(j)=yright(i,j)*log(pd1c(j))+(K-yright(i,j))*log(1-pd1c(j));
                partleftcCand(j)=yleft(i,j)*log(pd1ccand(j))+(K-yleft(i,j))*log(1-pd1ccand(j));
                partrightcCand(j)=yright(i,j)*log(pd1ccand(j))+(K-yright(i,j))*log(1-pd1ccand(j));
              }
              if(partbothc(j)==partbothc(j)){
                v+=partbothc(j);
              }
              if(partleftc(j)==partleftc(j)){
                v+=partleftc(j);
              }
              if(partrightc(j)==partrightc(j)){
                v+=partrightc(j);
              }
              if(partbothcCand(j)==partbothcCand(j)){
                vcand+=partbothcCand(j);
              }
              if(partleftcCand(j)==partleftcCand(j)){
                vcand+=partleftcCand(j);
              }
              if(partrightcCand(j)==partrightcCand(j)){
                vcand+=partrightcCand(j);
              }
            }
            rand=Rcpp::runif(1);
            lldiff=vcand-v;
          }
          if((rand(0)<exp(lldiff))|(z[i]==0)){
            s(i,0)=ScandX(0);
            s(i,1)=ScandY(0);
            for(int j=0; j<J; j++){
              D(i,j) = dtmp(j);
              lamd1(i,j) = lamd1candc(j);
              lamd2(i,j) = lamd2candc(j);
            }
          }
        }
      }

      //Record output
      if(((iter+1)>nburn)&((iter+1) % nthin==0)){
        sxout(iteridx,_)= s(_,0);
        syout(iteridx,_)= s(_,1);
        zout(iteridx,_)= z;
        ID_Lout(iteridx,_)=ID_L;
        ID_Rout(iteridx,_)=ID_R;
        out(iteridx,0)=lam01;
        out(iteridx,1)=lam02;
        out(iteridx,2)=sigma;
        out(iteridx,3)=N;
        iteridx=iteridx+1;
      }
    }

    List to_return(8);
    to_return[0] = out;
    to_return[1] = sxout;
    to_return[2] = syout;
    to_return[3] = ID_Lout;
    to_return[4] = ID_Rout;
    to_return[5] = zout;
    to_return[6] = lamd1;
    to_return[7] = lamd2;
    return to_return;
}

////////////////////////2side model with tf////////////////////////////////
//likelihood calculation
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fullliktf(IntegerVector z,NumericMatrix lamd1,NumericMatrix lamd2,NumericMatrix yboth, NumericMatrix yleft,
               NumericMatrix yright, NumericMatrix X,IntegerVector K) {
  //Preallocate
  int M = z.size();
  int J=X.nrow();
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix partboth(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);
  double v=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    if(z[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(i,j));
        if(X(j,2)==2){ //if double trap
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          pd2(i,j)=1-exp(-lamd2(i,j));
          partboth(i,j)=yboth(i,j)*log(pd2(i,j))+(K(j)-yboth(i,j))*log(1-pd2(i,j));
          partleft(i,j)=yleft(i,j)*log(pd12(i,j))+(K(j)-yleft(i,j))*log(1-pd12(i,j));
          partright(i,j)=yright(i,j)*log(pd12(i,j))+(K(j)-yright(i,j))*log(1-pd12(i,j));
        }else{ //if single trap
          partboth(i,j)=0;
          partleft(i,j)=yleft(i,j)*log(pd1(i,j))+(K(j)-yleft(i,j))*log(1-pd1(i,j));
          partright(i,j)=yright(i,j)*log(pd1(i,j))+(K(j)-yright(i,j))*log(1-pd1(i,j));
        }
        if(partboth(i,j)==partboth(i,j)){
          v+=partboth(i,j);
        }
        if(partleft(i,j)==partleft(i,j)){
          v+=partleft(i,j);
        }
        if(partright(i,j)==partright(i,j)){
          v+=partright(i,j);
        }
      }
    }
  }
  double to_return;
  to_return = v;
  return to_return;
}
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMCtf(double lam01,double lam02, double sigma,NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright, IntegerVector z,
            NumericMatrix X,IntegerVector K,int Kmax,NumericMatrix D,int Nfixed, IntegerVector knownvector,IntegerVector ID_L,
            IntegerVector ID_R,int swap,double swaptol,IntegerVector left,IntegerVector right,NumericMatrix s,NumericVector psi,
            NumericVector xlim,NumericVector ylim,bool useverts, NumericMatrix vertices,double proplam01, double proplam02, double propsigma, double propsx,
            double propsy, int niter, int nburn, int nthin, LogicalVector updates) {
  RNGScope scope;
  int M = z.size();
  int J=yleft.ncol();

  //Preallocate for update lam01 lam02 sigma
  double likcurr;
  double liknew;
  double lam01cand;
  double lam02cand;
  double sigmacand;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd1cand;
  NumericMatrix lamd2cand;
  //Preallocate side swapping structures
  //Swapping structures
  IntegerVector IDs=seq_len(M);
  int guy1=0;
  int guy2=0;
  int swapin;
  int swapout;
  int idx;
  double ncand;
  double ncand2;
  int match;
  double jumpprob;
  double backprob;
  IntegerVector possible;
  IntegerVector unpossible;
  IntegerVector newID(M);
  NumericMatrix ylefttmp(M,J);
  NumericMatrix yrighttmp(M,J);  //could use same structure but dont want to get confused

  //ll structures
  NumericMatrix prepart(2,J);
  NumericMatrix postpart(2,J);
  NumericMatrix pd1b(2,J);
  NumericMatrix pd12b(2,J);
  NumericVector swapped(2);
  double llpre;
  double llpost;

  //MH structures
  double lldiff=0;

  //Housekeeping structures
  NumericVector guycounts(M);
  LogicalVector zeroguys(M);

  //Preallocate for Psi update
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix pd1traps(M,J);
  NumericMatrix pd2traps(M,J);
  NumericMatrix pbar1(M,J);
  NumericMatrix pbar2(M,J);
  NumericMatrix pbar0(M,J);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericMatrix partboth(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);
  int N;

  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector pd1c(J);
  NumericVector pd12c(J);
  NumericVector pd2c(J);
  NumericVector pd1ccand(J);
  NumericVector pd12ccand(J);
  NumericVector pd2ccand(J);
  NumericVector partbothc(J);
  NumericVector partleftc(J);
  NumericVector partrightc(J);
  NumericVector partbothcCand(J);
  NumericVector partleftcCand(J);
  NumericVector partrightcCand(J);
  NumericVector dtmp(J);
  NumericVector lamd1candc(J);
  NumericVector lamd2candc(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double vcand;
  double v;
  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,4);
  NumericMatrix sxout(nstore,M);
  NumericMatrix syout(nstore,M);
  NumericMatrix zout(nstore,M);
  NumericMatrix ID_Lout(nstore,M);
  NumericMatrix ID_Rout(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamd1=calclamd(lam01,sigma,D);
  NumericMatrix lamd2=calclamd(lam02,sigma,D);
  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    //update lam01, lam02, and sigma

    likcurr=fullliktf( z, lamd1, lamd2, yboth,  yleft, yright,  X, K);
    liknew=0;
    lam01cand=0;
    lam02cand=0;
    sigmacand=0;

    //Update lam01
    if(updates(0)){
      rand=Rcpp::rnorm(1,lam01,proplam01);
      if(rand(0) > 0){
        lam01cand=rand(0);
        lamd1cand=calclamd(lam01cand,sigma,D);
        liknew=fullliktf( z, lamd1cand, lamd2, yboth,  yleft, yright,  X, K);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam01=lam01cand;
          lamd1=lamd1cand;
          likcurr=liknew;
        }
      }
    }
    //Update lam02
    if(updates(1)){
      rand=Rcpp::rnorm(1,lam02,proplam02);
      if(rand(0) > 0){
        lam02cand=rand(0);
        lamd2cand=calclamd(lam02cand,sigma,D);
        liknew=fullliktf( z, lamd1, lamd2cand, yboth,  yleft, yright,  X, K);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam02=lam02cand;
          lamd2=lamd2cand;
          likcurr=liknew;
        }
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamd1cand=calclamd(lam01,sigmacand,D);
      lamd2cand=calclamd(lam02,sigmacand,D);
      liknew=fullliktf( z, lamd1cand, lamd2cand, yboth,  yleft, yright,  X, K);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamd1=lamd1cand;
        lamd2=lamd2cand;
        likcurr=liknew;
      }
    }

    //  Update left sides
    if(updates(2)){
      //Build mapL
      IntegerMatrix mapL(M,2);
      mapL(_,1)=ID_L;
      mapL(_,0)=IDs;
      IntegerMatrix candmapL=Rcpp::clone(mapL);
      for (int i=0; i<M; i++) {
        if(z[candmapL(i,1)-1]==0){
          candmapL(i,1)= -1;
        }
        if(i<Nfixed){
          candmapL(i,0)= -1;
          candmapL(i,1)= -1;
        }
        if(z(i)==0){
          candmapL(i,0)= -1;
        }
      }

      //Find outcandsL and incandsL
      IntegerVector outcandsL=IDs[candmapL(_,1)>0];
      IntegerVector incandsL=IDs[candmapL(_,0)>0];
      int insizeL=incandsL.size();
      int outsizeL=outcandsL.size();

      //Structures to store distances
      NumericVector Dl(incandsL.size());
      NumericVector trashL(insizeL);

      //Build left swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsL(outsizeL+1);
      double x;
      for (int i=0; i<(outsizeL+1); i++) {
        binsL(i)=i;
        x=binsL(i)/outsizeL;
        binsL(i)=x;
      }

      ///////////////// //swap left sides//////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeL+1); i++) {
          if(binsL(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsL(idx-1);

        //find swapout
        swapout=mapL(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeL; i++) {
          Dl(i) = sqrt( pow( s(swapout-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        possible=incandsL[Dl < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsL2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsL2(i)=i;
            x=binsL2(i)/ncand;
            binsL2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsL.size(); i++) {
          trashL(i) = sqrt( pow( s(swapin-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        unpossible=incandsL[trashL < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapL(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_L
        newID=Rcpp::clone(ID_L);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new left data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            ylefttmp(newID(i)-1,j)=0;
            for (int k=0; k<Kmax; k++) {
              ylefttmp(newID(i)-1,j)+=left(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for left data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              prepart(i,j)=yleft(swapped(i)-1,j)*log(pd12(i,j))+(K(j)-yleft(swapped(i)-1,j))*log(1-pd12(i,j));
              postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd12(i,j))+(K(j)-ylefttmp(swapped(i)-1,j))*log(1-pd12(i,j));
            }else{ //if single trap
              prepart(i,j)=yleft(swapped(i)-1,j)*log(pd1(i,j))+(K(j)-yleft(swapped(i)-1,j))*log(1-pd1(i,j));
              postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd1(i,j))+(K(j)-ylefttmp(swapped(i)-1,j))*log(1-pd1(i,j));
            }
            if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
              llpre+=prepart(i,j);
            }
            if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
              llpost+=postpart(i,j);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yleft=Rcpp::clone(ylefttmp);
          ID_L=Rcpp::clone(newID);
          mapL(guy1-1,1)=swapin;
          mapL(guy2-1,1)=swapout;
        }

      }
    }
    //Update Right sides
    if(updates(3)){
      //Build mapR
      IntegerMatrix mapR(M,2);
      mapR(_,1)=ID_R;
      mapR(_,0)=IDs;
      IntegerMatrix candmapR=Rcpp::clone(mapR);
      for (int i=0; i<M; i++) {
        if(z[candmapR(i,1)-1]==0){
          candmapR(i,1)= -1;
        }
        if(i<Nfixed){
          candmapR(i,0)= -1;
          candmapR(i,1)= -1;
        }
        if(z(i)==0){
          candmapR(i,0)= -1;
        }
      }

      //Find outcandsR and incandsR
      IntegerVector outcandsR=IDs[candmapR(_,1)>0];
      IntegerVector incandsR=IDs[candmapR(_,0)>0];
      int insizeR=incandsR.size();
      int outsizeR=outcandsR.size();

      //Structures to store distances
      NumericVector Dr(incandsR.size());
      NumericVector trashR(insizeR);

      //Build right swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsR(outsizeR+1);
      double x;
      for (int i=0; i<(outsizeR+1); i++) {
        binsR(i)=i;
        x=binsR(i)/outsizeR;
        binsR(i)=x;
      }

      ///////////////// //swap right sides/////////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeR+1); i++) {
          if(binsR(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsR(idx-1);

        //find swapout
        swapout=mapR(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeR; i++) {
          Dr(i) = sqrt( pow( s(swapout-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        possible=incandsR[Dr < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsR2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsR2(i)=i;
            x=binsR2(i)/ncand;
            binsR2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsR.size(); i++) {
          trashR(i) = sqrt( pow( s(swapin-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        unpossible=incandsR[trashR < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapR(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_R
        newID=Rcpp::clone(ID_R);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new right data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            yrighttmp(newID(i)-1,j)=0;
            for (int k=0; k<Kmax; k++) {
              yrighttmp(newID(i)-1,j)+=right(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for right side data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++){
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              prepart(i,j)=yright(swapped(i)-1,j)*log(pd12(i,j))+(K(j)-yright(swapped(i)-1,j))*log(1-pd12(i,j));
              postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd12(i,j))+(K(j)-yrighttmp(swapped(i)-1,j))*log(1-pd12(i,j));
            }else{ //if single trap
              prepart(i,j)=yright(swapped(i)-1,j)*log(pd1(i,j))+(K(j)-yright(swapped(i)-1,j))*log(1-pd1(i,j));
              postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd1(i,j))+(K(j)-yrighttmp(swapped(i)-1,j))*log(1-pd1(i,j));
            }
            if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
              llpre+=prepart(i,j);
            }
            if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
              llpost+=postpart(i,j);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yright=Rcpp::clone(yrighttmp);
          ID_R=Rcpp::clone(newID);
          mapR(guy1-1,1)=swapin;
          mapR(guy2-1,1)=swapout;
        }

      } //end swapping rights
    }

    //Recalculate zeroguys
    for (int i=0; i<M; i++) {
      guycounts(i)=0;
      for(int j=0; j<J; j++) {
        guycounts(i)+=yboth(i,j)+yright(i,j)+yleft(i,j);
      }
      zeroguys(i)=guycounts(i)==0;
    }

    //Update Psi
    //  Calculate probability of no capture and update z
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(i,j));
        if(X(j,2)==2){ //if double trap
          pd2(i,j)=1-exp(-lamd2(i,j));
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          pd1traps(i,j)=pd12(i,j);
          pd2traps(i,j)=pd2(i,j);
        }else{ //if single trap
          pd1traps(i,j)=pd1(i,j);
          pd2traps(i,j)=0;
        }
        pbar1(i,j)=pow(1-pd1traps(i,j),2*K(j));
        pbar2(i,j)=pow(1-pd2traps(i,j),K(j));
        pbar0(i,j)=pbar1(i,j)*pbar2(i,j);
        prob0(i)*=pbar0(i,j);
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0)&zeroguys(i);
      if(swappable(i)){
        rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
    }
    //Calculate current likelihood
    v=0;
    for(int i=0; i<M; i++) {
      if(z[i]==1){
        for(int j=0; j<J; j++){
          if(X(j,2)==2){ //if double trap
            partboth(i,j)=yboth(i,j)*log(pd2(i,j))+(K(j)-yboth(i,j))*log(1-pd2(i,j));
            partleft(i,j)=yleft(i,j)*log(pd12(i,j))+(K(j)-yleft(i,j))*log(1-pd12(i,j));
            partright(i,j)=yright(i,j)*log(pd12(i,j))+(K(j)-yright(i,j))*log(1-pd12(i,j));
          }else{ //if single trap
            partboth(i,j)=0;
            partleft(i,j)=yleft(i,j)*log(pd1(i,j))+(K(j)-yleft(i,j))*log(1-pd1(i,j));
            partright(i,j)=yright(i,j)*log(pd1(i,j))+(K(j)-yright(i,j))*log(1-pd1(i,j));
          }
          if(partboth(i,j)==partboth(i,j)){
            v+=partboth(i,j);
          }
          if(partleft(i,j)==partleft(i,j)){
            v+=partleft(i,j);
          }
          if(partright(i,j)==partright(i,j)){
            v+=partright(i,j);
          }
        }
      }
    }
    //update psi
    N=sum(z);
    psi=Rcpp::rbeta(1, 1 + N, 1 + M - N);
    likcurr=v;

    //Update Activity Centers
    for(int i=0; i<M; i++) {
      ScandX=Rcpp::rnorm(1,s(i,0),propsx);
      ScandY=Rcpp::rnorm(1,s(i,1),propsy);
      if(useverts==FALSE){
        inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      }else{
        inbox=inoutCpp(ScandX,ScandY,vertices);
      }
      if(inbox(0)){
        dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
        lamd1candc=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
        lamd2candc=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
        if(z[i]==1){ //if in pop, otherwise keep update
          v=0;
          vcand=0;
          //Calculate likelihood for original and candidate
          for(int j=0; j<J; j++){
            pd1c(j)=1-exp(-lamd1(i,j));
            pd1ccand(j)=1-exp(-lamd1candc(j));
            if(X(j,2)==2){ //if double trap
              pd12c(j)=2*pd1c(j)-pd1c(j)*pd1c(j);
              pd12ccand(j)=2*pd1ccand(j)-pd1ccand(j)*pd1ccand(j);
              pd2c(j)=1-exp(-lamd2(i,j));
              pd2ccand(j)=1-exp(-lamd2candc(j));
              partbothc(j)=yboth(i,j)*log(pd2c(j))+(K(j)-yboth(i,j))*log(1-pd2c(j));
              partleftc(j)=yleft(i,j)*log(pd12c(j))+(K(j)-yleft(i,j))*log(1-pd12c(j));
              partrightc(j)=yright(i,j)*log(pd12c(j))+(K(j)-yright(i,j))*log(1-pd12c(j));
              partbothcCand(j)=yboth(i,j)*log(pd2ccand(j))+(K(j)-yboth(i,j))*log(1-pd2ccand(j));
              partleftcCand(j)=yleft(i,j)*log(pd12ccand(j))+(K(j)-yleft(i,j))*log(1-pd12ccand(j));
              partrightcCand(j)=yright(i,j)*log(pd12ccand(j))+(K(j)-yright(i,j))*log(1-pd12ccand(j));
            }else{ //if single trap
              partbothc(j)=0;
              partleftc(j)=yleft(i,j)*log(pd1c(j))+(K(j)-yleft(i,j))*log(1-pd1c(j));
              partrightc(j)=yright(i,j)*log(pd1c(j))+(K(j)-yright(i,j))*log(1-pd1c(j));
              partleftcCand(j)=yleft(i,j)*log(pd1ccand(j))+(K(j)-yleft(i,j))*log(1-pd1ccand(j));
              partrightcCand(j)=yright(i,j)*log(pd1ccand(j))+(K(j)-yright(i,j))*log(1-pd1ccand(j));
            }
            if(partbothc(j)==partbothc(j)){
              v+=partbothc(j);
            }
            if(partleftc(j)==partleftc(j)){
              v+=partleftc(j);
            }
            if(partrightc(j)==partrightc(j)){
              v+=partrightc(j);
            }
            if(partbothcCand(j)==partbothcCand(j)){
              vcand+=partbothcCand(j);
            }
            if(partleftcCand(j)==partleftcCand(j)){
              vcand+=partleftcCand(j);
            }
            if(partrightcCand(j)==partrightcCand(j)){
              vcand+=partrightcCand(j);
            }
          }
          rand=Rcpp::runif(1);
          lldiff=vcand-v;
        }
        if((rand(0)<exp(lldiff))|(z[i]==0)){
          s(i,0)=ScandX(0);
          s(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D(i,j) = dtmp(j);
            lamd1(i,j) = lamd1candc(j);
            lamd2(i,j) = lamd2candc(j);
          }
        }
      }
    }

    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      sxout(iteridx,_)= s(_,0);
      syout(iteridx,_)= s(_,1);
      zout(iteridx,_)= z;
      ID_Lout(iteridx,_)=ID_L;
      ID_Rout(iteridx,_)=ID_R;
      out(iteridx,0)=lam01;
      out(iteridx,1)=lam02;
      out(iteridx,2)=sigma;
      out(iteridx,3)=N;
      iteridx=iteridx+1;
    }
  }

  List to_return(6);
  to_return[0] = out;
  to_return[1] = sxout;
  to_return[2] = syout;
  to_return[3] = ID_Lout;
  to_return[4] = ID_Rout;
  to_return[5] = zout;
  return to_return;
}

///////////////////////2side model full tf functionality//////////////////////
//likelihood calculation

//[[Rcpp::export]]
double fullliktf2(IntegerVector z,NumericMatrix lamd1,NumericMatrix lamd2,arma::cube ybothC, arma::cube yleftC,
                  arma::cube yrightC, NumericMatrix X,IntegerMatrix tf) {
  //Preallocate
  int M = z.size();
  int J=X.nrow();
  int K=tf.ncol();
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  arma::cube partboth(M,J,K);
  arma::cube partleft(M,J,K);
  arma::cube partright(M,J,K);
  double v=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    if(z[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(i,j));
        if(X(j,2)==2){//if double trap
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          pd2(i,j)=1-exp(-lamd2(i,j));
          for(int k=0; k<K; k++){
            if(tf(j,k)==2){ //if two cams on
              partboth(i,j,k)=ybothC(i,j,k)*log(pd2(i,j))+(1-ybothC(i,j,k))*log(1-pd2(i,j));
              partleft(i,j,k)=yleftC(i,j,k)*log(pd12(i,j))+(1-yleftC(i,j,k))*log(1-pd12(i,j));
              partright(i,j,k)=yrightC(i,j,k)*log(pd12(i,j))+(1-yrightC(i,j,k))*log(1-pd12(i,j));
            }else if(tf(j,k)==1){//if one cam on
              partboth(i,j,k)=0;
              partleft(i,j,k)=yleftC(i,j,k)*log(pd1(i,j))+(1-yleftC(i,j,k))*log(1-pd1(i,j));
              partright(i,j,k)=yrightC(i,j,k)*log(pd1(i,j))+(1-yrightC(i,j,k))*log(1-pd1(i,j));
            }else{//if no cams on
              partboth(i,j,k)=0;
              partleft(i,j,k)=0;
              partright(i,j,k)=0;
            }
            if(partboth(i,j,k)==partboth(i,j,k)){
              v+=partboth(i,j,k);
            }
            if(partleft(i,j,k)==partleft(i,j,k)){
              v+=partleft(i,j,k);
            }
            if(partright(i,j,k)==partright(i,j,k)){
              v+=partright(i,j,k);
            }
          }
        }else{ //if single trap
          for(int k=0; k<K; k++){
            if(tf(j,k)==1){//if 1 cam on
              partboth(i,j,k)=0;
              partleft(i,j,k)=yleftC(i,j,k)*log(pd1(i,j))+(1-yleftC(i,j,k))*log(1-pd1(i,j));
              partright(i,j,k)=yrightC(i,j,k)*log(pd1(i,j))+(1-yrightC(i,j,k))*log(1-pd1(i,j));
            }else{//if no cams on
              partboth(i,j,k)=0;
              partleft(i,j,k)=0;
              partright(i,j,k)=0;
            }
            if(partboth(i,j,k)==partboth(i,j,k)){
              v+=partboth(i,j,k);
            }
            if(partleft(i,j,k)==partleft(i,j,k)){
              v+=partleft(i,j,k);
            }
            if(partright(i,j,k)==partright(i,j,k)){
              v+=partright(i,j,k);
            }
          }
        }
      }
    }
  }
  double to_return;
  to_return = v;
  return to_return;
}



using namespace Rcpp;
// [[Rcpp::export]]
List MCMCtf2(double lam01,double lam02, double sigma,NumericVector yboth, NumericVector yleft, NumericVector yright, IntegerVector z,
             NumericMatrix X,IntegerMatrix tf,NumericMatrix D,int Nfixed, IntegerVector knownvector,IntegerVector ID_L,
             IntegerVector ID_R,int swap,double swaptol,IntegerVector left,IntegerVector right,NumericMatrix s,NumericVector psi,
             NumericVector xlim,NumericVector ylim,bool useverts, NumericMatrix vertices,double proplam01, double proplam02, double propsigma, double propsx,
             double propsy, int niter, int nburn, int nthin, LogicalVector updates) {
  RNGScope scope;
  int M = z.size();
  int J=tf.nrow();
  int K=tf.ncol();

  //put data in cubes
  IntegerVector arrayDims = yboth.attr("dim");
  arma::cube ybothC(yboth.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  arma::cube yleftC(yleft.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  arma::cube yrightC(yright.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  //Preallocate for update lam01 lam02 sigma
  double likcurr;
  double liknew;
  double lam01cand;
  double lam02cand;
  double sigmacand;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd1cand;
  NumericMatrix lamd2cand;
  //Preallocate side swapping structures
  //Swapping structures
  IntegerVector IDs=seq_len(M);
  int guy1=0;
  int guy2=0;
  int swapin;
  int swapout;
  int idx;
  double ncand;
  double ncand2;
  int match;
  double jumpprob;
  double backprob;
  IntegerVector possible;
  IntegerVector unpossible;
  IntegerVector newID(M);
  arma::cube ylefttmp(M,J,K);
  arma::cube yrighttmp(M,J,K);  //could use same structure but dont want to get confused

  //ll structures
  arma::cube prepart(2,J,K);
  arma::cube postpart(2,J,K);
  arma::cube pd1b(2,J,K);
  arma::cube pd12b(2,J,K);
  NumericVector swapped(2);
  double llpre;
  double llpost;

  //MH structures
  double lldiff=0;

  //Housekeeping structures
  NumericVector guycounts(M);
  LogicalVector zeroguys(M);

  //Preallocate for Psi update
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  arma::cube pd1traps(M,J,K);
  arma::cube pd2traps(M,J,K);
  arma::cube pbar1(M,J,K);
  arma::cube pbar2(M,J,K);
  arma::cube pbar0(M,J,K);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  arma::cube partboth(M,J,K);
  arma::cube partleft(M,J,K);
  arma::cube partright(M,J,K);
  int N;

  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector pd1c(J);
  NumericVector pd12c(J);
  NumericVector pd2c(J);
  NumericVector pd1ccand(J);
  NumericVector pd12ccand(J);
  NumericVector pd2ccand(J);
  NumericMatrix partbothc(J,K);
  NumericMatrix partleftc(J,K);
  NumericMatrix partrightc(J,K);
  NumericMatrix partbothcCand(J,K);
  NumericMatrix partleftcCand(J,K);
  NumericMatrix partrightcCand(J,K);
  NumericVector dtmp(J);
  NumericVector lamd1candc(J);
  NumericVector lamd2candc(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double vcand;
  double v;
  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,4);
  NumericMatrix sxout(nstore,M);
  NumericMatrix syout(nstore,M);
  NumericMatrix zout(nstore,M);
  NumericMatrix ID_Lout(nstore,M);
  NumericMatrix ID_Rout(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamd1=calclamd(lam01,sigma,D);
  NumericMatrix lamd2=calclamd(lam02,sigma,D);
  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    //update lam01, lam02, and sigma

    likcurr=fullliktf2( z, lamd1, lamd2, ybothC,  yleftC, yrightC,  X, tf);
    //Update lam01
    if(updates(0)){
      rand=Rcpp::rnorm(1,lam01,proplam01);
      if(rand(0) > 0){
        lam01cand=rand(0);
        lamd1cand=calclamd(lam01cand,sigma,D);
        liknew=fullliktf2( z, lamd1cand, lamd2, ybothC,  yleftC, yrightC,  X, tf);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam01=lam01cand;
          lamd1=lamd1cand;
          likcurr=liknew;
        }
      }
    }
    //Update lam02
    if(updates(1)){
      rand=Rcpp::rnorm(1,lam02,proplam02);
      if(rand(0) > 0){
        lam02cand=rand(0);
        lamd2cand=calclamd(lam02cand,sigma,D);
        liknew=fullliktf2( z, lamd1, lamd2cand, ybothC,  yleftC, yrightC,  X, tf);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam02=lam02cand;
          lamd2=lamd2cand;
          likcurr=liknew;
        }
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamd1cand=calclamd(lam01,sigmacand,D);
      lamd2cand=calclamd(lam02,sigmacand,D);
      liknew=fullliktf2( z, lamd1cand, lamd2cand, ybothC,  yleftC, yrightC,  X, tf);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamd1=lamd1cand;
        lamd2=lamd2cand;
        likcurr=liknew;
      }
    }

    //  Update left sides
    if(updates(2)){
      //Build mapL
      IntegerMatrix mapL(M,2);
      mapL(_,1)=ID_L;
      mapL(_,0)=IDs;
      IntegerMatrix candmapL=Rcpp::clone(mapL);
      for (int i=0; i<M; i++) {
        if(z[candmapL(i,1)-1]==0){
          candmapL(i,1)= -1;
        }
        if(i<Nfixed){
          candmapL(i,0)= -1;
          candmapL(i,1)= -1;
        }
        if(z(i)==0){
          candmapL(i,0)= -1;
        }
      }

      //Find outcandsL and incandsL
      IntegerVector outcandsL=IDs[candmapL(_,1)>0];
      IntegerVector incandsL=IDs[candmapL(_,0)>0];
      int insizeL=incandsL.size();
      int outsizeL=outcandsL.size();

      //Structures to store distances
      NumericVector Dl(incandsL.size());
      NumericVector trashL(insizeL);

      //Build left swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsL(outsizeL+1);
      double x;
      for (int i=0; i<(outsizeL+1); i++) {
        binsL(i)=i;
        x=binsL(i)/outsizeL;
        binsL(i)=x;
      }

      ///////////////// //swap left sides//////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeL+1); i++) {
          if(binsL(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsL(idx-1);

        //find swapout
        swapout=mapL(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeL; i++) {
          Dl(i) = sqrt( pow( s(swapout-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        possible=incandsL[Dl < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsL2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsL2(i)=i;
            x=binsL2(i)/ncand;
            binsL2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsL.size(); i++) {
          trashL(i) = sqrt( pow( s(swapin-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        unpossible=incandsL[trashL < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapL(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_L
        newID=Rcpp::clone(ID_L);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new left data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            for (int k=0; k<K; k++) {
              ylefttmp(newID(i)-1,j,k)=left(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for left data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              for(int k=0; k<K; k++){
                if(tf(j,k)==2){ //if two cams on
                  prepart(i,j,k)=yleftC(swapped(i)-1,j,k)*log(pd12(i,j))+(1-yleftC(swapped(i)-1,j,k))*log(1-pd12(i,j));
                  postpart(i,j,k)=ylefttmp(swapped(i)-1,j,k)*log(pd12(i,j))+(1-ylefttmp(swapped(i)-1,j,k))*log(1-pd12(i,j));
                }else if(tf(j,k)==1){//if one cam on
                  prepart(i,j,k)=yleftC(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yleftC(swapped(i)-1,j,k))*log(1-pd1(i,j));
                  postpart(i,j,k)=ylefttmp(swapped(i)-1,j,k)*log(pd1(i,j))+(1-ylefttmp(swapped(i)-1,j,k))*log(1-pd1(i,j));
                }else{//if no cam on
                  prepart(i,j,k)=0;
                  postpart(i,j,k)=0;
                }
                if((prepart(i,j,k)==prepart(i,j,k))&(!std::isinf(prepart(i,j,k)))){
                  llpre+=prepart(i,j,k);
                }
                if((postpart(i,j,k)==postpart(i,j,k))&(!std::isinf(postpart(i,j,k)))){
                  llpost+=postpart(i,j,k);
                }
              }
            }else{ //if single trap
              for(int k=0; k<K; k++){
                if(tf(j,k)==1){ //if one cam on
                  prepart(i,j,k)=yleftC(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yleftC(swapped(i)-1,j,k))*log(1-pd1(i,j));
                  postpart(i,j,k)=ylefttmp(swapped(i)-1,j,k)*log(pd1(i,j))+(1-ylefttmp(swapped(i)-1,j,k))*log(1-pd1(i,j));
                }else{//if no cam on
                  prepart(i,j,k)=yleftC(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yleftC(swapped(i)-1,j,k))*log(1-pd1(i,j));
                  postpart(i,j,k)=ylefttmp(swapped(i)-1,j,k)*log(pd1(i,j))+(1-ylefttmp(swapped(i)-1,j,k))*log(1-pd1(i,j));
                }
                if((prepart(i,j,k)==prepart(i,j,k))&(!std::isinf(prepart(i,j,k)))){
                  llpre+=prepart(i,j,k);
                }
                if((postpart(i,j,k)==postpart(i,j,k))&(!std::isinf(postpart(i,j,k)))){
                  llpost+=postpart(i,j,k);
                }
              }
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yleftC=ylefttmp;
          ID_L=Rcpp::clone(newID);
          mapL(guy1-1,1)=swapin;
          mapL(guy2-1,1)=swapout;
        }

      }
    }
    //Update Right sides
    if(updates(3)){
      //Build mapR
      IntegerMatrix mapR(M,2);
      mapR(_,1)=ID_R;
      mapR(_,0)=IDs;
      IntegerMatrix candmapR=Rcpp::clone(mapR);
      for (int i=0; i<M; i++) {
        if(z[candmapR(i,1)-1]==0){
          candmapR(i,1)= -1;
        }
        if(i<Nfixed){
          candmapR(i,0)= -1;
          candmapR(i,1)= -1;
        }
        if(z(i)==0){
          candmapR(i,0)= -1;
        }
      }

      //Find outcandsR and incandsR
      IntegerVector outcandsR=IDs[candmapR(_,1)>0];
      IntegerVector incandsR=IDs[candmapR(_,0)>0];
      int insizeR=incandsR.size();
      int outsizeR=outcandsR.size();

      //Structures to store distances
      NumericVector Dr(incandsR.size());
      NumericVector trashR(insizeR);

      //Build right swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsR(outsizeR+1);
      double x;
      for (int i=0; i<(outsizeR+1); i++) {
        binsR(i)=i;
        x=binsR(i)/outsizeR;
        binsR(i)=x;
      }

      ///////////////// //swap right sides/////////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeR+1); i++) {
          if(binsR(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsR(idx-1);

        //find swapout
        swapout=mapR(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeR; i++) {
          Dr(i) = sqrt( pow( s(swapout-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        possible=incandsR[Dr < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsR2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsR2(i)=i;
            x=binsR2(i)/ncand;
            binsR2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsR.size(); i++) {
          trashR(i) = sqrt( pow( s(swapin-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        unpossible=incandsR[trashR < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapR(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_R
        newID=Rcpp::clone(ID_R);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new right data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            for (int k=0; k<K; k++) {
              yrighttmp(newID(i)-1,j,k)=right(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for right side data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              for(int k=0; k<K; k++){
                if(tf(j,k)==2){ //if two cams on
                  prepart(i,j,k)=yrightC(swapped(i)-1,j,k)*log(pd12(i,j))+(1-yrightC(swapped(i)-1,j,k))*log(1-pd12(i,j));
                  postpart(i,j,k)=yrighttmp(swapped(i)-1,j,k)*log(pd12(i,j))+(1-yrighttmp(swapped(i)-1,j,k))*log(1-pd12(i,j));
                }else if(tf(j,k)==1){//if one cam on
                  prepart(i,j,k)=yrightC(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yrightC(swapped(i)-1,j,k))*log(1-pd1(i,j));
                  postpart(i,j,k)=yrighttmp(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yrighttmp(swapped(i)-1,j,k))*log(1-pd1(i,j));
                }else{//if no cam on
                  prepart(i,j,k)=0;
                  postpart(i,j,k)=0;
                }
                if((prepart(i,j,k)==prepart(i,j,k))&(!std::isinf(prepart(i,j,k)))){
                  llpre+=prepart(i,j,k);
                }
                if((postpart(i,j,k)==postpart(i,j,k))&(!std::isinf(postpart(i,j,k)))){
                  llpost+=postpart(i,j,k);
                }
              }
            }else{ //if single trap
              for(int k=0; k<K; k++){
                if(tf(j,k)==1){ //if one cam on
                  prepart(i,j,k)=yrightC(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yrightC(swapped(i)-1,j,k))*log(1-pd1(i,j));
                  postpart(i,j,k)=yrighttmp(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yrighttmp(swapped(i)-1,j,k))*log(1-pd1(i,j));
                }else{//if no cam on
                  prepart(i,j,k)=yrightC(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yrightC(swapped(i)-1,j,k))*log(1-pd1(i,j));
                  postpart(i,j,k)=yrighttmp(swapped(i)-1,j,k)*log(pd1(i,j))+(1-yrighttmp(swapped(i)-1,j,k))*log(1-pd1(i,j));
                }
                if((prepart(i,j,k)==prepart(i,j,k))&(!std::isinf(prepart(i,j,k)))){
                  llpre+=prepart(i,j,k);
                }
                if((postpart(i,j,k)==postpart(i,j,k))&(!std::isinf(postpart(i,j,k)))){
                  llpost+=postpart(i,j,k);
                }
              }
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yrightC=yrighttmp;
          ID_R=Rcpp::clone(newID);
          mapR(guy1-1,1)=swapin;
          mapR(guy2-1,1)=swapout;
        }

      } //end swapping rights
    }

    //Recalculate zeroguys
    for (int i=0; i<M; i++) {
      guycounts(i)=0;
      for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
          guycounts(i)+=ybothC(i,j,k)+yrightC(i,j,k)+yleftC(i,j,k);
        }
      }
      zeroguys(i)=guycounts(i)==0;
    }

    //Update Psi
    //  Calculate probability of no capture and update z
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(i,j));
        if(X(j,2)==2){ //if double trap
          pd2(i,j)=1-exp(-lamd2(i,j));
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          for(int k=0; k<K; k++){
            if(tf(j,k)==2){ //if two cams on
              pd1traps(i,j,k)=pd12(i,j);
              pd2traps(i,j,k)=pd2(i,j);
            }else if(tf(j,k)==1){//if one cam on
              pd1traps(i,j,k)=pd1(i,j);
              pd2traps(i,j,k)=0;
            }else{//if no cams on
              pd1traps(i,j,k)=0;
              pd2traps(i,j,k)=0;
            }
            pbar1(i,j,k)=pow(1-pd1traps(i,j,k),2);
            pbar2(i,j,k)=1-pd2traps(i,j,k);
            pbar0(i,j,k)=pbar1(i,j,k)*pbar2(i,j,k);
            prob0(i)*=pbar0(i,j,k);
          }
        }else{ //if single trap
          for(int k=0; k<K; k++){
            if(tf(j,k)==1){//if 1 cam on
              pd1traps(i,j,k)=pd1(i,j);
              pd2traps(i,j,k)=0;
            }else{//if no cams on
              pd1traps(i,j,k)=0;
              pd2traps(i,j,k)=0;
            }
            pbar1(i,j,k)=pow(1-pd1traps(i,j,k),2);
            //pbar2(i,j,k)=1;
            pbar0(i,j,k)=pbar1(i,j,k);
            prob0(i)*=pbar0(i,j,k);
          }
        }
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0)&zeroguys(i);
      if(swappable(i)){
        rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
    }
    N=sum(z);
    psi=Rcpp::rbeta(1, 1 + N, 1 + M - N);

    //Update Activity Centers
    for(int i=0; i<M; i++) {
      ScandX=Rcpp::rnorm(1,s(i,0),propsx);
      ScandY=Rcpp::rnorm(1,s(i,1),propsy);
      if(useverts==FALSE){
        inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      }else{
        inbox=inoutCpp(ScandX,ScandY,vertices);
      }
      if(inbox(0)){
        dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
        lamd1candc=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
        lamd2candc=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
        v=0;
        vcand=0;
        if(z[i]==1){ //if in pop, otherwise keep update
          //Calculate likelihood for original and candidate
          for(int j=0; j<J; j++){
            pd1c(j)=1-exp(-lamd1(i,j));
            pd1ccand(j)=1-exp(-lamd1candc(j));
            if(X(j,2)==2){ //if double trap
              pd12c(j)=2*pd1c(j)-pd1c(j)*pd1c(j);
              pd12ccand(j)=2*pd1ccand(j)-pd1ccand(j)*pd1ccand(j);
              pd2c(j)=1-exp(-lamd2(i,j));
              pd2ccand(j)=1-exp(-lamd2candc(j));
              for(int k=0; k<K; k++){
                if(tf(j,k)==2){ //if two cams on
                  partbothc(j,k)=ybothC(i,j,k)*log(pd2c(j))+(1-ybothC(i,j,k))*log(1-pd2c(j));
                  partleftc(j,k)=yleftC(i,j,k)*log(pd12c(j))+(1-yleftC(i,j,k))*log(1-pd12c(j));
                  partrightc(j,k)=yrightC(i,j,k)*log(pd12c(j))+(1-yrightC(i,j,k))*log(1-pd12c(j));
                  partbothcCand(j,k)=ybothC(i,j,k)*log(pd2ccand(j))+(1-ybothC(i,j,k))*log(1-pd2ccand(j));
                  partleftcCand(j,k)=yleftC(i,j,k)*log(pd12ccand(j))+(1-yleftC(i,j,k))*log(1-pd12ccand(j));
                  partrightcCand(j,k)=yrightC(i,j,k)*log(pd12ccand(j))+(1-yrightC(i,j,k))*log(1-pd12ccand(j));
                }else if(tf(j,k)==1){//if one cam on
                  partbothc(j,k)=0;
                  partleftc(j,k)=yleftC(i,j,k)*log(pd1c(j))+(1-yleftC(i,j,k))*log(1-pd1c(j));
                  partrightc(j,k)=yrightC(i,j,k)*log(pd1c(j))+(1-yrightC(i,j,k))*log(1-pd1c(j));
                  partbothcCand(j,k)=0;
                  partleftcCand(j,k)=yleftC(i,j,k)*log(pd1ccand(j))+(1-yleftC(i,j,k))*log(1-pd1ccand(j));
                  partrightcCand(j,k)=yrightC(i,j,k)*log(pd1ccand(j))+(1-yrightC(i,j,k))*log(1-pd1ccand(j));
                }else{//if no cams on
                  partbothc(j,k)=0;
                  partleftc(j,k)=0;
                  partrightc(j,k)=0;
                  partbothcCand(j,k)=0;
                  partleftcCand(j,k)=0;
                  partrightcCand(j,k)=0;
                }
                if(partbothc(j,k)==partbothc(j,k)){
                  v+=partbothc(j,k);
                }
                if(partleftc(j,k)==partleftc(j,k)){
                  v+=partleftc(j,k);
                }
                if(partrightc(j,k)==partrightc(j,k)){
                  v+=partrightc(j,k);
                }
                if(partbothcCand(j,k)==partbothcCand(j,k)){
                  vcand+=partbothcCand(j,k);
                }
                if(partleftcCand(j,k)==partleftcCand(j,k)){
                  vcand+=partleftcCand(j,k);
                }
                if(partrightcCand(j,k)==partrightcCand(j,k)){
                  vcand+=partrightcCand(j,k);
                }
              }
            }else{ //if single trap
              for(int k=0; k<K; k++){
                if(tf(j,k)==1){//if 1 cam on
                  partbothc(j,k)=0;
                  partleftc(j,k)=yleftC(i,j,k)*log(pd1c(j))+(1-yleftC(i,j,k))*log(1-pd1c(j));
                  partrightc(j,k)=yrightC(i,j,k)*log(pd1c(j))+(1-yrightC(i,j,k))*log(1-pd1c(j));
                  partleftcCand(j,k)=yleftC(i,j,k)*log(pd1ccand(j))+(1-yleftC(i,j,k))*log(1-pd1ccand(j));
                  partrightcCand(j,k)=yrightC(i,j,k)*log(pd1ccand(j))+(1-yrightC(i,j,k))*log(1-pd1ccand(j));
                }else{//no cam on
                  partbothc(j,k)=0;
                  partleftc(j,k)=0;
                  partrightc(j,k)=0;
                  partleftcCand(j,k)=0;
                  partrightcCand(j,k)=0;
                }
                if(partbothc(j,k)==partbothc(j,k)){
                  v+=partbothc(j,k);
                }
                if(partleftc(j,k)==partleftc(j,k)){
                  v+=partleftc(j,k);
                }
                if(partrightc(j,k)==partrightc(j,k)){
                  v+=partrightc(j,k);
                }
                if(partbothcCand(j,k)==partbothcCand(j,k)){
                  vcand+=partbothcCand(j,k);
                }
                if(partleftcCand(j,k)==partleftcCand(j,k)){
                  vcand+=partleftcCand(j,k);
                }
                if(partrightcCand(j,k)==partrightcCand(j,k)){
                  vcand+=partrightcCand(j,k);
                }
              }
            }
          }
        }
        rand=Rcpp::runif(1);
        lldiff=vcand-v;
        if((rand(0)<exp(lldiff))|(z[i]==0)){
          s(i,0)=ScandX(0);
          s(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D(i,j) = dtmp(j);
            lamd1(i,j) = lamd1candc(j);
            lamd2(i,j) = lamd2candc(j);
          }
        }
      }
    }

    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      sxout(iteridx,_)= s(_,0);
      syout(iteridx,_)= s(_,1);
      zout(iteridx,_)= z;
      ID_Lout(iteridx,_)=ID_L;
      ID_Rout(iteridx,_)=ID_R;
      out(iteridx,0)=lam01;
      out(iteridx,1)=lam02;
      out(iteridx,2)=sigma;
      out(iteridx,3)=N;
      iteridx=iteridx+1;
    }
  }

  List to_return(6);
  to_return[0] = out;
  to_return[1] = sxout;
  to_return[2] = syout;
  to_return[3] = ID_Lout;
  to_return[4] = ID_Rout;
  to_return[5] = zout;
  return to_return;
}
/////////////////////Heuristic model////////////////////////
//likelihood calculation
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fulllikInd(IntegerVector z1,IntegerVector z2,NumericMatrix lamd11,NumericMatrix lamd21,NumericMatrix lamd12,
                  NumericMatrix lamd22,NumericMatrix yboth, NumericMatrix yleft,NumericMatrix yright, NumericMatrix X,int K,
                  int Nknown) {
  //Preallocate
  int M = z1.size();
  int J=X.nrow();
  NumericMatrix pd11(M,J);
  NumericMatrix pd11b(M,J);
  NumericMatrix pd21(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd12b(M,J);
  NumericMatrix pd22(M,J);
  NumericMatrix partboth1(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);
  NumericMatrix partboth2(M,J);

  double v=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    //Left side data set first
    if(z1[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pd11(i,j)=1-exp(-lamd11(i,j));
        pd21(i,j)=1-exp(-lamd21(i,j));
        if(X(j,2)==2){ //if double trap
          pd11b(i,j)=2*pd11(i,j)-pd11(i,j)*pd11(i,j);
          pd21(i,j)=1-exp(-lamd21(i,j));
          partboth1(i,j)=yboth(i,j)*log(pd21(i,j))+(K-yboth(i,j))*log(1-pd21(i,j));
          partleft(i,j)=yleft(i,j)*log(pd11b(i,j))+(K-yleft(i,j))*log(1-pd11b(i,j));
        }else{ //if single trap
          partboth1(i,j)=0;
          partleft(i,j)=yleft(i,j)*log(pd11(i,j))+(K-yleft(i,j))*log(1-pd11(i,j));
        }
        if(partboth1(i,j)==partboth1(i,j)){
          v+=partboth1(i,j);
        }
        if(partleft(i,j)==partleft(i,j)){
          v+=partleft(i,j);
        }
      }
    }
    //Now do right side data set
    if(z2[i]==1){ //if in right side pop
      for(int j=0; j<J; j++){
        pd12(i,j)=1-exp(-lamd12(i,j));
        pd22(i,j)=1-exp(-lamd22(i,j));
        if(X(j,2)==2){ //if double trap
          pd12b(i,j)=2*pd12(i,j)-pd12(i,j)*pd12(i,j);
          if(i>Nknown){ //don't count both guys twice
            pd22(i,j)=1-exp(-lamd22(i,j));
            partboth2(i,j)=yboth(i,j)*log(pd22(i,j))+(K-yboth(i,j))*log(1-pd22(i,j));
          }else{
            partboth2(i,j)=0;
          }
          partright(i,j)=yright(i,j)*log(pd12b(i,j))+(K-yright(i,j))*log(1-pd12b(i,j));
        }else{ //if single trap
          partboth2(i,j)=0;
          partright(i,j)=yright(i,j)*log(pd12(i,j))+(K-yright(i,j))*log(1-pd12(i,j));
        }
        if(partboth2(i,j)==partboth2(i,j)){
          v+=partboth2(i,j);
        }
        if(partright(i,j)==partright(i,j)){
          v+=partright(i,j);
        }
      }

    }
  }
  double to_return;
  to_return = v;
  return to_return;
}
//likelihood calculation2
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fulllikInd1D(int z1,int z2,int nth,int Nfixed,NumericVector lamd11,NumericVector lamd21,NumericVector lamd12,
                    NumericVector lamd22,NumericVector yboth, NumericVector yleft,NumericVector yright, NumericMatrix X,int K) {
  //Preallocate
  int J=X.nrow();
  NumericVector pd11(J);
  NumericVector pd11b(J);
  NumericVector pd21(J);
  NumericVector pd12(J);
  NumericVector pd12b(J);
  NumericVector pd22(J);
  NumericVector partboth1(J);
  NumericVector partleft(J);
  NumericVector partright(J);
  NumericVector partboth2(J);

  double v=0;
  //  Calculate likelihood
  //Left side data set first
  if(z1==1){ //if in pop
    for(int j=0; j<J; j++){
      pd11(j)=1-exp(-lamd11(j));
      pd21(j)=1-exp(-lamd21(j));
      if(X(j,2)==2){ //if double trap
        pd11b(j)=2*pd11(j)-pd11(j)*pd11(j);
        pd21(j)=1-exp(-lamd21(j));
        partboth1(j)=yboth(j)*log(pd21(j))+(K-yboth(j))*log(1-pd21(j));
        partleft(j)=yleft(j)*log(pd11b(j))+(K-yleft(j))*log(1-pd11b(j));
      }else{ //if single trap
        partboth1(j)=0;
        partleft(j)=yleft(j)*log(pd11(j))+(K-yleft(j))*log(1-pd11(j));
      }
      if(partboth1(j)==partboth1(j)){
        v+=partboth1(j);
      }
      if(partleft(j)==partleft(j)){
        v+=partleft(j);
      }
    }
  }
  //Now do right side data set
  if(z2==1){ //if in right side pop
    for(int j=0; j<J; j++){
      pd12(j)=1-exp(-lamd12(j));
      if(X(j,2)==2){ //if double trap
        pd12b(j)=2*pd12(j)-pd12(j)*pd12(j);
        if(nth>Nfixed){ //don't count both guys twice
          pd22(j)=1-exp(-lamd22(j));
          partboth2(j)=yboth(j)*log(pd22(j))+(K-yboth(j))*log(1-pd22(j));
        }else{
          partboth2(j)=0;
        }
        partright(j)=yright(j)*log(pd12b(j))+(K-yright(j))*log(1-pd12b(j));
      }else{ //if single trap
        partboth2(j)=0;
        partright(j)=yright(j)*log(pd12(j))+(K-yright(j))*log(1-pd12(j));
      }
      if(partboth2(j)==partboth2(j)){
        v+=partboth2(j);
      }
      if(partright(j)==partright(j)){
        v+=partright(j);
      }
    }

  }
  double to_return;
  to_return = v;
  return to_return;
}


//update lam01, lam02, and sigma

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateparmsInd(double lam01,double lam02, double sigma, NumericMatrix lamd11,NumericMatrix lamd21,NumericMatrix lamd12,
                    NumericMatrix lamd22,NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright, IntegerVector z1,
                    IntegerVector z2,NumericMatrix X,int K,NumericMatrix D1,NumericMatrix D2,int Nfixed) {
  RNGScope scope;
  //Preallocate
  //int M = z.size();
  //int J=X.nrow();
  double likcurr=fulllikInd(z1,z2,lamd11,lamd21,lamd12,lamd22,yboth, yleft,yright,X,K,Nfixed);
  double liknew=0;
  double lam01cand=0;
  double lam02cand=0;
  double sigmacand=0;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd11cand;
  NumericMatrix lamd21cand;
  NumericMatrix lamd12cand;
  NumericMatrix lamd22cand;
  double explldiff;

  //Update lam01
  rand=Rcpp::rnorm(1,lam01,0.05);
  if(rand(0) > 0){
    lam01cand=rand(0);
    lamd11cand=calclamd(lam01cand,sigma,D1);
    lamd12cand=calclamd(lam01cand,sigma,D2);
    liknew=fulllikInd(z1,z2,lamd11cand,lamd21,lamd12cand,lamd22,yboth, yleft,yright,X,K,Nfixed);
    rand2=Rcpp::runif(1);
    explldiff = exp(liknew-likcurr);
    if(rand2(0)<explldiff){
      lam01=lam01cand;
      lamd11=lamd11cand;
      lamd12=lamd12cand;
      likcurr=liknew;
    }
  }
  //Update lam02
  rand=Rcpp::rnorm(1,lam02,0.05);
  if(rand(0) > 0){
    lam02cand=rand(0);
    lamd21cand=calclamd(lam02cand,sigma,D1);
    lamd22cand=calclamd(lam02cand,sigma,D2);
    liknew=fulllikInd(z1,z2,lamd11,lamd21cand,lamd12,lamd22cand,yboth, yleft,yright,X,K,Nfixed);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(liknew-likcurr)){
      lam02=lam02cand;
      lamd21=lamd21cand;
      lamd22=lamd22cand;
      likcurr=liknew;
    }
  }
  //Update sigma
  rand=Rcpp::rnorm(1,sigma,0.1);
  if(rand(0) > 0){
    sigmacand=rand(0);
    lamd11cand=calclamd(lam01,sigmacand,D1);
    lamd21cand=calclamd(lam02,sigmacand,D1);
    lamd12cand=calclamd(lam01,sigmacand,D2);
    lamd22cand=calclamd(lam02,sigmacand,D2);
    liknew=fulllikInd(z1,z2,lamd11cand,lamd21cand,lamd12cand,lamd22cand,yboth, yleft,yright,X,K,Nfixed);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(liknew-likcurr)){
      sigma=sigmacand;
      lamd11=lamd11cand;
      lamd21=lamd21cand;
      lamd12=lamd12cand;
      lamd22=lamd22cand;
      likcurr=liknew;
    }
  }
  List to_return(8);
  to_return[0] = lam01;
  to_return[1] = lam02;
  to_return[2] = sigma;
  to_return[3] = lamd11;
  to_return[4] = lamd21;
  to_return[5] = lamd12;
  to_return[6] = lamd22;
  to_return[7] = likcurr;
  return to_return;
}
//Update psi,z
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updatePsiInd(IntegerVector z1, IntegerVector z2,IntegerVector capvector1, IntegerVector capvector2,NumericMatrix lamd11,
                  NumericMatrix lamd12,NumericMatrix lamd21,NumericMatrix lamd22, NumericMatrix yboth,NumericMatrix yleft,
                  NumericMatrix yright,NumericMatrix X, int K,double psi1,double psi2,int Nfixed) {
  RNGScope scope;
  //Preallocate
  int M = z1.size();
  int J=yleft.ncol();
  NumericMatrix pd11(M,J);
  NumericMatrix pd11b(M,J);
  NumericMatrix pd21(M,J);
  NumericMatrix pd11traps(M,J);
  NumericMatrix pd21traps(M,J);
  NumericMatrix pbar11(M,J);
  NumericMatrix pbar21(M,J);
  NumericMatrix pbar01(M,J);
  NumericVector prob01(M);
  NumericVector fc1(M);
  NumericMatrix pd12(M,J);
  NumericMatrix pd12b(M,J);
  NumericMatrix pd22(M,J);
  NumericMatrix pd12traps(M,J);
  NumericMatrix pd22traps(M,J);
  NumericMatrix pbar12(M,J);
  NumericMatrix pbar22(M,J);
  NumericMatrix pbar02(M,J);
  NumericVector prob02(M);
  NumericVector fc2(M);
  //  Calculate probability of no capture and update z1
  for(int i=0; i<M; i++) {
    prob01(i)=1;
    for(int j=0; j<J; j++){
      pd11(i,j)=1-exp(-lamd11(i,j));
      if(X(j,2)==2){ //if double trap
        pd21(i,j)=1-exp(-lamd21(i,j));
        pd11b(i,j)=2*pd11(i,j)-pd11(i,j)*pd11(i,j);
        pd11traps(i,j)=pd11b(i,j);
        pd21traps(i,j)=pd21(i,j);
      }else{ //if single trap
        pd11traps(i,j)=pd11(i,j);
        pd21traps(i,j)=0;
      }
      pbar11(i,j)=pow(1-pd11traps(i,j),K);
      pbar21(i,j)=pow(1-pd21traps(i,j),K);
      pbar01(i,j)=pbar11(i,j)*pbar21(i,j);
      prob01(i)*=pbar01(i,j);
    }
    fc1(i)=prob01(i)*psi1/(prob01(i)*psi1 + 1-psi1);
    if(capvector1(i)==0){
      NumericVector rand=Rcpp::rbinom(1,1,fc1(i));
      z1(i)=rand(0);
    }
  }
  //  Calculate probability of no capture and update z2
  for(int i=0; i<M; i++) {
    prob02(i)=1;
    for(int j=0; j<J; j++){
      pd12(i,j)=1-exp(-lamd12(i,j));
      if(X(j,2)==2){ //if double trap
        pd22(i,j)=1-exp(-lamd22(i,j));
        pd12b(i,j)=2*pd12(i,j)-pd12(i,j)*pd12(i,j);
        pd12traps(i,j)=pd12b(i,j);
        pd22traps(i,j)=pd22(i,j);
      }else{ //if single trap
        pd12traps(i,j)=pd12(i,j);
        pd22traps(i,j)=0;
      }
      pbar12(i,j)=pow(1-pd12traps(i,j),K);
      pbar22(i,j)=pow(1-pd22traps(i,j),K);
      pbar02(i,j)=pbar12(i,j)*pbar22(i,j);
      prob02(i)*=pbar02(i,j);
    }
    fc2(i)=prob02(i)*psi2/(prob02(i)*psi2 + 1-psi2);
    if(capvector2(i)==0){
      NumericVector rand=Rcpp::rbinom(1,1,fc2(i));
      z2(i)=rand(0);
    }
  }
  //Calculate current likelihood (no constants)
  double v=fulllikInd(z1,z2,lamd11,lamd21,lamd12,lamd22,yboth, yleft,yright,X,K,Nfixed);
  //update psis
  NumericVector psi1new=Rcpp::rbeta(1, 1 + sum(z1), 1 + M - sum(z1));
  NumericVector psi2new=Rcpp::rbeta(1, 1 + sum(z2), 1 + M - sum(z2));
  List to_return(5);
  to_return[0] = psi1new;
  to_return[1] = psi2new;
  to_return[2] = v;
  to_return[3] = z1;
  to_return[4] = z2;
  return to_return;
}


//Update activity centers
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateACsInd(IntegerVector z1,IntegerVector z2, NumericMatrix s1, NumericMatrix s2, IntegerVector xlim,IntegerVector ylim,
                  bool useverts, NumericMatrix vertices,
                  NumericMatrix D1,NumericMatrix D2, NumericMatrix lamd11,NumericMatrix lamd21, NumericMatrix lamd12,NumericMatrix lamd22,
                  double lam01, double lam02,double sigma, NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright,
                  NumericMatrix X,int K,int Nfixed) {
  RNGScope scope;
  //Preallocate
  //capture process
  int M = z1.size();
  int J=X.nrow();
  LogicalVector inbox(1);
  NumericVector rand(1);
  double explldiff;
  NumericVector dtmp(J);
  NumericVector lamd11cand(J);
  NumericVector lamd21cand(J);
  NumericVector lamd12cand(J);
  NumericVector lamd22cand(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double v;
  double vcand;
  NumericMatrix s1new=Rcpp::clone(s1);
  NumericMatrix D1new=Rcpp::clone(D1);
  NumericMatrix s2new=Rcpp::clone(s2);
  NumericMatrix D2new=Rcpp::clone(D2);
  NumericMatrix lamd11new=Rcpp::clone(lamd11);
  NumericMatrix lamd21new=Rcpp::clone(lamd21);
  NumericMatrix lamd12new=Rcpp::clone(lamd12);
  NumericMatrix lamd22new=Rcpp::clone(lamd22);

  //  Update both activity centers
  if(Nfixed>0){
    for(int i=0; i<(Nfixed-1); i++) {
      ScandX=Rcpp::rnorm(1,s1(i,0),.2);
      ScandY=Rcpp::rnorm(1,s1(i,1),.2);
      if(useverts==FALSE){
        inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      }else{
        inbox=inoutCpp(ScandX,ScandY,vertices);
      }
      if(inbox(0)){
          dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
          lamd11cand=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
          lamd21cand=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
          //Calculate likelihood for original and candidate
          v=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11(i,_), lamd21(i,_), lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
          vcand=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11cand,lamd21cand, lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
          rand=Rcpp::runif(1);
          explldiff=exp(vcand-v);
          if((rand(0)<explldiff)){
            s1new(i,0)=ScandX(0);
            s1new(i,1)=ScandY(0);
            s2new(i,0)=ScandX(0);
            s2new(i,1)=ScandY(0);
            for(int j=0; j<J; j++){
              D1new(i,j) = dtmp(j);
              D2new(i,j) = dtmp(j);
              lamd11new(i,j) = lamd11cand(j);
              lamd21new(i,j) = lamd21cand(j);
              lamd12new(i,j) = lamd11cand(j);
              lamd22new(i,j) = lamd21cand(j);
            }
          }
        }
      }
    }
    //Update Left
    for(int i=Nfixed; i<(M-1); i++) {
      if(z1[i]==1){
        ScandX=Rcpp::rnorm(1,s1(i,0),.2);
        ScandY=Rcpp::rnorm(1,s1(i,1),.2);
        if(useverts==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
            dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
            lamd11cand=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
            lamd21cand=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
            //Calculate likelihood for original and candidate
            v=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11(i,_), lamd21(i,_), lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
            vcand=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11cand,lamd21cand, lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
            rand=Rcpp::runif(1);
            explldiff=exp(vcand-v);
            if((rand(0)<explldiff)){
              s1new(i,0)=ScandX(0);
              s1new(i,1)=ScandY(0);
              for(int j=0; j<J; j++){
                D1new(i,j) = dtmp(j);
                lamd11new(i,j) = lamd11cand(j);
                lamd21new(i,j) = lamd21cand(j);
              }
            }
          }
        }
    }
    //Update Right
    for(int i=Nfixed; i<(M-1); i++) {
      if(z2[i]==1){
        ScandX=Rcpp::rnorm(1,s2(i,0),.2);
        ScandY=Rcpp::rnorm(1,s2(i,1),.2);
        if(useverts==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
            dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
            lamd12cand=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
            lamd22cand=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
            //Calculate likelihood for original and candidate
            v=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11(i,_), lamd21(i,_), lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
            vcand=fulllikInd1D( z1(i), z2(i), i, Nfixed,lamd11(i,_), lamd21(i,_), lamd12cand,lamd22cand, yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
            rand=Rcpp::runif(1);
            explldiff=exp(vcand-v);
            if((rand(0)<explldiff)){
              s2new(i,0)=ScandX(0);
              s2new(i,1)=ScandY(0);
              for(int j=0; j<J; j++){
                D2new(i,j) = dtmp(j);
                lamd12new(i,j) = lamd12cand(j);
                lamd22new(i,j) = lamd22cand(j);
              }
            }
          }
        }
    }

    List to_return(8);
    to_return[0] = s1new;
    to_return[1] = s2new;
    to_return[2] = D1new;
    to_return[3] = D2new;
    to_return[4] = lamd11new;
    to_return[5] = lamd21new;
    to_return[6] = lamd12new;
    to_return[7] = lamd22new;
    return to_return;
}


//WholeShebang
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMCInd(double lam01,double lam02, double sigma,
             NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright, IntegerVector z1,IntegerVector z2,
             NumericMatrix X,int K,NumericMatrix D1,NumericMatrix D2,int Nfixed, IntegerVector capvector1,
             IntegerVector capvector2,NumericMatrix s1,NumericMatrix s2,NumericVector psi1, NumericVector psi2, NumericVector xlim,
             NumericVector ylim,LogicalVector useverts, NumericMatrix vertices,double proplam01, double proplam02,double propsigma, double propsx,
             double propsy, int niter, int nburn, int nthin, LogicalVector updates) {
  RNGScope scope;
  int M = z1.size();
  int J=X.nrow();
  //Preallocate for update lam01 lam02 sigma
  //ll structures
  double likcurr=0;
  double liknew=0;
  double lam01cand=0;
  double lam02cand=0;
  double sigmacand=0;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd11cand;
  NumericMatrix lamd21cand;
  NumericMatrix lamd12cand;
  NumericMatrix lamd22cand;
  double explldiff;

  //Preallocate for Psi update
  NumericMatrix pd11(M,J);
  NumericMatrix pd11b(M,J);
  NumericMatrix pd21(M,J);
  NumericMatrix pd11traps(M,J);
  NumericMatrix pd21traps(M,J);
  NumericMatrix pbar11(M,J);
  NumericMatrix pbar21(M,J);
  NumericMatrix pbar01(M,J);
  NumericVector prob01(M);
  NumericVector fc1(M);
  NumericMatrix pd12(M,J);
  NumericMatrix pd12b(M,J);
  NumericMatrix pd22(M,J);
  NumericMatrix pd12traps(M,J);
  NumericMatrix pd22traps(M,J);
  NumericMatrix pbar12(M,J);
  NumericMatrix pbar22(M,J);
  NumericMatrix pbar02(M,J);
  NumericVector prob02(M);
  NumericVector fc2(M);
  int N;

  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector dtmp(J);
  NumericVector lamd11cand1d(J);
  NumericVector lamd21cand1d(J);
  NumericVector lamd12cand1d(J);
  NumericVector lamd22cand1d(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double vcand;
  double v;

  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,4);
  NumericMatrix s1xout(nstore,M);
  NumericMatrix s1yout(nstore,M);
  NumericMatrix z1out(nstore,M);
  NumericMatrix s2xout(nstore,M);
  NumericMatrix s2yout(nstore,M);
  NumericMatrix z2out(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamd11=calclamd(lam01,sigma,D1);
  NumericMatrix lamd12=calclamd(lam02,sigma,D1);
  NumericMatrix lamd21=calclamd(lam01,sigma,D2);
  NumericMatrix lamd22=calclamd(lam02,sigma,D2);
  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    ////////update lam01, lam02, and sigma
    likcurr=fulllikInd(z1,z2,lamd11,lamd21,lamd12,lamd22,yboth, yleft,yright,X,K,Nfixed);
    //Update lam01
    if(updates(0)){
      rand=Rcpp::rnorm(1,lam01,proplam01);
      if(rand(0) > 0){
        lam01cand=rand(0);
        lamd11cand=calclamd(lam01cand,sigma,D1);
        lamd12cand=calclamd(lam01cand,sigma,D2);
        liknew=fulllikInd(z1,z2,lamd11cand,lamd21,lamd12cand,lamd22,yboth, yleft,yright,X,K,Nfixed);
        rand2=Rcpp::runif(1);
        explldiff = exp(liknew-likcurr);
        if(rand2(0)<explldiff){
          lam01=lam01cand;
          lamd11=lamd11cand;
          lamd12=lamd12cand;
          likcurr=liknew;
        }
      }
    }
    //Update lam02
    if(updates(1)){
      rand=Rcpp::rnorm(1,lam02,proplam02);
      if(rand(0) > 0){
        lam02cand=rand(0);
        lamd21cand=calclamd(lam02cand,sigma,D1);
        lamd22cand=calclamd(lam02cand,sigma,D2);
        liknew=fulllikInd(z1,z2,lamd11,lamd21cand,lamd12,lamd22cand,yboth, yleft,yright,X,K,Nfixed);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam02=lam02cand;
          lamd21=lamd21cand;
          lamd22=lamd22cand;
          likcurr=liknew;
        }
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamd11cand=calclamd(lam01,sigmacand,D1);
      lamd21cand=calclamd(lam02,sigmacand,D1);
      lamd12cand=calclamd(lam01,sigmacand,D2);
      lamd22cand=calclamd(lam02,sigmacand,D2);
      liknew=fulllikInd(z1,z2,lamd11cand,lamd21cand,lamd12cand,lamd22cand,yboth, yleft,yright,X,K,Nfixed);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamd11=lamd11cand;
        lamd21=lamd21cand;
        lamd12=lamd12cand;
        lamd22=lamd22cand;
        likcurr=liknew;
      }
    }
    ///////update psi, z
    //  Calculate probability of no capture and update z1
    for(int i=0; i<M; i++) {
      prob01(i)=1;
      for(int j=0; j<J; j++){
        pd11(i,j)=1-exp(-lamd11(i,j));
        if(X(j,2)==2){ //if double trap
          pd21(i,j)=1-exp(-lamd21(i,j));
          pd11b(i,j)=2*pd11(i,j)-pd11(i,j)*pd11(i,j);
          pd11traps(i,j)=pd11b(i,j);
          pd21traps(i,j)=pd21(i,j);
        }else{ //if single trap
          pd11traps(i,j)=pd11(i,j);
          pd21traps(i,j)=0;
        }
        pbar11(i,j)=pow(1-pd11traps(i,j),K);
        pbar21(i,j)=pow(1-pd21traps(i,j),K);
        pbar01(i,j)=pbar11(i,j)*pbar21(i,j);
        prob01(i)*=pbar01(i,j);
      }
      fc1(i)=prob01(i)*psi1(0)/(prob01(i)*psi1(0) + 1-psi1(0));
      if(capvector1(i)==0){
        NumericVector rand=Rcpp::rbinom(1,1,fc1(i));
        z1(i)=rand(0);
      }
    }
    //  Calculate probability of no capture and update z2
    for(int i=0; i<M; i++) {
      prob02(i)=1;
      for(int j=0; j<J; j++){
        pd12(i,j)=1-exp(-lamd12(i,j));
        if(X(j,2)==2){ //if double trap
          pd22(i,j)=1-exp(-lamd22(i,j));
          pd12b(i,j)=2*pd12(i,j)-pd12(i,j)*pd12(i,j);
          pd12traps(i,j)=pd12b(i,j);
          pd22traps(i,j)=pd22(i,j);
        }else{ //if single trap
          pd12traps(i,j)=pd12(i,j);
          pd22traps(i,j)=0;
        }
        pbar12(i,j)=pow(1-pd12traps(i,j),K);
        pbar22(i,j)=pow(1-pd22traps(i,j),K);
        pbar02(i,j)=pbar12(i,j)*pbar22(i,j);
        prob02(i)*=pbar02(i,j);
      }
      fc2(i)=prob02(i)*psi2(0)/(prob02(i)*psi2(0) + 1-psi2(0));
      if(capvector2(i)==0){
        NumericVector rand=Rcpp::rbinom(1,1,fc2(i));
        z2(i)=rand(0);
      }
    }
    //Calculate current likelihood
    likcurr=fulllikInd(z1,z2,lamd11,lamd21,lamd12,lamd22,yboth, yleft,yright,X,K,Nfixed);
    //update psis
    psi1=Rcpp::rbeta(1, 1 + sum(z1), 1 + M - sum(z1));
    psi2=Rcpp::rbeta(1, 1 + sum(z2), 1 + M - sum(z2));
    N=(sum(z1)+sum(z2))/2;

    //////////Update ACs
    //  Update both activity centers
    if(Nfixed>0){
      for(int i=0; i<(Nfixed-1); i++) {
        ScandX=Rcpp::rnorm(1,s1(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s1(i,1),propsy);
        if(useverts(0)==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
            dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
            lamd11cand1d=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
            lamd21cand1d=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
            //Calculate likelihood for original and candidate
            v=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11(i,_), lamd21(i,_), lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
            vcand=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11cand1d,lamd21cand1d, lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
            rand=Rcpp::runif(1);
            explldiff=exp(vcand-v);
            if((rand(0)<explldiff)){
              s1(i,0)=ScandX(0);
              s1(i,1)=ScandY(0);
              s2(i,0)=ScandX(0);
              s2(i,1)=ScandY(0);
              for(int j=0; j<J; j++){
                D1(i,j) = dtmp(j);
                D2(i,j) = dtmp(j);
                lamd11(i,j) = lamd11cand1d(j);
                lamd21(i,j) = lamd21cand1d(j);
                lamd12(i,j) = lamd11cand1d(j);
                lamd22(i,j) = lamd21cand1d(j);
              }
            }
          }
        }
      }
      //Update Left
      for(int i=Nfixed; i<(M-1); i++) {
        if(z1[i]==1){
          ScandX=Rcpp::rnorm(1,s1(i,0),propsx);
          ScandY=Rcpp::rnorm(1,s1(i,1),propsy);
          if(useverts(0)==FALSE){
            inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
          }else{
            inbox=inoutCpp(ScandX,ScandY,vertices);
          }
          if(inbox(0)){
              dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
              lamd11cand1d=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
              lamd21cand1d=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
              //Calculate likelihood for original and candidate
              v=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11(i,_), lamd21(i,_), lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
              vcand=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11cand1d,lamd21cand1d, lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
              rand=Rcpp::runif(1);
              explldiff=exp(vcand-v);
              if((rand(0)<explldiff)){
                s1(i,0)=ScandX(0);
                s1(i,1)=ScandY(0);
                for(int j=0; j<J; j++){
                  D1(i,j) = dtmp(j);
                  lamd11(i,j) = lamd11cand1d(j);
                  lamd21(i,j) = lamd21cand1d(j);
                }
              }
            }
          }
        }
        //Update Right
        for(int i=Nfixed; i<(M-1); i++) {
          ScandX=Rcpp::rnorm(1,s2(i,0),propsx);
          ScandY=Rcpp::rnorm(1,s2(i,1),propsy);
          if(useverts(0)==FALSE){
            inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
          }else{
            inbox=inoutCpp(ScandX,ScandY,vertices);
          }
          if(inbox(0)){
                dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
                lamd12cand1d=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
                lamd22cand1d=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
                //Calculate likelihood for original and candidate
                v=fulllikInd1D( z1(i), z2(i), i, Nfixed, lamd11(i,_), lamd21(i,_), lamd12(i,_),lamd22(i,_), yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
                vcand=fulllikInd1D( z1(i), z2(i), i, Nfixed,lamd11(i,_), lamd21(i,_), lamd12cand1d,lamd22cand1d, yboth(i,_),  yleft(i,_), yright(i,_),  X, K);
                rand=Rcpp::runif(1);
                explldiff=exp(vcand-v);
                if((rand(0)<explldiff)){
                  s2(i,0)=ScandX(0);
                  s2(i,1)=ScandY(0);
                  for(int j=0; j<J; j++){
                    D2(i,j) = dtmp(j);
                    lamd12(i,j) = lamd12cand1d(j);
                    lamd22(i,j) = lamd22cand1d(j);
                  }
                }
              }
            }
        //Record output
        if(((iter+1)>nburn)&((iter+1) % nthin==0)){
          s1xout(iteridx,_)= s1(_,0);
          s1yout(iteridx,_)= s1(_,1);
          s2xout(iteridx,_)= s2(_,0);
          s2yout(iteridx,_)= s2(_,1);
          z1out(iteridx,_)= z1;
          z2out(iteridx,_)= z2;
          out(iteridx,0)=lam01;
          out(iteridx,1)=lam02;
          out(iteridx,2)=sigma;
          out(iteridx,3)=N;
          iteridx=iteridx+1;
        }
    }
    List to_return(7);
    to_return[0] = out;
    to_return[1] = s1xout;
    to_return[2] = s1yout;
    to_return[3] = s2xout;
    to_return[4] = s2yout;
    to_return[5] = z1out;
    to_return[6] = z2out;
    return to_return;
}
/////////////////////Heuristic estimator V2//////////////////////////

//likelihood calculation
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fulllikInd2(IntegerVector z1,IntegerVector z2,NumericMatrix lamdL,NumericMatrix lamdR,
                   NumericMatrix yleft,NumericMatrix yright, NumericMatrix X,int K) {
  //Preallocate
  int M = z1.size();
  int J=X.nrow();
  NumericMatrix pdL(M,J);
  NumericMatrix pdLb(M,J);
  NumericMatrix pdR(M,J);
  NumericMatrix pdRb(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);

  double v=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    //Left side data set first
    if(z1[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pdL(i,j)=1-exp(-lamdL(i,j));
        if(X(j,2)==2){ //if double trap
          pdLb(i,j)=2*pdL(i,j)-pdL(i,j)*pdL(i,j);
          partleft(i,j)=yleft(i,j)*log(pdLb(i,j))+(K-yleft(i,j))*log(1-pdLb(i,j));
        }else{ //if single trap
          partleft(i,j)=yleft(i,j)*log(pdL(i,j))+(K-yleft(i,j))*log(1-pdL(i,j));
        }
        if(partleft(i,j)==partleft(i,j)){
          v+=partleft(i,j);
        }
      }
    }
    //Now do right side
    if(z2[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pdR(i,j)=1-exp(-lamdR(i,j));
        if(X(j,2)==2){ //if double trap
          pdRb(i,j)=2*pdR(i,j)-pdR(i,j)*pdR(i,j);
          partright(i,j)=yright(i,j)*log(pdRb(i,j))+(K-yright(i,j))*log(1-pdRb(i,j));
        }else{ //if single trap
          partright(i,j)=yright(i,j)*log(pdR(i,j))+(K-yright(i,j))*log(1-pdR(i,j));
        }
        if(partright(i,j)==partright(i,j)){
          v+=partright(i,j);
        }
      }
    }
  }
  double to_return;
  to_return = v;
  return to_return;
}
//likelihood calculation2
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fulllikInd1D2(int z1,int z2,NumericVector lamdL,NumericVector lamdR,NumericVector yleft,NumericVector yright, NumericMatrix X,int K) {
  //Preallocate
  int J=X.nrow();
  NumericVector pdL(J);
  NumericVector pdLb(J);
  NumericVector pdR(J);
  NumericVector pdRb(J);
  NumericVector partleft(J);
  NumericVector partright(J);
  double v=0;
  //  Calculate likelihood
  //Left side data set first
  if(z1==1){ //if in pop
    for(int j=0; j<J; j++){
      pdL(j)=1-exp(-lamdL(j));
      if(X(j,2)==2){ //if double trap
        pdLb(j)=2*pdL(j)-pdL(j)*pdL(j);
        partleft(j)=yleft(j)*log(pdLb(j))+(K-yleft(j))*log(1-pdLb(j));
      }else{ //if single trap
        partleft(j)=yleft(j)*log(pdL(j))+(K-yleft(j))*log(1-pdL(j));
      }
      if(partleft(j)==partleft(j)){
        v+=partleft(j);
      }
    }
  }
  //Now do right side data set
  if(z2==1){ //if in right side pop
    for(int j=0; j<J; j++){
      pdR(j)=1-exp(-lamdR(j));
      if(X(j,2)==2){ //if double trap
        pdRb(j)=2*pdR(j)-pdR(j)*pdR(j);
        partright(j)=yright(j)*log(pdRb(j))+(K-yright(j))*log(1-pdRb(j));
      }else{ //if single trap
        partright(j)=yright(j)*log(pdR(j))+(K-yright(j))*log(1-pdR(j));
      }
      if(partright(j)==partright(j)){
        v+=partright(j);
      }
    }
  }
  double to_return;
  to_return = v;
  return to_return;
}

//update lam01 and sigma

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateparmsInd2(double lam01,double sigma, NumericMatrix lamdL,NumericMatrix lamdR,NumericMatrix yleft,
                     NumericMatrix yright, IntegerVector z1,IntegerVector z2,NumericMatrix X,int K,NumericMatrix D1,
                     NumericMatrix D2) {
  RNGScope scope;
  //Preallocate
  //int M = z.size();
  //int J=X.nrow();
  double likcurr=fulllikInd2(z1,z2,lamdL,lamdR,yleft,yright,X,K);
  double liknew=0;
  double lam01cand=0;
  double sigmacand=0;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamdLcand;
  NumericMatrix lamdRcand;
  double explldiff;

  //Update lam01
  rand=Rcpp::rnorm(1,lam01,0.05);
  if(rand(0) > 0){
    lam01cand=rand(0);
    lamdLcand=calclamd(lam01cand,sigma,D1);
    lamdRcand=calclamd(lam01cand,sigma,D2);
    liknew=fulllikInd2(z1,z2,lamdLcand,lamdRcand, yleft,yright,X,K);
    rand2=Rcpp::runif(1);
    explldiff = exp(liknew-likcurr);
    if(rand2(0)<explldiff){
      lam01=lam01cand;
      lamdL=lamdLcand;
      lamdR=lamdRcand;
      likcurr=liknew;
    }
  }
  //Update sigma
  rand=Rcpp::rnorm(1,sigma,0.1);
  if(rand(0) > 0){
    sigmacand=rand(0);
    lamdLcand=calclamd(lam01,sigmacand,D1);
    lamdRcand=calclamd(lam01,sigmacand,D2);
    liknew=fulllikInd2(z1,z2,lamdLcand,lamdRcand,yleft,yright,X,K);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(liknew-likcurr)){
      sigma=sigmacand;
      lamdL=lamdLcand;
      lamdR=lamdRcand;
      likcurr=liknew;
    }
  }
  List to_return(5);
  to_return[0] = lam01;
  to_return[1] = sigma;
  to_return[2] = lamdL;
  to_return[3] = lamdR;
  to_return[4] = likcurr;
  return to_return;
}
//Update psi,z
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updatePsiInd2(IntegerVector z1, IntegerVector z2,IntegerVector capvector1, IntegerVector capvector2,NumericMatrix lamdL,
                   NumericMatrix lamdR,NumericMatrix yleft,NumericMatrix yright,NumericMatrix X, int K,double psi1,double psi2) {
  RNGScope scope;
  //Preallocate
  int M = z1.size();
  int J=yleft.ncol();
  NumericMatrix pdL(M,J);
  NumericMatrix pdLb(M,J);
  NumericMatrix pdLtraps(M,J);
  NumericMatrix pbar0L(M,J);
  NumericVector prob0L(M);
  NumericVector fc1(M);
  NumericMatrix pdR(M,J);
  NumericMatrix pdRb(M,J);
  NumericMatrix pdRtraps(M,J);
  NumericMatrix pbar0R(M,J);
  NumericVector prob0R(M);
  NumericVector fc2(M);
  //  Calculate probability of no capture and update z1
  for(int i=0; i<M; i++) {
    prob0L(i)=1;
    for(int j=0; j<J; j++){
      pdL(i,j)=1-exp(-lamdL(i,j));
      if(X(j,2)==2){ //if double trap
        pdLb(i,j)=2*pdL(i,j)-pdL(i,j)*pdL(i,j);
        pdLtraps(i,j)=pdLb(i,j);
      }else{ //if single trap
        pdLtraps(i,j)=pdL(i,j);
      }
      pbar0L(i,j)=pow(1-pdLtraps(i,j),K);
      prob0L(i)*=pbar0L(i,j);
    }
    fc1(i)=prob0L(i)*psi1/(prob0L(i)*psi1 + 1-psi1);
    if(capvector1(i)==0){
      NumericVector rand=Rcpp::rbinom(1,1,fc1(i));
      z1(i)=rand(0);
    }
  }
  //  Calculate probability of no capture and update z2
  for(int i=0; i<M; i++) {
    prob0R(i)=1;
    for(int j=0; j<J; j++){
      pdR(i,j)=1-exp(-lamdR(i,j));
      if(X(j,2)==2){ //if double trap
        pdRb(i,j)=2*pdR(i,j)-pdR(i,j)*pdR(i,j);
        pdRtraps(i,j)=pdRb(i,j);
      }else{ //if single trap
        pdRtraps(i,j)=pdR(i,j);
      }
      pbar0R(i,j)=pow(1-pdRtraps(i,j),K);
      prob0R(i)*=pbar0R(i,j);
    }
    fc2(i)=prob0R(i)*psi2/(prob0R(i)*psi2 + 1-psi2);
    if(capvector2(i)==0){
      NumericVector rand=Rcpp::rbinom(1,1,fc2(i));
      z2(i)=rand(0);
    }
  }
  //Calculate current likelihood (no constants)
  double v=fulllikInd2(z1,z2,lamdL,lamdR,yleft,yright,X,K);
  //update psis
  NumericVector psi1new=Rcpp::rbeta(1, 1 + sum(z1), 1 + M - sum(z1));
  NumericVector psi2new=Rcpp::rbeta(1, 1 + sum(z2), 1 + M - sum(z2));
  List to_return(5);
  to_return[0] = psi1new;
  to_return[1] = psi2new;
  to_return[2] = v;
  to_return[3] = z1;
  to_return[4] = z2;
  return to_return;
}


//Update activity centers
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List updateACsInd2(IntegerVector z1,IntegerVector z2, NumericMatrix s1, NumericMatrix s2, IntegerVector xlim,IntegerVector ylim,
                   NumericMatrix D1,NumericMatrix D2, NumericMatrix lamdL,NumericMatrix lamdR,double lam01,double sigma,
                   NumericMatrix yleft, NumericMatrix yright,NumericMatrix X,int K) {
  RNGScope scope;
  //Preallocate
  //capture process
  int M = z1.size();
  int J=X.nrow();
  LogicalVector inbox(1);
  NumericVector rand(1);
  double explldiff;
  NumericVector dtmp(J);
  NumericVector lamdLcand(J);
  NumericVector lamdRcand(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double v;
  double vcand;
  NumericMatrix s1new=Rcpp::clone(s1);
  NumericMatrix D1new=Rcpp::clone(D1);
  NumericMatrix s2new=Rcpp::clone(s2);
  NumericMatrix D2new=Rcpp::clone(D2);
  NumericMatrix lamdLnew=Rcpp::clone(lamdL);
  NumericMatrix lamdRnew=Rcpp::clone(lamdR);


  //Update Left
  for(int i=0; i<M; i++) {
    if(z1[i]==1){
      ScandX=Rcpp::rnorm(1,s1(i,0),.2);
      ScandY=Rcpp::rnorm(1,s1(i,1),.2);
      inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      if(inbox(0)){
        dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
        lamdLcand=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
        //Calculate likelihood for original and candidate
        v=fulllikInd1D2( z1(i), z2(i), lamdL(i,_),lamdR(i,_),yleft(i,_), yright(i,_),  X, K);
        vcand=fulllikInd1D2( z1(i), z2(i), lamdLcand,lamdR(i,_),yleft(i,_), yright(i,_),  X, K);
        rand=Rcpp::runif(1);
        explldiff=exp(vcand-v);
        if((rand(0)<explldiff)){
          s1new(i,0)=ScandX(0);
          s1new(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D1new(i,j) = dtmp(j);
            lamdLnew(i,j) = lamdLcand(j);
          }
        }
      }
    }
  }
  //Update Right
  for(int i=0; i<M; i++) {
    if(z2[i]==1){
      ScandX=Rcpp::rnorm(1,s2(i,0),.2);
      ScandY=Rcpp::rnorm(1,s2(i,1),.2);
      inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      if(inbox(0)){
        dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
        lamdRcand=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
        //Calculate likelihood for original and candidate
        v=fulllikInd1D2( z1(i), z2(i), lamdL(i,_), lamdR(i,_),yleft(i,_), yright(i,_),  X, K);
        vcand=fulllikInd1D2( z1(i), z2(i),lamdL(i,_),lamdRcand,yleft(i,_), yright(i,_),  X, K);
        rand=Rcpp::runif(1);
        explldiff=exp(vcand-v);
        if((rand(0)<explldiff)){
          s2new(i,0)=ScandX(0);
          s2new(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D2new(i,j) = dtmp(j);
            lamdRnew(i,j) = lamdRcand(j);
          }
        }
      }
    }
  }

  List to_return(8);
  to_return[0] = s1new;
  to_return[1] = s2new;
  to_return[2] = D1new;
  to_return[3] = D2new;
  to_return[4] = lamdLnew;
  to_return[5] = lamdRnew;
  return to_return;
}


//WholeShebang
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMCInd2(double lam01, double sigma,
              NumericMatrix yleft, NumericMatrix yright, IntegerVector z1,IntegerVector z2,NumericMatrix X,int K,
              NumericMatrix D1,NumericMatrix D2, IntegerVector capvector1, IntegerVector capvector2,
              NumericMatrix s1,NumericMatrix s2,NumericVector psi1, NumericVector psi2, NumericVector xlim,
              NumericVector ylim,LogicalVector useverts, NumericMatrix vertices,double proplam01,double propsigma, double propsx,
              double propsy,int niter, int nburn, int nthin) {
  RNGScope scope;
  int M = z1.size();
  int J=X.nrow();
  //Preallocate for update lam01  sigma
  //ll structures
  double likcurr=0;
  double liknew=0;
  double lam01cand=0;
  double sigmacand=0;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamdLcand;
  NumericMatrix lamdRcand;
  double explldiff;

  //Preallocate for Psi update
  NumericMatrix pdL(M,J);
  NumericMatrix pdLb(M,J);
  NumericMatrix pdLtraps(M,J);
  NumericMatrix pbar0L(M,J);
  NumericVector prob0L(M);
  NumericVector fc1(M);
  NumericMatrix pdR(M,J);
  NumericMatrix pdRb(M,J);
  NumericMatrix pdRtraps(M,J);
  NumericMatrix pbar0R(M,J);
  NumericVector prob0R(M);
  NumericVector fc2(M);
  int N1;
  int N2;

  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector dtmp(J);
  NumericVector lamdLcand1d(J);
  NumericVector lamdRcand1d(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double vcand;
  double v;

  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,4);
  NumericMatrix s1xout(nstore,M);
  NumericMatrix s1yout(nstore,M);
  NumericMatrix z1out(nstore,M);
  NumericMatrix s2xout(nstore,M);
  NumericMatrix s2yout(nstore,M);
  NumericMatrix z2out(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamdL=calclamd(lam01,sigma,D1);
  NumericMatrix lamdR=calclamd(lam01,sigma,D2);
  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    ////////update lam01 and sigma
    likcurr=fulllikInd2(z1,z2,lamdL,lamdR, yleft,yright,X,K);
    //Update lam01
    rand=Rcpp::rnorm(1,lam01,proplam01);
    if(rand(0) > 0){
      lam01cand=rand(0);
      lamdLcand=calclamd(lam01cand,sigma,D1);
      lamdRcand=calclamd(lam01cand,sigma,D2);
      liknew=fulllikInd2(z1,z2,lamdLcand,lamdRcand, yleft,yright,X,K);
      rand2=Rcpp::runif(1);
      explldiff = exp(liknew-likcurr);
      if(rand2(0)<explldiff){
        lam01=lam01cand;
        lamdL=lamdLcand;
        lamdR=lamdRcand;
        likcurr=liknew;
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamdLcand=calclamd(lam01,sigmacand,D1);
      lamdRcand=calclamd(lam01,sigmacand,D2);
      liknew=fulllikInd2(z1,z2,lamdLcand,lamdRcand,yleft,yright,X,K);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamdL=lamdLcand;
        lamdR=lamdRcand;
        likcurr=liknew;
      }
    }
    ///////update psi, z
    //  Calculate probability of no capture and update z1
    for(int i=0; i<M; i++) {
      prob0L(i)=1;
      for(int j=0; j<J; j++){
        pdL(i,j)=1-exp(-lamdL(i,j));
        if(X(j,2)==2){ //if double trap
          pdLb(i,j)=2*pdL(i,j)-pdL(i,j)*pdL(i,j);
          pdLtraps(i,j)=pdLb(i,j);
        }else{ //if single trap
          pdLtraps(i,j)=pdL(i,j);
        }
        pbar0L(i,j)=pow(1-pdLtraps(i,j),K);
        prob0L(i)*=pbar0L(i,j);
      }
      fc1(i)=prob0L(i)*psi1(0)/(prob0L(i)*psi1(0) + 1-psi1(0));
      if(capvector1(i)==0){
        NumericVector rand=Rcpp::rbinom(1,1,fc1(i));
        z1(i)=rand(0);
      }
    }
    //  Calculate probability of no capture and update z2
    for(int i=0; i<M; i++) {
      prob0R(i)=1;
      for(int j=0; j<J; j++){
        pdR(i,j)=1-exp(-lamdR(i,j));
        if(X(j,2)==2){ //if double trap
          pdRb(i,j)=2*pdR(i,j)-pdR(i,j)*pdR(i,j);
          pdRtraps(i,j)=pdRb(i,j);
        }else{ //if single trap
          pdRtraps(i,j)=pdR(i,j);
        }
        pbar0R(i,j)=pow(1-pdRtraps(i,j),K);
        prob0R(i)*=pbar0R(i,j);
      }
      fc2(i)=prob0R(i)*psi2(0)/(prob0R(i)*psi2(0) + 1-psi2(0));
      if(capvector2(i)==0){
        NumericVector rand=Rcpp::rbinom(1,1,fc2(i));
        z2(i)=rand(0);
      }
    }
    //Calculate current likelihood (no constants)
    v=fulllikInd2(z1,z2,lamdL,lamdR,yleft,yright,X,K);
    //update psis
    psi1=Rcpp::rbeta(1, 1 + sum(z1), 1 + M - sum(z1));
    psi2=Rcpp::rbeta(1, 1 + sum(z2), 1 + M - sum(z2));
    N1=sum(z1);
    N2=sum(z2);

    //////////Update ACs
    //Update Left
    for(int i=0; i<M; i++) {
      if(z1[i]==1){
        ScandX=Rcpp::rnorm(1,s1(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s1(i,1),propsy);
        if(useverts(0)==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
          dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
          lamdLcand1d=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
          //Calculate likelihood for original and candidate
          v=fulllikInd1D2( z1(i), z2(i), lamdL(i,_),lamdR(i,_),yleft(i,_), yright(i,_),  X, K);
          vcand=fulllikInd1D2( z1(i), z2(i), lamdLcand1d,lamdR(i,_),yleft(i,_), yright(i,_),  X, K);
          rand=Rcpp::runif(1);
          explldiff=exp(vcand-v);
          if((rand(0)<explldiff)){
            s1(i,0)=ScandX(0);
            s1(i,1)=ScandY(0);
            for(int j=0; j<J; j++){
              D1(i,j) = dtmp(j);
              lamdL(i,j) = lamdLcand1d(j);
            }
          }
        }
      }
    }
    //Update Right
    for(int i=0; i<M; i++) {
      if(z2[i]==1){
        ScandX=Rcpp::rnorm(1,s2(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s2(i,1),propsy);
        if(useverts(0)==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
          dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
          lamdRcand1d=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
          //Calculate likelihood for original and candidate
          v=fulllikInd1D2( z1(i), z2(i), lamdL(i,_), lamdR(i,_),yleft(i,_), yright(i,_),  X, K);
          vcand=fulllikInd1D2( z1(i), z2(i),lamdL(i,_),lamdRcand1d,yleft(i,_), yright(i,_),  X, K);
          rand=Rcpp::runif(1);
          explldiff=exp(vcand-v);
          if((rand(0)<explldiff)){
            s2(i,0)=ScandX(0);
            s2(i,1)=ScandY(0);
            for(int j=0; j<J; j++){
              D2(i,j) = dtmp(j);
              lamdR(i,j) = lamdRcand1d(j);
            }
          }
        }
      }
    }
    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      s1xout(iteridx,_)= s1(_,0);
      s1yout(iteridx,_)= s1(_,1);
      s2xout(iteridx,_)= s2(_,0);
      s2yout(iteridx,_)= s2(_,1);
      z1out(iteridx,_)= z1;
      z2out(iteridx,_)= z2;
      out(iteridx,0)=lam01;
      out(iteridx,1)=sigma;
      out(iteridx,2)=N1;
      out(iteridx,3)=N2;
      iteridx=iteridx+1;
    }
  }
  List to_return(7);
  to_return[0] = out;
  to_return[1] = s1xout;
  to_return[2] = s1yout;
  to_return[3] = s2xout;
  to_return[4] = s2yout;
  to_return[5] = z1out;
  to_return[6] = z2out;
  return to_return;
}

/////////////////////Heuristic estimator V3//////////////////////////

//likelihood calculation
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fulllikInd3(IntegerVector z1,IntegerVector z2,IntegerVector z3,NumericMatrix lamdB,NumericMatrix lamdL,
                   NumericMatrix lamdR,NumericMatrix yboth,NumericMatrix yleft,NumericMatrix yright, NumericMatrix X,int K) {
  //Preallocate
  int M = z1.size();
  int J=X.nrow();
  NumericMatrix pdB(M,J);
  NumericMatrix pdL(M,J);
  NumericMatrix pdLb(M,J);
  NumericMatrix pdR(M,J);
  NumericMatrix pdRb(M,J);
  NumericMatrix partboth(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);

  double v=0;
  //  Calculate likelihood
  for(int i=0; i<M; i++) {
    //Both data set first
    if(z1[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          pdB(i,j)=1-exp(-lamdB(i,j));
          partboth(i,j)=yboth(i,j)*log(pdB(i,j))+(K-yboth(i,j))*log(1-pdB(i,j));
        }else{ //if single trap
          partboth(i,j)=0;
        }
        if(partboth(i,j)==partboth(i,j)){
          v+=partboth(i,j);
        }
      }
    }
    //Left side data set next
    if(z2[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pdL(i,j)=1-exp(-lamdL(i,j));
        if(X(j,2)==2){ //if double trap
          pdLb(i,j)=2*pdL(i,j)-pdL(i,j)*pdL(i,j);
          partleft(i,j)=yleft(i,j)*log(pdLb(i,j))+(K-yleft(i,j))*log(1-pdLb(i,j));
        }else{ //if single trap
          partleft(i,j)=yleft(i,j)*log(pdL(i,j))+(K-yleft(i,j))*log(1-pdL(i,j));
        }
        if(partleft(i,j)==partleft(i,j)){
          v+=partleft(i,j);
        }
      }
    }
    //Now do right side
    if(z3[i]==1){ //if in pop
      for(int j=0; j<J; j++){
        pdR(i,j)=1-exp(-lamdR(i,j));
        if(X(j,2)==2){ //if double trap
          pdRb(i,j)=2*pdR(i,j)-pdR(i,j)*pdR(i,j);
          partright(i,j)=yright(i,j)*log(pdRb(i,j))+(K-yright(i,j))*log(1-pdRb(i,j));
        }else{ //if single trap
          partright(i,j)=yright(i,j)*log(pdR(i,j))+(K-yright(i,j))*log(1-pdR(i,j));
        }
        if(partright(i,j)==partright(i,j)){
          v+=partright(i,j);
        }
      }
    }
  }
  double to_return;
  to_return = v;
  return to_return;
}
//likelihood calculation3
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double fulllikInd1D3(int z1,int z2,int z3,NumericVector lamdB,NumericVector lamdL,NumericVector lamdR,
                     NumericVector yboth,NumericVector yleft,NumericVector yright, NumericMatrix X,int K) {
  //Preallocate
  int J=X.nrow();
  NumericVector pdB(J);
  NumericVector pdL(J);
  NumericVector pdLb(J);
  NumericVector pdR(J);
  NumericVector pdRb(J);
  NumericVector partboth(J);
  NumericVector partleft(J);
  NumericVector partright(J);
  double v=0;
  //  Calculate likelihood
  //both side data set first
  if(z1==1){ //if in pop
    for(int j=0; j<J; j++){
      if(X(j,2)==2){ //if double trap
        pdB(j)=1-exp(-lamdB(j));
        partboth(j)=yboth(j)*log(pdB(j))+(K-yboth(j))*log(1-pdB(j));
      }else{ //if single trap
        partboth(j)=0;
      }
      if(partboth(j)==partboth(j)){
        v+=partboth(j);
      }
    }
  }
  //Left side data set first
  if(z2==1){ //if in pop
    for(int j=0; j<J; j++){
      pdL(j)=1-exp(-lamdL(j));
      if(X(j,2)==2){ //if double trap
        pdLb(j)=2*pdL(j)-pdL(j)*pdL(j);
        partleft(j)=yleft(j)*log(pdLb(j))+(K-yleft(j))*log(1-pdLb(j));
      }else{ //if single trap
        partleft(j)=yleft(j)*log(pdL(j))+(K-yleft(j))*log(1-pdL(j));
      }
      if(partleft(j)==partleft(j)){
        v+=partleft(j);
      }
    }
  }
  //Now do right side data set
  if(z3==1){ //if in right side pop
    for(int j=0; j<J; j++){
      pdR(j)=1-exp(-lamdR(j));
      if(X(j,2)==2){ //if double trap
        pdRb(j)=2*pdR(j)-pdR(j)*pdR(j);
        partright(j)=yright(j)*log(pdRb(j))+(K-yright(j))*log(1-pdRb(j));
      }else{ //if single trap
        partright(j)=yright(j)*log(pdR(j))+(K-yright(j))*log(1-pdR(j));
      }
      if(partright(j)==partright(j)){
        v+=partright(j);
      }
    }
  }
  double to_return;
  to_return = v;
  return to_return;
}


//WholeShebang
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMCInd3(double lam01, double lam02, double sigma,NumericMatrix yboth,NumericMatrix yleft, NumericMatrix yright,
              IntegerVector z1,IntegerVector z2,IntegerVector z3,NumericMatrix X,int K,NumericMatrix D1,
              NumericMatrix D2,NumericMatrix D3, IntegerVector capvector1, IntegerVector capvector2,IntegerVector capvector3,
              NumericMatrix s1,NumericMatrix s2,NumericMatrix s3,NumericVector psi1,NumericVector psi2,NumericVector psi3,
              NumericVector xlim, NumericVector ylim,LogicalVector useverts, NumericMatrix vertices,double proplam01,
              double proplam02,double propsigma,double propsx,double propsy,int niter, int nburn, int nthin) {
  RNGScope scope;
  int M = z1.size();
  int J=X.nrow();
  //Preallocate for update lam01  sigma
  //ll structures
  double likcurr=0;
  double liknew=0;
  double lam01cand=0;
  double lam02cand=0;
  double sigmacand=0;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamdBcand;
  NumericMatrix lamdLcand;
  NumericMatrix lamdRcand;
  double explldiff;

  //Preallocate for Psi update
  NumericMatrix pdB(M,J);
  NumericMatrix pdBtraps(M,J);
  NumericMatrix pbar0B(M,J);
  NumericVector prob0B(M);
  NumericVector fc1(M);
  NumericMatrix pdL(M,J);
  NumericMatrix pdLb(M,J);
  NumericMatrix pdLtraps(M,J);
  NumericMatrix pbar0L(M,J);
  NumericVector prob0L(M);
  NumericVector fc2(M);
  NumericMatrix pdR(M,J);
  NumericMatrix pdRb(M,J);
  NumericMatrix pdRtraps(M,J);
  NumericMatrix pbar0R(M,J);
  NumericVector prob0R(M);
  NumericVector fc3(M);
  int N1;
  int N2;
  int N3;

  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector dtmp(J);
  NumericVector lamdBcand1d(J);
  NumericVector lamdLcand1d(J);
  NumericVector lamdRcand1d(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  double vcand;
  double v;

  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,6);
  NumericMatrix s1xout(nstore,M);
  NumericMatrix s1yout(nstore,M);
  NumericMatrix z1out(nstore,M);
  NumericMatrix s2xout(nstore,M);
  NumericMatrix s2yout(nstore,M);
  NumericMatrix z2out(nstore,M);
  NumericMatrix s3xout(nstore,M);
  NumericMatrix s3yout(nstore,M);
  NumericMatrix z3out(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamdB=calclamd(lam02,sigma,D1);
  NumericMatrix lamdL=calclamd(lam01,sigma,D2);
  NumericMatrix lamdR=calclamd(lam01,sigma,D3);
  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    ////////update lam01 and sigma
    likcurr=fulllikInd3(z1,z2,z3,lamdB,lamdL,lamdR,yboth,yleft,yright,X,K);
    //Update lam01
    rand=Rcpp::rnorm(1,lam01,proplam01);
    if(rand(0) > 0){
      lam01cand=rand(0);
      lamdLcand=calclamd(lam01cand,sigma,D2);
      lamdRcand=calclamd(lam01cand,sigma,D3);
      liknew=fulllikInd3(z1,z2,z3,lamdB,lamdLcand,lamdRcand,yboth, yleft,yright,X,K);
      rand2=Rcpp::runif(1);
      explldiff = exp(liknew-likcurr);
      if(rand2(0)<explldiff){
        lam01=lam01cand;
        lamdL=lamdLcand;
        lamdR=lamdRcand;
        likcurr=liknew;
      }
    }
    //Update lam02
    rand=Rcpp::rnorm(1,lam02,proplam02);
    if(rand(0) > 0){
      lam02cand=rand(0);
      lamdBcand=calclamd(lam02cand,sigma,D1);
      liknew=fulllikInd3(z1,z2,z3,lamdBcand,lamdL,lamdR,yboth, yleft,yright,X,K);
      rand2=Rcpp::runif(1);
      explldiff = exp(liknew-likcurr);
      if(rand2(0)<explldiff){
        lam02=lam02cand;
        lamdB=lamdBcand;
        likcurr=liknew;
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamdBcand=calclamd(lam02,sigmacand,D1);
      lamdLcand=calclamd(lam01,sigmacand,D2);
      lamdRcand=calclamd(lam01,sigmacand,D3);
      liknew=fulllikInd3(z1,z2,z3,lamdBcand,lamdLcand,lamdRcand,yboth,yleft,yright,X,K);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamdB=lamdBcand;
        lamdL=lamdLcand;
        lamdR=lamdRcand;
        likcurr=liknew;
      }
    }
    ///////update psi, z
    //  Calculate probability of no capture and update z1
    for(int i=0; i<M; i++) {
      prob0B(i)=1;
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          pdBtraps(i,j)=1-exp(-lamdB(i,j));;
        }else{ //if single trap
          pdBtraps(i,j)=0;
        }
        pbar0B(i,j)=pow(1-pdBtraps(i,j),K);
        prob0B(i)*=pbar0B(i,j);
      }
      fc1(i)=prob0B(i)*psi1(0)/(prob0B(i)*psi1(0) + 1-psi1(0));
      if(capvector1(i)==0){
        NumericVector rand=Rcpp::rbinom(1,1,fc1(i));
        z1(i)=rand(0);
      }
    }
    //  Calculate probability of no capture and update z2
    for(int i=0; i<M; i++) {
      prob0L(i)=1;
      for(int j=0; j<J; j++){
        pdL(i,j)=1-exp(-lamdL(i,j));
        if(X(j,2)==2){ //if double trap
          pdLb(i,j)=2*pdL(i,j)-pdL(i,j)*pdL(i,j);
          pdLtraps(i,j)=pdLb(i,j);
        }else{ //if single trap
          pdLtraps(i,j)=pdL(i,j);
        }
        pbar0L(i,j)=pow(1-pdLtraps(i,j),K);
        prob0L(i)*=pbar0L(i,j);
      }
      fc2(i)=prob0L(i)*psi2(0)/(prob0L(i)*psi2(0) + 1-psi2(0));
      if(capvector2(i)==0){
        NumericVector rand=Rcpp::rbinom(1,1,fc2(i));
        z2(i)=rand(0);
      }
    }
    //  Calculate probability of no capture and update z3
    for(int i=0; i<M; i++) {
      prob0R(i)=1;
      for(int j=0; j<J; j++){
        pdR(i,j)=1-exp(-lamdR(i,j));
        if(X(j,2)==2){ //if double trap
          pdRb(i,j)=2*pdR(i,j)-pdR(i,j)*pdR(i,j);
          pdRtraps(i,j)=pdRb(i,j);
        }else{ //if single trap
          pdRtraps(i,j)=pdR(i,j);
        }
        pbar0R(i,j)=pow(1-pdRtraps(i,j),K);
        prob0R(i)*=pbar0R(i,j);
      }
      fc3(i)=prob0R(i)*psi3(0)/(prob0R(i)*psi3(0) + 1-psi3(0));
      if(capvector3(i)==0){
        NumericVector rand=Rcpp::rbinom(1,1,fc3(i));
        z3(i)=rand(0);
      }
    }
    //Calculate current likelihood (no constants)
    v=fulllikInd3(z1,z2,z3,lamdB,lamdL,lamdR,yboth,yleft,yright,X,K);
    //update psis
    psi1=Rcpp::rbeta(1, 1 + sum(z1), 1 + M - sum(z1));
    psi2=Rcpp::rbeta(1, 1 + sum(z2), 1 + M - sum(z2));
    psi3=Rcpp::rbeta(1, 1 + sum(z3), 1 + M - sum(z3));
    N1=sum(z1);
    N2=sum(z2);
    N3=sum(z3);
    //////////Update ACs
    //Update Both
    for(int i=0; i<M; i++) {
      if(z1[i]==1){
        ScandX=Rcpp::rnorm(1,s1(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s1(i,1),propsy);
        if(useverts(0)==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
          dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
          lamdBcand1d=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
          //Calculate likelihood for original and candidate
          v=fulllikInd1D3( z1(i), z2(i),z3(i), lamdB(i,_), lamdL(i,_),lamdR(i,_),yboth(i,_),yleft(i,_), yright(i,_),  X, K);
          vcand=fulllikInd1D3( z1(i), z2(i),z3(i), lamdBcand1d,lamdL(i,_),lamdR(i,_),yboth(i,_),yleft(i,_), yright(i,_),  X, K);
          rand=Rcpp::runif(1);
          explldiff=exp(vcand-v);
          if((rand(0)<explldiff)){
            s1(i,0)=ScandX(0);
            s1(i,1)=ScandY(0);
            for(int j=0; j<J; j++){
              D1(i,j) = dtmp(j);
              lamdB(i,j) = lamdBcand1d(j);
            }
          }
        }
      }
    }
    //Update Left
    for(int i=0; i<M; i++) {
      if(z2[i]==1){
        ScandX=Rcpp::rnorm(1,s2(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s2(i,1),propsy);
        if(useverts(0)==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
          dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
          lamdLcand1d=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
          //Calculate likelihood for original and candidate
          v=fulllikInd1D3( z1(i), z2(i),z3(i),lamdB(i,_), lamdL(i,_),lamdR(i,_),yboth(i,_),yleft(i,_), yright(i,_),  X, K);
          vcand=fulllikInd1D3( z1(i), z2(i),z3(i),lamdB(i,_), lamdLcand1d,lamdR(i,_),yboth(i,_),yleft(i,_), yright(i,_),  X, K);
          rand=Rcpp::runif(1);
          explldiff=exp(vcand-v);
          if((rand(0)<explldiff)){
            s2(i,0)=ScandX(0);
            s2(i,1)=ScandY(0);
            for(int j=0; j<J; j++){
              D2(i,j) = dtmp(j);
              lamdL(i,j) = lamdLcand1d(j);
            }
          }
        }
      }
    }
    //Update Right
    for(int i=0; i<M; i++) {
      if(z3[i]==1){
        ScandX=Rcpp::rnorm(1,s3(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s3(i,1),propsy);
        if(useverts(0)==FALSE){
          inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
        }else{
          inbox=inoutCpp(ScandX,ScandY,vertices);
        }
        if(inbox(0)){
          dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
          lamdRcand1d=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
          //Calculate likelihood for original and candidate
          v=fulllikInd1D3( z1(i), z2(i),z3(i), lamdB(i,_),lamdL(i,_), lamdR(i,_),yboth(i,_),yleft(i,_), yright(i,_),  X, K);
          vcand=fulllikInd1D3( z1(i), z2(i),z3(i),lamdB(i,_),lamdL(i,_),lamdRcand1d,yboth(i,_),yleft(i,_), yright(i,_),  X, K);
          rand=Rcpp::runif(1);
          explldiff=exp(vcand-v);
          if((rand(0)<explldiff)){
            s3(i,0)=ScandX(0);
            s3(i,1)=ScandY(0);
            for(int j=0; j<J; j++){
              D3(i,j) = dtmp(j);
              lamdR(i,j) = lamdRcand1d(j);
            }
          }
        }
      }
    }
    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      s1xout(iteridx,_)= s1(_,0);
      s1yout(iteridx,_)= s1(_,1);
      s2xout(iteridx,_)= s2(_,0);
      s2yout(iteridx,_)= s2(_,1);
      s3xout(iteridx,_)= s3(_,0);
      s3yout(iteridx,_)= s3(_,1);
      z1out(iteridx,_)= z1;
      z2out(iteridx,_)= z2;
      z3out(iteridx,_)= z3;
      out(iteridx,0)=lam01;
      out(iteridx,1)=lam02;
      out(iteridx,2)=sigma;
      out(iteridx,3)=N1;
      out(iteridx,4)=N2;
      out(iteridx,5)=N3;
      iteridx=iteridx+1;
    }
  }
  List to_return(10);
  to_return[0] = out;
  to_return[1] = s1xout;
  to_return[2] = s1yout;
  to_return[3] = s2xout;
  to_return[4] = s2yout;
  to_return[5] = s3xout;
  to_return[6] = s3yout;
  to_return[7] = z1out;
  to_return[8] = z2out;
  to_return[9] = z3out;
  return to_return;
}

/////////////////////////////Inhomogenous Density SCR, discrete state space//////////////////////////
//update betas
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List upbeta(double beta0, double beta1, NumericMatrix grid, IntegerVector scell, double propbeta0,double propbeta1,double EN, int N,double cellArea, IntegerVector z) {
  RNGScope scope;
  int M = z.size();
  double ENcand;
  double betacand;
  double llbeta;
  double llbetacand;
  NumericVector rand;
  NumericVector rand2;
  int npix=grid.nrow();
  //Update beta0
  rand=Rcpp::rnorm(1,beta0,propbeta0);
  ENcand=0;
  betacand=rand(0);
  for(int i=0; i<npix; i++) {
    if(grid(i,2)==1){
      ENcand+=exp(betacand)*cellArea;
    }else{
      ENcand+=exp(betacand+beta1)*cellArea;
    }
  }
  llbeta=0;
  for(int i=0; i<M; i++) {
    if(grid(scell(i)-1,2)==1){
      llbeta+=(beta0-log(EN))*z(i);
    }else{
      llbeta+=(beta0+beta1-log(EN))*z(i);
    }
  }
  llbeta=llbeta+N*log(EN/M)+(M-N)*log(1-EN/M);
  if(ENcand < M){
    llbetacand=0;
    for(int i=0; i<M; i++) {
      if(grid(scell(i)-1,2)==1){
        llbetacand+=(betacand-log(ENcand))*z(i);
      }else{
        llbetacand+=(betacand+beta1-log(ENcand))*z(i);
      }
    }
    llbetacand=llbetacand+N*log(ENcand/M)+(M-N)*log(1-ENcand/M);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(llbetacand-llbeta)){
      beta0=betacand;
      EN=ENcand;
      llbeta=llbetacand;
    }
  }
  //Update beta1
  rand=Rcpp::rnorm(1,beta1,propbeta1);
  ENcand=0;
  betacand=rand(0);
  for(int i=0; i<npix; i++) {
    if(grid(i,2)==1){
      ENcand+=exp(beta0)*cellArea;
    }else{
      ENcand+=exp(beta0+betacand)*cellArea;
    }
  }
  if(ENcand < M){
    llbetacand=0;
    for(int i=0; i<M; i++) {
      if(grid(scell(i)-1,2)==1){
        llbetacand+=(beta0-log(ENcand))*z(i);
      }else{
        llbetacand+=(beta0+betacand-log(ENcand))*z(i);
      }
    }
    llbetacand=llbetacand+N*log(ENcand/M)+(M-N)*log(1-ENcand/M);
    rand2=Rcpp::runif(1);
    if(rand2(0)<exp(llbetacand-llbeta)){
      beta1=betacand;
      EN=ENcand;
      llbeta=llbetacand;
    }
  }
  List to_return(3);
  to_return[0] = beta0;
  to_return[1] = beta1;
  to_return[2] = ENcand;
  return to_return;
}
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMC1b(double lam0, double sigma,double beta0, double beta1, NumericMatrix y, IntegerVector z, NumericMatrix X,int K,NumericMatrix D, int N,
            IntegerVector knownvector, NumericMatrix s,IntegerVector scell,NumericVector psi, NumericMatrix grid,double cellArea, double EN,
            double proplam0, double propsigma,double propbeta0,double propbeta1,double propsx,double propsy, int niter, int nburn, int nthin) {
  RNGScope scope;
  int M = z.size();
  //Preallocate for update lam0  sigma, betas
  double likcurr;
  double liknew;
  double lam0cand;
  double sigmacand;
  double ENcand;
  double betacand;
  double llbeta;
  double llbetacand;
  double explldiff;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamdcand;
  int J=y.ncol();
  int npix=grid.nrow();

  //Preallocate for z update
  NumericMatrix pd(M,J);
  NumericMatrix pbar(M,J);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericMatrix v1(M,J);

  //Preallocate for updating activity centers
  NumericVector pdc(J);
  NumericVector pdccand(J);
  NumericVector vc(J);
  NumericVector vcCand(J);
  NumericVector dtmp(J);
  NumericVector dtmp2(npix);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  NumericVector Scandcell(1);
  double mindist;
  double vcand;
  double priors;
  double priorscand;
  NumericVector lamdcandc(J);

  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,6);
  NumericMatrix sxout(nstore,M);
  NumericMatrix syout(nstore,M);
  NumericMatrix zout(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamd=calclamd(lam0,sigma,D);

  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    likcurr=SCRlik( z, lamd, y, K);
    liknew=0;
    lam0cand=0;
    sigmacand=0;
    //Update lam0
    rand=Rcpp::rnorm(1,lam0,proplam0);
    if(rand(0) > 0){
      lam0cand=rand(0);
      lamdcand=calclamd(lam0cand,sigma,D);
      liknew=SCRlik( z, lamdcand, y, K);
      rand2=Rcpp::runif(1);
      explldiff = exp(liknew-likcurr);
      if(rand2(0)<explldiff){
        lam0=lam0cand;
        lamd=lamdcand;
        likcurr=liknew;
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamdcand=calclamd(lam0,sigmacand,D);
      liknew=SCRlik( z, lamdcand, y, K);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamd=lamdcand;
        likcurr=liknew;
      }
    }
    //Update z and N
    //  Calculate probability of no capture and update z
    N=0;
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        pd(i,j)=1-exp(-lamd(i,j));
        pbar(i,j)=pow(1-pd(i,j),K);
        prob0(i)*=pbar(i,j);
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0);
      if(swappable(i)){
        NumericVector rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
      N=N+z(i);
    }
    //Calculate current likelihood
    double v=0;
    for(int i=0; i<M; i++) {
      if(z[i]==1){
        for(int j=0; j<J; j++){
          v1(i,j)=y(i,j)*log(pd(i,j))+(K-y(i,j))*log(1-pd(i,j));
          if(v1(i,j)==v1(i,j)){
            v+=v1(i,j);
          }
        }
      }
    }
    //Update beta0
    rand=Rcpp::rnorm(1,beta0,propbeta0);
    ENcand=0;
    betacand=rand(0);
    for(int i=0; i<npix; i++) {
      if(grid(i,2)==1){
        ENcand+=exp(betacand)*cellArea;
      }else{
        ENcand+=exp(betacand+beta1)*cellArea;
      }
    }
    llbeta=0;
    for(int i=0; i<M; i++) {
      if(grid(scell(i)-1,2)==1){
        llbeta+=(beta0-log(EN))*z(i);
      }else{
        llbeta+=(beta0+beta1-log(EN))*z(i);
      }
    }
    llbeta=llbeta+N*log(EN/M)+(M-N)*log(1-EN/M);
    if(ENcand < M){
      llbetacand=0;
      for(int i=0; i<M; i++) {
        if(grid(scell(i)-1,2)==1){
          llbetacand+=(betacand-log(ENcand))*z(i);
        }else{
          llbetacand+=(betacand+beta1-log(ENcand))*z(i);
        }
      }
      llbetacand=llbetacand+N*log(ENcand/M)+(M-N)*log(1-ENcand/M);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(llbetacand-llbeta)){
        beta0=betacand;
        EN=ENcand;
        llbeta=llbetacand;
      }
    }
    //Update beta1
    rand=Rcpp::rnorm(1,beta1,propbeta1);
    ENcand=0;
    betacand=rand(0);
    for(int i=0; i<npix; i++) {
      if(grid(i,2)==1){
        ENcand+=exp(beta0)*cellArea;
      }else{
        ENcand+=exp(beta0+betacand)*cellArea;
      }
    }
    if(ENcand < M){
      llbetacand=0;
      for(int i=0; i<M; i++) {
        if(grid(scell(i)-1,2)==1){
          llbetacand+=(beta0-log(ENcand))*z(i);
        }else{
          llbetacand+=(beta0+betacand-log(ENcand))*z(i);
        }
      }
      llbetacand=llbetacand+N*log(ENcand/M)+(M-N)*log(1-ENcand/M);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(llbetacand-llbeta)){
        beta1=betacand;
        EN=ENcand;
        llbeta=llbetacand;
      }
    }
    psi=EN/M;

    //Update Activity Centers
    for(int i=0; i<M; i++) {
      if(z(i)==1){ //if in pop
        ScandX=Rcpp::rnorm(1,s(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s(i,1),propsy);
        //snap location to nearest grid point
        mindist=10000000;
        Scandcell=0;
        for(int j=0; j<npix; j++){
          dtmp2(j)=pow( pow(ScandX(0) - grid(j,0), 2.0) + pow( ScandY(0) - grid(j,1), 2.0), 0.5 );
          if(dtmp2(j)<mindist){
            mindist=dtmp2(j);
            Scandcell=j+1;
          }
        }
        //update snapped location
        ScandX(0)=grid(Scandcell(0)-1,0);
        ScandY(0)=grid(Scandcell(0)-1,1);
        dtmp=pow( pow( rep(ScandX(0),J) - X(_,0), 2.0) + pow( rep(ScandY(0),J) - X(_,1), 2.0), 0.5 );
        lamdcandc=lam0*exp(-dtmp*dtmp/(2*sigma*sigma));
        v=0;
        vcand=0;
        //Calculate likelihood for original and candidate
        for(int j=0; j<J; j++){
          pdc(j)=1-exp(-lamd(i,j));
          pdccand(j)=1-exp(-lamdcandc(j));
          vc(j)=y(i,j)*log(pdc(j))+(K-y(i,j))*log(1-pdc(j));
          vcCand(j)=y(i,j)*log(pdccand(j))+(K-y(i,j))*log(1-pdccand(j));
          if(vc(j)==vc(j)){
            v+=vc(j);
          }
          if(vcCand(j)==vcCand(j)){
            vcand+=vcCand(j);
          }
        }
        //Calculate prior for original and candidate
        if(grid(scell(i)-1,2)==2){
          priors=beta0 + beta1;
        }else{
          priors=beta0;
        }
        if(grid(Scandcell(0)-1,2)==2){
          priorscand=beta0 + beta1;
        }else{
          priorscand=beta0;
        }
        rand=Rcpp::runif(1);
        explldiff=exp((vcand+priorscand)-(v+priors));
        if((rand(0)<explldiff)){
          scell(i)=Scandcell(0);
          s(i,0)=ScandX(0);
          s(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D(i,j) = dtmp(j);
            lamd(i,j) = lamdcandc(j);
          }
        }
      }
    }
    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      sxout(iteridx,_)= s(_,0);
      syout(iteridx,_)= s(_,1);
      zout(iteridx,_)= z;
      out(iteridx,0)=lam0;
      out(iteridx,1)=sigma;
      out(iteridx,2)=N;
      out(iteridx,3)=EN;
      out(iteridx,4)=beta0;
      out(iteridx,5)=beta1;
      iteridx=iteridx+1;
    }
  }
  List to_return(4);
  to_return[0] = out;
  to_return[1] = sxout;
  to_return[2] = syout;
  to_return[3] = zout;
  return to_return;
}

////////////////////////2side model with inhomogenous density////////////////////////////////
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMC2(double lam01,double lam02, double sigma,double beta0, double beta1,
          NumericMatrix yboth, NumericMatrix yleft, NumericMatrix yright, IntegerVector z, NumericMatrix X,int K,
          NumericMatrix D,int Nfixed, IntegerVector knownvector,IntegerVector ID_L, IntegerVector ID_R,int swap,double swaptol,
          IntegerVector left,IntegerVector right,NumericMatrix s,IntegerVector scell,NumericVector psi,  NumericMatrix grid,double cellArea,
          double EN, double proplam01, double proplam02,double propsigma, double propbeta0, double propbeta1,double propsx, double propsy,
          int niter, int nburn, int nthin, LogicalVector updates) {
  RNGScope scope;
  int M = z.size();
  int npix=grid.nrow();
  // Preallocate for update lam01 lam02 sigma
  double likcurr;
  double liknew;
  double lam01cand;
  double lam02cand;
  double sigmacand;
  double ENcand;
  double betacand;
  double llbeta;
  double llbetacand;
  double explldiff;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd1cand;
  NumericMatrix lamd2cand;
  //Preallocate side swapping structures
  IntegerVector IDs=seq_len(M);
  int guy1=0;
  int guy2=0;
  int swapin;
  int swapout;
  int idx;
  double ncand;
  double ncand2;
  int match;
  double jumpprob;
  double backprob;
  IntegerVector possible;
  IntegerVector unpossible;
  IntegerVector newID(M);
  int J=yleft.ncol();
  NumericMatrix ylefttmp(M,J);
  NumericMatrix yrighttmp(M,J);  //could use same structure but dont want to get confused

  //ll structures
  NumericMatrix prepart(2,J);
  NumericMatrix postpart(2,J);
  NumericMatrix pd1b(2,J);
  NumericMatrix pd12b(2,J);
  NumericVector swapped(2);
  double llpre;
  double llpost;

  //MH structures
  double lldiff=0;

  //Housekeeping structures
  NumericVector guycounts(M);
  LogicalVector zeroguys(M);

  //Preallocate for z, N update
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix pd1traps(M,J);
  NumericMatrix pd2traps(M,J);
  NumericMatrix pbar1(M,J);
  NumericMatrix pbar2(M,J);
  NumericMatrix pbar0(M,J);
  NumericVector prob0(M);
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericMatrix partboth(M,J);
  NumericMatrix partleft(M,J);
  NumericMatrix partright(M,J);
  int N;

  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector pd1c(J);
  NumericVector pd12c(J);
  NumericVector pd2c(J);
  NumericVector pd1ccand(J);
  NumericVector pd12ccand(J);
  NumericVector pd2ccand(J);
  NumericVector partbothc(J);
  NumericVector partleftc(J);
  NumericVector partrightc(J);
  NumericVector partbothcCand(J);
  NumericVector partleftcCand(J);
  NumericVector partrightcCand(J);
  NumericVector dtmp(J);
  NumericVector dtmp2(npix);
  NumericVector lamd1candc(J);
  NumericVector lamd2candc(J);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  NumericVector Scandcell(1);

  double mindist;
  double vcand;
  double priors;
  double priorscand;
  double v;
  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,7);
  NumericMatrix sxout(nstore,M);
  NumericMatrix syout(nstore,M);
  NumericMatrix zout(nstore,M);
  NumericMatrix ID_Lout(nstore,M);
  NumericMatrix ID_Rout(nstore,M);
  int iteridx=0;
  //calc lamds  D too? Currently an input
  NumericMatrix lamd1=calclamd(lam01,sigma,D);
  NumericMatrix lamd2=calclamd(lam02,sigma,D);
  //start for loop here
  int iter;
  for(iter=0; iter<niter; iter++){
    //update lam01, lam02, and sigma

    likcurr=fulllik( z, lamd1, lamd2, yboth,  yleft, yright,  X, K);
    liknew=0;
    lam01cand=0;
    lam02cand=0;
    sigmacand=0;

    //Update lam01
    if(updates(0)){
      rand=Rcpp::rnorm(1,lam01,proplam01);
      if(rand(0) > 0){
        lam01cand=rand(0);
        lamd1cand=calclamd(lam01cand,sigma,D);
        liknew=fulllik( z, lamd1cand, lamd2, yboth,  yleft, yright,  X, K);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam01=lam01cand;
          lamd1=lamd1cand;
          likcurr=liknew;
        }
      }
    }
    //Update lam02
    if(updates(1)){
      rand=Rcpp::rnorm(1,lam02,proplam02);
      if(rand(0) > 0){
        lam02cand=rand(0);
        lamd2cand=calclamd(lam02cand,sigma,D);
        liknew=fulllik( z, lamd1, lamd2cand, yboth,  yleft, yright,  X, K);
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(liknew-likcurr)){
          lam02=lam02cand;
          lamd2=lamd2cand;
          likcurr=liknew;
        }
      }
    }
    //Update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lamd1cand=calclamd(lam01,sigmacand,D);
      lamd2cand=calclamd(lam02,sigmacand,D);
      liknew=fulllik( z, lamd1cand, lamd2cand, yboth,  yleft, yright,  X, K);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(liknew-likcurr)){
        sigma=sigmacand;
        lamd1=lamd1cand;
        lamd2=lamd2cand;
        likcurr=liknew;
      }
    }

    //Update z, N
    //  Calculate probability of no capture and update z and N
    N=0;
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        pd1(i,j)=1-exp(-lamd1(i,j));
        if(X(j,2)==2){ //if double trap
          pd2(i,j)=1-exp(-lamd2(i,j));
          pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
          pd1traps(i,j)=pd12(i,j);
          pd2traps(i,j)=pd2(i,j);
        }else{ //if single trap
          pd1traps(i,j)=pd1(i,j);
          pd2traps(i,j)=0;
        }
        pbar1(i,j)=pow(1-pd1traps(i,j),2*K);
        pbar2(i,j)=pow(1-pd2traps(i,j),K);
        pbar0(i,j)=pbar1(i,j)*pbar2(i,j);
        prob0(i)*=pbar0(i,j);
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0)&zeroguys(i);
      if(swappable(i)){
        rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
      N=N+z(i);
    }
    //Calculate current likelihood
    v=0;
    for(int i=0; i<M; i++) {
      if(z[i]==1){
        for(int j=0; j<J; j++){
          if(X(j,2)==2){ //if double trap
            partboth(i,j)=yboth(i,j)*log(pd2(i,j))+(K-yboth(i,j))*log(1-pd2(i,j));
            partleft(i,j)=yleft(i,j)*log(pd12(i,j))+(K-yleft(i,j))*log(1-pd12(i,j));
            partright(i,j)=yright(i,j)*log(pd12(i,j))+(K-yright(i,j))*log(1-pd12(i,j));
          }else{ //if single trap
            partboth(i,j)=0;
            partleft(i,j)=yleft(i,j)*log(pd1(i,j))+(K-yleft(i,j))*log(1-pd1(i,j));
            partright(i,j)=yright(i,j)*log(pd1(i,j))+(K-yright(i,j))*log(1-pd1(i,j));
          }
          if(partboth(i,j)==partboth(i,j)){
            v+=partboth(i,j);
          }
          if(partleft(i,j)==partleft(i,j)){
            v+=partleft(i,j);
          }
          if(partright(i,j)==partright(i,j)){
            v+=partright(i,j);
          }
        }
      }
    }
    likcurr=v;

    //Update beta0
    rand=Rcpp::rnorm(1,beta0,propbeta0);
    ENcand=0;
    betacand=rand(0);
    for(int i=0; i<npix; i++) {
      if(grid(i,2)==1){
        ENcand+=exp(betacand)*cellArea;
      }else{
        ENcand+=exp(betacand+beta1)*cellArea;
      }
    }
    llbeta=0;
    for(int i=0; i<M; i++) {
      if(grid(scell(i)-1,2)==1){
        llbeta+=(beta0-log(EN))*z(i);
      }else{
        llbeta+=(beta0+beta1-log(EN))*z(i);
      }
    }
    llbeta=llbeta+N*log(EN/M)+(M-N)*log(1-EN/M);
    if(ENcand < M){
      llbetacand=0;
      for(int i=0; i<M; i++) {
        if(grid(scell(i)-1,2)==1){
          llbetacand+=(betacand-log(ENcand))*z(i);
        }else{
          llbetacand+=(betacand+beta1-log(ENcand))*z(i);
        }
      }
      llbetacand=llbetacand+N*log(ENcand/M)+(M-N)*log(1-ENcand/M);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(llbetacand-llbeta)){
        beta0=betacand;
        EN=ENcand;
        llbeta=llbetacand;
      }
    }
    //Update beta1
    rand=Rcpp::rnorm(1,beta1,propbeta1);
    ENcand=0;
    betacand=rand(0);
    for(int i=0; i<npix; i++) {
      if(grid(i,2)==1){
        ENcand+=exp(beta0)*cellArea;
      }else{
        ENcand+=exp(beta0+betacand)*cellArea;
      }
    }
    if(ENcand < M){
      llbetacand=0;
      for(int i=0; i<M; i++) {
        if(grid(scell(i)-1,2)==1){
          llbetacand+=(beta0-log(ENcand))*z(i);
        }else{
          llbetacand+=(beta0+betacand-log(ENcand))*z(i);
        }
      }
      llbetacand=llbetacand+N*log(ENcand/M)+(M-N)*log(1-ENcand/M);
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp(llbetacand-llbeta)){
        beta1=betacand;
        EN=ENcand;
        llbeta=llbetacand;
      }
    }
    psi=EN/M;

    //  Update left sides
    if(updates(2)){
      //Build mapL
      IntegerMatrix mapL(M,2);
      mapL(_,1)=ID_L;
      mapL(_,0)=IDs;
      IntegerMatrix candmapL=Rcpp::clone(mapL);
      for (int i=0; i<M; i++) {
        if(z[candmapL(i,1)-1]==0){
          candmapL(i,1)= -1;
        }
        if(i<Nfixed){
          candmapL(i,0)= -1;
          candmapL(i,1)= -1;
        }
        if(z(i)==0){
          candmapL(i,0)= -1;
        }
      }

      //Find outcandsL and incandsL
      IntegerVector outcandsL=IDs[candmapL(_,1)>0];
      IntegerVector incandsL=IDs[candmapL(_,0)>0];
      int insizeL=incandsL.size();
      int outsizeL=outcandsL.size();

      //Structures to store distances
      NumericVector Dl(incandsL.size());
      NumericVector trashL(insizeL);

      //Build left swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsL(outsizeL+1);
      double x;
      for (int i=0; i<(outsizeL+1); i++) {
        binsL(i)=i;
        x=binsL(i)/outsizeL;
        binsL(i)=x;
      }

      ///////////////// //swap left sides//////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeL+1); i++) {
          if(binsL(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsL(idx-1);

        //find swapout
        swapout=mapL(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeL; i++) {
          Dl(i) = sqrt( pow( s(swapout-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        possible=incandsL[Dl < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsL2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsL2(i)=i;
            x=binsL2(i)/ncand;
            binsL2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsL.size(); i++) {
          trashL(i) = sqrt( pow( s(swapin-1,0) - s(incandsL(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsL(i)-1,1), 2.0) );
        }
        unpossible=incandsL[trashL < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapL(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_L
        newID=Rcpp::clone(ID_L);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new left data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            ylefttmp(newID(i)-1,j)=0;
            for (int k=0; k<K; k++) {
              ylefttmp(newID(i)-1,j)+=left(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for left data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              prepart(i,j)=yleft(swapped(i)-1,j)*log(pd12(i,j))+(K-yleft(swapped(i)-1,j))*log(1-pd12(i,j));
              postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd12(i,j))+(K-ylefttmp(swapped(i)-1,j))*log(1-pd12(i,j));
            }else{ //if single trap
              prepart(i,j)=yleft(swapped(i)-1,j)*log(pd1(i,j))+(K-yleft(swapped(i)-1,j))*log(1-pd1(i,j));
              postpart(i,j)=ylefttmp(swapped(i)-1,j)*log(pd1(i,j))+(K-ylefttmp(swapped(i)-1,j))*log(1-pd1(i,j));
            }
            if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
              llpre+=prepart(i,j);
            }
            if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
              llpost+=postpart(i,j);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yleft=Rcpp::clone(ylefttmp);
          ID_L=Rcpp::clone(newID);
          mapL(guy1-1,1)=swapin;
          mapL(guy2-1,1)=swapout;
        }

      }
    }
    //Update Right sides
    if(updates(3)){
      //Build mapR
      IntegerMatrix mapR(M,2);
      mapR(_,1)=ID_R;
      mapR(_,0)=IDs;
      IntegerMatrix candmapR=Rcpp::clone(mapR);
      for (int i=0; i<M; i++) {
        if(z[candmapR(i,1)-1]==0){
          candmapR(i,1)= -1;
        }
        if(i<Nfixed){
          candmapR(i,0)= -1;
          candmapR(i,1)= -1;
        }
        if(z(i)==0){
          candmapR(i,0)= -1;
        }
      }

      //Find outcandsR and incandsR
      IntegerVector outcandsR=IDs[candmapR(_,1)>0];
      IntegerVector incandsR=IDs[candmapR(_,0)>0];
      int insizeR=incandsR.size();
      int outsizeR=outcandsR.size();

      //Structures to store distances
      NumericVector Dr(incandsR.size());
      NumericVector trashR(insizeR);

      //Build right swapout bins for sample function
      //swapin bins vary in size so must be calculated in for loop
      NumericVector binsR(outsizeR+1);
      double x;
      for (int i=0; i<(outsizeR+1); i++) {
        binsR(i)=i;
        x=binsR(i)/outsizeR;
        binsR(i)=x;
      }

      ///////////////// //swap right sides/////////////////
      for (int l=0; l<swap; l++) {
        //Find guy 1
        rand=Rcpp::runif(1);
        idx=0;
        for (int i=0; i<(outsizeR+1); i++) {
          if(binsR(i)<rand(0)){
            idx=idx+1;
          }
        }
        guy1=outcandsR(idx-1);

        //find swapout
        swapout=mapR(guy1-1,1);

        //calculate distances and find possible switches
        for (int i=0; i<insizeR; i++) {
          Dr(i) = sqrt( pow( s(swapout-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapout-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        possible=incandsR[Dr < swaptol];
        ncand=possible.size();
        jumpprob=1/ncand;

        //swap in
        if(ncand>1){
          NumericVector binsR2(ncand+1);
          rand=Rcpp::runif(1);
          idx=0;
          for (int i=0; i<(ncand+1); i++) {
            binsR2(i)=i;
            x=binsR2(i)/ncand;
            binsR2(i)=x;
            if(x<rand(0)){
              idx=idx+1;
            }
          }
          swapin=possible(idx-1);
        }else{
          swapin=possible(0);
        }
        for (int i=0; i<incandsR.size(); i++) {
          trashR(i) = sqrt( pow( s(swapin-1,0) - s(incandsR(i)-1,0), 2.0) + pow( s(swapin-1,1) - s(incandsR(i)-1,1), 2.0) );
        }
        unpossible=incandsR[trashR < swaptol];
        ncand2=unpossible.size();
        backprob=1/ncand2;
        //find guy2
        for (int i=0; i<M; i++) {
          match=mapR(i,1);
          if(match==swapin){
            guy2=i+1;
          }
        }
        //update ID_R
        newID=Rcpp::clone(ID_R);
        newID(guy1-1)=swapin;
        newID(guy2-1)=swapout;

        //Make new right data set
        idx=0;
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            yrighttmp(newID(i)-1,j)=0;
            for (int k=0; k<K; k++) {
              yrighttmp(newID(i)-1,j)+=right(idx);
              idx+=1;
            }
          }
        }
        //Calculate likelihoods pre and post switch. Only need to calc for right side data
        swapped(0)=swapout;
        swapped(1)=swapin;
        llpre=0;
        llpost=0;
        for(int i=0; i<2; i++){
          for(int j=0; j<J; j++){
            pd1(i,j)=1-exp(-lamd1(swapped(i)-1,j));
            if(X(j,2)==2){ //if double trap
              pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
              prepart(i,j)=yright(swapped(i)-1,j)*log(pd12(i,j))+(K-yright(swapped(i)-1,j))*log(1-pd12(i,j));
              postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd12(i,j))+(K-yrighttmp(swapped(i)-1,j))*log(1-pd12(i,j));
            }else{ //if single trap
              prepart(i,j)=yright(swapped(i)-1,j)*log(pd1(i,j))+(K-yright(swapped(i)-1,j))*log(1-pd1(i,j));
              postpart(i,j)=yrighttmp(swapped(i)-1,j)*log(pd1(i,j))+(K-yrighttmp(swapped(i)-1,j))*log(1-pd1(i,j));
            }
            if((prepart(i,j)==prepart(i,j))&(!std::isinf(prepart(i,j)))){
              llpre+=prepart(i,j);
            }
            if((postpart(i,j)==postpart(i,j))&(!std::isinf(postpart(i,j)))){
              llpost+=postpart(i,j);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        lldiff=llpost-llpre;
        if(rand(0)<(exp(lldiff)*(backprob/jumpprob))){
          // likcurr+=lldiff;
          yright=Rcpp::clone(yrighttmp);
          ID_R=Rcpp::clone(newID);
          mapR(guy1-1,1)=swapin;
          mapR(guy2-1,1)=swapout;
        }

      } //end swapping rights
    }

    //Recalculate zeroguys
    for (int i=0; i<M; i++) {
      guycounts(i)=0;
      for(int j=0; j<J; j++) {
        guycounts(i)+=yboth(i,j)+yright(i,j)+yleft(i,j);
      }
      zeroguys(i)=guycounts(i)==0;
    }

    //Update Activity Centers
    for(int i=0; i<M; i++) {
      if(z[i]==1){ //if in pop
        ScandX=Rcpp::rnorm(1,s(i,0),propsx);
        ScandY=Rcpp::rnorm(1,s(i,1),propsy);
        //snap location to nearest grid point
        mindist=10000000;
        Scandcell=0;
        for(int j=0; j<npix; j++){
          dtmp2(j)=pow( pow(ScandX(0) - grid(j,0), 2.0) + pow( ScandY(0) - grid(j,1), 2.0), 0.5 );
          if(dtmp2(j)<mindist){
            mindist=dtmp2(j);
            Scandcell=j+1;
          }
        }
        //update snapped location
        ScandX(0)=grid(Scandcell(0)-1,0);
        ScandY(0)=grid(Scandcell(0)-1,1);
        //business as usual
        dtmp=pow( pow( rep(ScandX,J) - X(_,0), 2.0) + pow( rep(ScandY,J) - X(_,1), 2.0), 0.5 );
        lamd1candc=lam01*exp(-dtmp*dtmp/(2*sigma*sigma));
        lamd2candc=lam02*exp(-dtmp*dtmp/(2*sigma*sigma));
        v=0;
        vcand=0;
        //Calculate likelihood for original and candidate
        for(int j=0; j<J; j++){
          pd1c(j)=1-exp(-lamd1(i,j));
          pd1ccand(j)=1-exp(-lamd1candc(j));
          if(X(j,2)==2){ //if double trap
            pd12c(j)=2*pd1c(j)-pd1c(j)*pd1c(j);
            pd12ccand(j)=2*pd1ccand(j)-pd1ccand(j)*pd1ccand(j);
            pd2c(j)=1-exp(-lamd2(i,j));
            pd2ccand(j)=1-exp(-lamd2candc(j));
            partbothc(j)=yboth(i,j)*log(pd2c(j))+(K-yboth(i,j))*log(1-pd2c(j));
            partleftc(j)=yleft(i,j)*log(pd12c(j))+(K-yleft(i,j))*log(1-pd12c(j));
            partrightc(j)=yright(i,j)*log(pd12c(j))+(K-yright(i,j))*log(1-pd12c(j));
            partbothcCand(j)=yboth(i,j)*log(pd2ccand(j))+(K-yboth(i,j))*log(1-pd2ccand(j));
            partleftcCand(j)=yleft(i,j)*log(pd12ccand(j))+(K-yleft(i,j))*log(1-pd12ccand(j));
            partrightcCand(j)=yright(i,j)*log(pd12ccand(j))+(K-yright(i,j))*log(1-pd12ccand(j));
          }else{ //if single trap
            partbothc(j)=0;
            partleftc(j)=yleft(i,j)*log(pd1c(j))+(K-yleft(i,j))*log(1-pd1c(j));
            partrightc(j)=yright(i,j)*log(pd1c(j))+(K-yright(i,j))*log(1-pd1c(j));
            partleftcCand(j)=yleft(i,j)*log(pd1ccand(j))+(K-yleft(i,j))*log(1-pd1ccand(j));
            partrightcCand(j)=yright(i,j)*log(pd1ccand(j))+(K-yright(i,j))*log(1-pd1ccand(j));
          }
          if(partbothc(j)==partbothc(j)){
            v+=partbothc(j);
          }
          if(partleftc(j)==partleftc(j)){
            v+=partleftc(j);
          }
          if(partrightc(j)==partrightc(j)){
            v+=partrightc(j);
          }
          if(partbothcCand(j)==partbothcCand(j)){
            vcand+=partbothcCand(j);
          }
          if(partleftcCand(j)==partleftcCand(j)){
            vcand+=partleftcCand(j);
          }
          if(partrightcCand(j)==partrightcCand(j)){
            vcand+=partrightcCand(j);
          }
        }
        //Calculate prior for original and candidate
        if(grid(scell(i)-1,2)==2){
          priors=beta0 + beta1;
        }else{
          priors=beta0;
        }
        if(grid(Scandcell(0)-1,2)==2){
          priorscand=beta0 + beta1;
        }else{
          priorscand=beta0;
        }
        rand=Rcpp::runif(1);
        explldiff=exp((vcand+priorscand)-(v+priors));
        if((rand(0)<explldiff)){
          scell(i)=Scandcell(0);
          s(i,0)=ScandX(0);
          s(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D(i,j) = dtmp(j);
            lamd1(i,j) = lamd1candc(j);
            lamd2(i,j) = lamd2candc(j);
            }
        }
      }
    }

    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      sxout(iteridx,_)= s(_,0);
      syout(iteridx,_)= s(_,1);
      zout(iteridx,_)= z;
      ID_Lout(iteridx,_)=ID_L;
      ID_Rout(iteridx,_)=ID_R;
      out(iteridx,0)=lam01;
      out(iteridx,1)=lam02;
      out(iteridx,2)=sigma;
      out(iteridx,3)=N;
      out(iteridx,4)=EN;
      out(iteridx,5)=beta0;
      out(iteridx,6)=beta1;
      iteridx=iteridx+1;
    }
  }

  List to_return(6);
  to_return[0] = out;
  to_return[1] = sxout;
  to_return[2] = syout;
  to_return[3] = ID_Lout;
  to_return[4] = ID_Rout;
  to_return[5] = zout;
  return to_return;
}


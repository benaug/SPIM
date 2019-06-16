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

////////////////////////////////////////Basic SCR model///////////////////////////////////////////////
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMC1(double lam0, double sigma,NumericMatrix y, NumericMatrix lamd, IntegerVector z, NumericMatrix X,int K,NumericMatrix D, IntegerVector knownvector,
           NumericMatrix s,NumericVector psi, NumericVector xlim,NumericVector ylim,bool useverts, NumericMatrix vertices, double proplam0,
           double propsigma,double propsx,double propsy, int niter, int nburn, int nthin,
           int obstype,IntegerVector tf,bool storeLatent) {
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
  int nstore2=nstore;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,3);
  if(!storeLatent){
    nstore2=1;
  }
  NumericMatrix sxout(nstore2,M);
  NumericMatrix syout(nstore2,M);
  NumericMatrix zout(nstore2,M);
  int iteridx=0;
  //Starting log likelihood obs mod
  llysum=0;
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      if(obstype==1){
        pd(i,j)=1-exp(-lamd(i,j));
        if(z(i)==0){
          ll_y_curr(i,j)=0;
        }else{
          ll_y_curr(i,j)=R::dbinom(y(i,j),tf(j),pd(i,j),TRUE);
        }
      }else{
        if(z(i)==0){
          ll_y_curr(i,j)=0;
        }else{
          ll_y_curr(i,j)=R::dpois(y(i,j),tf(j)*lamd(i,j),TRUE);
        }
      }	  
      llysum+=ll_y_curr(i,j); 
    }
  }
  //start MCMC here
  int iter;
  for(iter=0; iter<niter; iter++){
    // Need to resum the ll_y on each iter
    llysum=0;
    for(int i=0; i<M; i++) {
      for(int j=0; j<J; j++){
        llysum+=ll_y_curr(i,j);
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
            if(z(i)==0){
              ll_y_cand(i,j)=0;
            }else{
              ll_y_cand(i,j)=R::dbinom(y(i,j),tf(j),pdcand(i,j),TRUE);
            }
          }else{
            if(z(i)==0){
              ll_y_cand(i,j)=0;
            }else{
              ll_y_cand(i,j)=R::dpois(y(i,j),tf(j)*lamdcand(i,j),TRUE);
            }
          }
          llycandsum+=ll_y_cand(i,j);
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
            if(z(i)==0){
              ll_y_cand(i,j)=0;
            }else{
              ll_y_cand(i,j)=R::dbinom(y(i,j),tf(j),pdcand(i,j),TRUE);
            }
          }else{
            if(z(i)==0){
              ll_y_cand(i,j)=0;
            }else{
              ll_y_cand(i,j)=R::dpois(y(i,j),tf(j)*lamdcand(i,j),TRUE);
            }
          }
          llycandsum+=ll_y_cand(i,j);
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
          if(z(i)==0){
            ll_y_curr(i,j)=0;
          }else{
            ll_y_curr(i,j)=R::dbinom(y(i,j),tf(j),pd(i,j),TRUE);
          }
        }else{
          if(z(i)==0){
            ll_y_curr(i,j)=0;
          }else{
            ll_y_curr(i,j)=R::dpois(y(i,j),tf(j)*lamd(i,j),TRUE);
          }
        }	  
        llysum+=ll_y_curr(i,j); 
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
        }
        if(obstype==1){
          for(int j=0; j<J; j++){
            pdcand(i,j)=1-exp(-lamdcand(i,j));
            if(z(i)==0){
              ll_y_cand(i,j)=0;
            }else{
              ll_y_cand(i,j)=R::dbinom(y(i,j),tf(j),pdcand(i,j),TRUE);
            }
          }
        }else{
          for(int j=0; j<J; j++){
            if(z(i)==0){
              ll_y_curr(i,j)=0;
            }else{
              ll_y_curr(i,j)=R::dpois(y(i,j),tf(j)*lamdcand(i,j),TRUE);
            }
          }
        }
        for(int j=0; j<J; j++){
          llycandsum+=ll_y_cand(i,j);
          llysum+=ll_y_curr(i,j);
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
      if(storeLatent){
        sxout(iteridx,_)= s(_,0);
        syout(iteridx,_)= s(_,1);
        zout(iteridx,_)= z;
      }
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

////////////////////////////////////////Side matching model///////////////////////////////////////////
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List MCMC2side(double lam01,double lam02, double sigma,NumericMatrix lamd1,NumericMatrix lamd2,
               NumericMatrix y_both, NumericMatrix y_left_true, NumericMatrix y_right_true,
               NumericMatrix y_left_obs, NumericMatrix y_right_obs,
               IntegerVector z, NumericMatrix X,IntegerVector tf,
               NumericMatrix D,int Nfixed, IntegerVector knownvector,IntegerVector ID_L, IntegerVector ID_R,int swap,double swaptol,
               NumericMatrix s,NumericVector psi, NumericVector xlim,NumericVector ylim,
               bool useverts, NumericMatrix vertices, double proplam01, double proplam02, double propsigma, double propsx, double propsy,
               int niter, int nburn, int nthin, LogicalVector updates,bool storeLatent) {
  RNGScope scope;
  int M = z.size();
  int J=y_both.ncol();
  //Preallocate for update lam01 lam02 sigma
  double  lam01cand;
  double  lam02cand;
  double  sigmacand;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd1cand(M,J);
  NumericMatrix lamd2cand(M,J);
  NumericMatrix pd1cand(M,J);
  NumericMatrix pd12cand(M,J);
  NumericMatrix pd2cand(M,J);
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix  ll_y_both_curr(M,J);
  NumericMatrix  ll_y_left_curr(M,J);
  NumericMatrix  ll_y_right_curr(M,J);
  NumericMatrix  ll_y_both_cand(M,J);
  NumericMatrix  ll_y_left_cand(M,J);
  NumericMatrix  ll_y_right_cand(M,J);
  double llysum;
  double llycandsum;
  double lly1sum;
  double lly1candsum;
  double lly2sum;
  double lly2candsum;
  //Preallocate side swapping structures
  //Swapping structures
  IntegerVector IDs=seq_len(M);
  int guy1=0;
  int guy2=0;
  IntegerVector swapped(2);
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
  NumericMatrix y_left_tmp(M,J);
  NumericMatrix y_right_tmp(M,J);  //could use same structure but dont want to get confused
  
  //Housekeeping structures
  NumericVector guycounts(M);
  LogicalVector zeroguys(M);
  
  //Preallocate for Psi update
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericVector prob0(M);
  NumericMatrix pbar1(M,J);
  NumericMatrix pbar2(M,J);
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
  NumericMatrix out(nstore,5);
  int nstore2=nstore;
  if(!storeLatent){
    nstore2=1;
  }
  NumericMatrix sxout(nstore2,M);
  NumericMatrix syout(nstore2,M);
  NumericMatrix zout(nstore2,M);
  NumericMatrix ID_Lout(nstore2,M);
  NumericMatrix ID_Rout(nstore2,M);
  int iter=0;
  int iteridx=0;
  //Starting log likelihood obs mod
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      pd1(i,j)=1-exp(-lamd1(i,j));
      pd2(i,j)=1-exp(-lamd2(i,j));
      if(X(j,2)==2){ //if double trap
        pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
        ll_y_both_curr(i,j)=R::dbinom(y_both(i,j),tf(j),z(i)*pd2(i,j),TRUE);
        ll_y_left_curr(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd12(i,j),TRUE);
        ll_y_right_curr(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd12(i,j),TRUE);
      }else{//single trap
        ll_y_both_curr(i,j)=0;
        ll_y_left_curr(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd1(i,j),TRUE);
        ll_y_right_curr(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd1(i,j),TRUE);
      }
    }
  }
  
  //start for loop here
  for(iter=0; iter<niter; iter++){
    //Sum single and both side LL on each iteration
    lly1sum=0;
    lly2sum=0;
    for(int i=0; i<M; i++) {
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          lly2sum+=ll_y_both_curr(i,j); 
        }
        lly1sum+=ll_y_left_curr(i,j); 
        lly1sum+=ll_y_right_curr(i,j); 
      }
    }
    //Update lam01
    if(updates(0)){
      rand=Rcpp::rnorm(1,lam01,proplam01);
      if(rand(0) > 0){
        lam01cand=rand(0);
        lly1candsum=0;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd1cand(i,j)=lam01cand*exp(-D(i,j)*D(i,j)/(2*sigma*sigma));
            pd1cand(i,j)=1-exp(-lamd1cand(i,j));
            if(X(j,2)==2){ //if double trap
              pd12cand(i,j)=2*pd1cand(i,j)-pd1cand(i,j)*pd1cand(i,j);
              ll_y_left_cand(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd12cand(i,j),TRUE);
              ll_y_right_cand(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd12cand(i,j),TRUE);
            }else{//single trap
              ll_y_left_cand(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd1cand(i,j),TRUE);
              ll_y_right_cand(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd1cand(i,j),TRUE);
            }
            lly1candsum+=ll_y_left_cand(i,j);
            lly1candsum+=ll_y_right_cand(i,j);
          }
        }
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(lly1candsum-lly1sum)){
          lam01=lam01cand;
          lly1sum=lly1candsum;
          for(int i=0; i<M; i++) {
            for(int j=0; j<J; j++){
              lamd1(i,j)=lamd1cand(i,j);
              pd1(i,j)=pd1cand(i,j);
              pd12(i,j)=pd12cand(i,j);
              ll_y_left_curr(i,j)=ll_y_left_cand(i,j);
              ll_y_right_curr(i,j)=ll_y_right_cand(i,j);
            }
          }
        }
      }
    }
    //Update lam02
    if(updates(1)){
      rand=Rcpp::rnorm(1,lam02,proplam02);
      if(rand(0) > 0){
        lam02cand=rand(0);
        lly2candsum=0;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd2cand(i,j)=lam02cand*exp(-D(i,j)*D(i,j)/(2*sigma*sigma));
            pd2cand(i,j)=1-exp(-lamd2cand(i,j));
            if(X(j,2)==2){ //if double trap
              ll_y_both_cand(i,j)=R::dbinom(y_both(i,j),tf(j),z(i)*pd2cand(i,j),TRUE);
              lly2candsum+=ll_y_both_cand(i,j);
            }
          }
        }
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(lly2candsum-lly2sum)){
          lly2sum=lly2candsum;
          lam02=lam02cand;
          for(int i=0; i<M; i++) {
            for(int j=0; j<J; j++){
              lamd2(i,j)=lamd2cand(i,j);
              pd2(i,j)=pd2cand(i,j);
              ll_y_both_curr(i,j)=ll_y_both_cand(i,j);
            }
          }
        }
      }
    }
    //update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lly1candsum=0;
      lly2candsum=0;
      for(int i=0; i<M; i++) {
        for(int j=0; j<J; j++){
          lamd1cand(i,j)=lam01*exp(-D(i,j)*D(i,j)/(2*sigmacand*sigmacand));
          lamd2cand(i,j)=lam02*exp(-D(i,j)*D(i,j)/(2*sigmacand*sigmacand));
          pd1cand(i,j)=1-exp(-lamd1cand(i,j));
          pd2cand(i,j)=1-exp(-lamd2cand(i,j));
          if(X(j,2)==2){ //if double trap
            ll_y_both_cand(i,j)=R::dbinom(y_both(i,j),tf(j),z(i)*pd2cand(i,j),TRUE);
            lly2candsum+=ll_y_both_cand(i,j);
            pd12cand(i,j)=2*pd1cand(i,j)-pd1cand(i,j)*pd1cand(i,j);
            ll_y_left_cand(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd12cand(i,j),TRUE);
            ll_y_right_cand(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd12cand(i,j),TRUE);
          }else{//single trap
            ll_y_left_cand(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd1cand(i,j),TRUE);
            ll_y_right_cand(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd1cand(i,j),TRUE);
          }
          lly1candsum+=ll_y_left_cand(i,j);
          lly1candsum+=ll_y_right_cand(i,j);
        }
      }
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp((lly1candsum+lly2candsum)-(lly1sum+lly2sum))){
        sigma=sigmacand;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd1(i,j)=lamd1cand(i,j);
            lamd2(i,j)=lamd2cand(i,j);
            pd1(i,j)=pd1cand(i,j);
            pd12(i,j)=pd12cand(i,j);
            pd2(i,j)=pd2cand(i,j);
            ll_y_both_curr(i,j)=ll_y_both_cand(i,j);
            ll_y_left_curr(i,j)=ll_y_left_cand(i,j);
            ll_y_right_curr(i,j)=ll_y_right_cand(i,j);
          }
        }
      }
    }
    //  Update left sides
    if(updates(2)){
      //   //Build mapL
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

        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            y_left_tmp(newID(i)-1,j)=y_left_obs(i,j);
          }
        }
        swapped(0)=swapout;
        swapped(1)=swapin;
        //Calculate likelihoods pre and post switch. Only need to calc for left data
        llysum=0;
        llycandsum=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            if(X(j,2)==2){ //if double trap
              ll_y_left_cand(swapped(i)-1,j)=R::dbinom(y_left_tmp(swapped(i)-1,j),tf(j),z(i)*pd12(swapped(i)-1,j),TRUE);
            }else{//single trap
              ll_y_left_cand(swapped(i)-1,j)=R::dbinom(y_left_tmp(swapped(i)-1,j),tf(j),z(i)*pd1(swapped(i)-1,j),TRUE);
            }
            llycandsum+=ll_y_left_cand(swapped(i)-1,j);
            llysum+=ll_y_left_curr(swapped(i)-1,j);
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        if(rand(0)<(exp(llycandsum-llysum)*(backprob/jumpprob))){
          y_left_true=Rcpp::clone(y_left_tmp);
          ID_L=Rcpp::clone(newID);
          mapL(guy1-1,1)=swapin;
          mapL(guy2-1,1)=swapout;
          for(int i=0; i<2; i++) {
            for(int j=0; j<J; j++){
              ll_y_left_curr(swapped(i)-1,j)=ll_y_left_cand(swapped(i)-1,j);
            }
          }
        }
      }
    }
    //Update Right sides
    if(updates(3)){
      //   //Build mapR
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

      ///////////////// //swap right sides//////////////
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

        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            y_right_tmp(newID(i)-1,j)=y_right_obs(i,j);
          }
        }
        swapped(0)=swapout;
        swapped(1)=swapin;
        //Calculate likelihoods pre and post switch. Only need to calc for right data
        llysum=0;
        llycandsum=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            if(X(j,2)==2){ //if double trap
              ll_y_right_cand(swapped(i)-1,j)=R::dbinom(y_right_tmp(swapped(i)-1,j),tf(j),z(i)*pd12(swapped(i)-1,j),TRUE);
            }else{//single trap
              ll_y_right_cand(swapped(i)-1,j)=R::dbinom(y_right_tmp(swapped(i)-1,j),tf(j),z(i)*pd1(swapped(i)-1,j),TRUE);
            }
            llycandsum+=ll_y_right_cand(swapped(i)-1,j);
            llysum+=ll_y_right_curr(swapped(i)-1,j);
          }
        }

        //MH step
        rand=Rcpp::runif(1);
        if(rand(0)<(exp(llycandsum-llysum)*(backprob/jumpprob))){
          y_right_true=Rcpp::clone(y_right_tmp);
          ID_R=Rcpp::clone(newID);
          mapR(guy1-1,1)=swapin;
          mapR(guy2-1,1)=swapout;
          for(int i=0; i<2; i++) {
            for(int j=0; j<J; j++){
              ll_y_right_curr(swapped(i)-1,j)=ll_y_right_cand(swapped(i)-1,j);
            }
          }
        }
      }
    }
    //Recalculate zeroguys
    for (int i=0; i<M; i++) {
      guycounts(i)=0;
      for(int j=0; j<J; j++) {
        guycounts(i)+=y_both(i,j)+y_left_true(i,j)+y_right_true(i,j);
      }
      zeroguys(i)=guycounts(i)==0;
    }
    //Update Psi
    // Calculate probability of no capture and update z
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          pbar1(i,j)=pow(1-pd12(i,j),2*tf(j));
          pbar2(i,j)=pow(1-pd2(i,j),tf(j));
          prob0(i)*=pbar1(i,j);
          prob0(i)*=pbar2(i,j);
        }else{//single trap
          pbar1(i,j)=pow(1-pd1(i,j),2*tf(j));
          prob0(i)*=pbar1(i,j);
        }
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0)&zeroguys(i);
      if(swappable(i)){
        NumericVector rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
    }
    
    //Calculate current likelihoods after z update
    for(int i=0; i<M; i++) {
      for(int j=0; j<J; j++){
        if(z(i)==0){
          ll_y_both_curr(i,j)=0;
          ll_y_left_curr(i,j)=0;
          ll_y_right_curr(i,j)=0;
        }else{
          if(X(j,2)==2){ //if double trap
            ll_y_both_curr(i,j)=R::dbinom(y_both(i,j),tf(j),z(i)*pd2(i,j),TRUE);
            ll_y_left_curr(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd12(i,j),TRUE);
            ll_y_right_curr(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd12(i,j),TRUE);
          }else{//single trap
            ll_y_both_curr(i,j)=0;
            ll_y_left_curr(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd1(i,j),TRUE);
            ll_y_right_curr(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd1(i,j),TRUE);
          }
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
        inbox=(ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      }else{
        inbox=inoutCpp(ScandX,ScandY,vertices);
      }
      if(inbox(0)){
        llysum=0;
        llycandsum=0;
        for(int j=0; j<J; j++){
          dtmp(j)=pow(pow(ScandX(0)-X(j,0), 2.0)+pow(ScandY(0)-X(j,1),2.0),0.5);
          lamd1cand(i,j)=lam01*exp(-dtmp(j)*dtmp(j)/(2*sigma*sigma));
          lamd2cand(i,j)=lam02*exp(-dtmp(j)*dtmp(j)/(2*sigma*sigma));
          pd1cand(i,j)=1-exp(-lamd1cand(i,j));
          pd2cand(i,j)=1-exp(-lamd2cand(i,j));
          if(X(j,2)==2){ //if double trap
            ll_y_both_cand(i,j)=R::dbinom(y_both(i,j),tf(j),z(i)*pd2cand(i,j),TRUE);
            llycandsum+=ll_y_both_cand(i,j);
            llysum+=ll_y_both_curr(i,j);
            pd12cand(i,j)=2*pd1cand(i,j)-pd1cand(i,j)*pd1cand(i,j);
            ll_y_left_cand(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd12cand(i,j),TRUE);
            ll_y_right_cand(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd12cand(i,j),TRUE);
          }else{//single trap
            ll_y_left_cand(i,j)=R::dbinom(y_left_true(i,j),tf(j),z(i)*pd1cand(i,j),TRUE);
            ll_y_right_cand(i,j)=R::dbinom(y_right_true(i,j),tf(j),z(i)*pd1cand(i,j),TRUE);
          }
          llycandsum+=ll_y_left_cand(i,j);
          llycandsum+=ll_y_right_cand(i,j);
          llysum+=ll_y_left_curr(i,j);
          llysum+=ll_y_right_curr(i,j);
        }
        rand=Rcpp::runif(1);
        if((rand(0)<exp(llycandsum-llysum))){
          s(i,0)=ScandX(0);
          s(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D(i,j) = dtmp(j);
            lamd1(i,j) = lamd1cand(i,j);
            lamd2(i,j) = lamd2cand(i,j);
            pd1(i,j) = pd1cand(i,j);
            pd12(i,j) = pd12cand(i,j);
            pd2(i,j) = pd2cand(i,j);
            ll_y_both_curr(i,j) = ll_y_both_cand(i,j);
            ll_y_left_curr(i,j) = ll_y_left_cand(i,j);
            ll_y_right_curr(i,j) = ll_y_right_cand(i,j);
          }
        }
      }
    }
    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      if(storeLatent){
        sxout(iteridx,_)= s(_,0);
        syout(iteridx,_)= s(_,1);
        zout(iteridx,_)= z;
        ID_Lout(iteridx,_)=ID_L;
        ID_Rout(iteridx,_)=ID_R;
      }
      out(iteridx,0)=lam01;
      out(iteridx,1)=lam02;
      out(iteridx,2)=sigma;
      out(iteridx,3)=N;
      out(iteridx,4)=psi(0);
      iteridx=iteridx+1;
    }
  }
  
  List to_return(9);
  to_return[0] = out;
  to_return[1] = sxout;
  to_return[2] = syout;
  to_return[3] = ID_Lout;
  to_return[4] = ID_Rout;
  to_return[5] = zout;
  to_return[6] = IDs;
  to_return[7] = lly1sum;
  to_return[8] = lly2sum;
  return to_return;
}

//2 side model with 2D trap file
//[[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h> included in openSCR file
using namespace Rcpp;
using namespace arma;
// #include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::export]]// [[Rcpp::export]]
List MCMCtf2(double lam01,double lam02,double sigma, NumericMatrix lamd1,
             NumericMatrix lamd2,NumericMatrix y_both, 
             arma::cube y_left_true, arma::cube y_right_true,
             arma::cube y_left_obs, arma::cube y_right_obs,IntegerVector z,
             NumericMatrix X,IntegerMatrix tf1,IntegerVector tf2,NumericMatrix D,int Nfixed, IntegerVector knownvector,IntegerVector ID_L,
             IntegerVector ID_R,int swap,double swaptol,NumericMatrix s,NumericVector psi,
             NumericVector xlim,NumericVector ylim,bool useverts, NumericMatrix vertices,double proplam01, double proplam02, double propsigma, double propsx,
             double propsy, int niter, int nburn, int nthin, LogicalVector updates,bool storeLatent) {
  RNGScope scope;
  int M = size(y_left_obs)[0];
  int J = size(y_left_obs)[1];
  int K = size(y_left_obs)[2];
  //Preallocate for update lam01 lam02 sigma
  double  lam01cand;
  double  lam02cand;
  double  sigmacand;
  NumericVector rand;
  NumericVector rand2;
  NumericMatrix lamd1cand(M,J);
  NumericMatrix lamd2cand(M,J);
  NumericMatrix pd1cand(M,J);
  NumericMatrix pd12cand(M,J);
  NumericMatrix pd2cand(M,J);
  NumericMatrix pd1(M,J);
  NumericMatrix pd12(M,J);
  NumericMatrix pd2(M,J);
  NumericMatrix  ll_y_both_curr(M,J);
  arma::cube  ll_y_left_curr=zeros<cube>(M,J,K);
  arma::cube  ll_y_right_curr=zeros<cube>(M,J,K);
  NumericMatrix  ll_y_both_cand(M,J);
  arma::cube  ll_y_left_cand=zeros<cube>(M,J,K);
  arma::cube  ll_y_right_cand=zeros<cube>(M,J,K);
  double llysum;
  double llycandsum;
  double lly1sum;
  double lly1candsum;
  double lly2sum;
  double lly2candsum;
  //Preallocate side swapping structures
  //Swapping structures
  IntegerVector IDs=seq_len(M);
  int guy1=0;
  int guy2=0;
  IntegerVector swapped(2);
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
  arma::cube  y_left_tmp=zeros<cube>(M,J,K);
  arma::cube  y_right_tmp=zeros<cube>(M,J,K);
  
  //Housekeeping structures
  NumericVector guycounts(M);
  LogicalVector zeroguys(M);
  
  //Preallocate for Psi update
  NumericVector fc(M);
  LogicalVector swappable(M);
  NumericVector prob0(M);
  arma::cube  pbar1=zeros<cube>(M,J,K);
  NumericMatrix pbar2(M,J);
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
  NumericMatrix out(nstore,5);
  int nstore2=nstore;
  if(!storeLatent){
    nstore2=1;
  }
  NumericMatrix sxout(nstore2,M);
  NumericMatrix syout(nstore2,M);
  NumericMatrix zout(nstore2,M);
  NumericMatrix ID_Lout(nstore2,M);
  NumericMatrix ID_Rout(nstore2,M);
  int iter=0;
  int iteridx=0;
  //Starting log likelihood obs mod
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      pd1(i,j)=1-exp(-lamd1(i,j));
      pd2(i,j)=1-exp(-lamd2(i,j));
      if(X(j,2)==2){ //if double trap
        pd12(i,j)=2*pd1(i,j)-pd1(i,j)*pd1(i,j);
        ll_y_both_curr(i,j)=R::dbinom(y_both(i,j),tf2(j),z(i)*pd2(i,j),TRUE);
      }
      for(int k=0; k<K; k++){
        if(tf1(j,k)==2){
          ll_y_left_curr(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd12(i,j),TRUE);
          ll_y_right_curr(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd12(i,j),TRUE);
        }else if(tf1(j,k)==1){
          ll_y_left_curr(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd1(i,j),TRUE);
          ll_y_right_curr(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd1(i,j),TRUE);
        }else{
          ll_y_left_curr(i,j,k)=0;
          ll_y_right_curr(i,j,k)=0;
        }
      }
    }
  }
  
  //start for loop here
  for(iter=0; iter<niter; iter++){
    //Sum single and both side LL on each iteration
    lly1sum=0;
    lly2sum=0;
    for(int i=0; i<M; i++) {
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          lly2sum+=ll_y_both_curr(i,j);
        }
        for(int k=0; k<K; k++){
          lly1sum+=ll_y_left_curr(i,j,k);
          lly1sum+=ll_y_right_curr(i,j,k);
        }
      }
    }
    //   //Update lam01
    if(updates(0)){
      rand=Rcpp::rnorm(1,lam01,proplam01);
      if(rand(0) > 0){
        lam01cand=rand(0);
        lly1candsum=0;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd1cand(i,j)=lam01cand*exp(-D(i,j)*D(i,j)/(2*sigma*sigma));
            pd1cand(i,j)=1-exp(-lamd1cand(i,j));
            if(X(j,2)==2){ //if double trap
              pd12cand(i,j)=2*pd1cand(i,j)-pd1cand(i,j)*pd1cand(i,j);
            }
            for(int k=0; k<K; k++){
              if(tf1(j,k)==2){
                ll_y_left_cand(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd12cand(i,j),TRUE);
                ll_y_right_cand(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd12cand(i,j),TRUE);
              }else if(tf1(j,k)==1){
                ll_y_left_cand(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd1cand(i,j),TRUE);
                ll_y_right_cand(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd1cand(i,j),TRUE);
              }else{
                ll_y_left_cand(i,j,k)=0;
                ll_y_right_cand(i,j,k)=0;
              }
              lly1candsum+=ll_y_left_cand(i,j,k);
              lly1candsum+=ll_y_right_cand(i,j,k);
            }
          }
        }
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(lly1candsum-lly1sum)){
          lam01=lam01cand;
          lly1sum=lly1candsum;
          for(int i=0; i<M; i++) {
            for(int j=0; j<J; j++){
              lamd1(i,j)=lamd1cand(i,j);
              pd1(i,j)=pd1cand(i,j);
              pd12(i,j)=pd12cand(i,j);
              for(int k=0; k<K; k++){
                ll_y_left_curr(i,j,k)=ll_y_left_cand(i,j,k);
                ll_y_right_curr(i,j,k)=ll_y_right_cand(i,j,k);
              }
            }
          }
        }
      }
    }
    //Update lam02
    if(updates(1)){
      rand=Rcpp::rnorm(1,lam02,proplam02);
      if(rand(0) > 0){
        lam02cand=rand(0);
        lly2candsum=0;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd2cand(i,j)=lam02cand*exp(-D(i,j)*D(i,j)/(2*sigma*sigma));
            pd2cand(i,j)=1-exp(-lamd2cand(i,j));
            if(X(j,2)==2){ //if double trap
              ll_y_both_cand(i,j)=R::dbinom(y_both(i,j),tf2(j),z(i)*pd2cand(i,j),TRUE);
              lly2candsum+=ll_y_both_cand(i,j);
            }
          }
        }
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(lly2candsum-lly2sum)){
          lly2sum=lly2candsum;
          lam02=lam02cand;
          for(int i=0; i<M; i++) {
            for(int j=0; j<J; j++){
              lamd2(i,j)=lamd2cand(i,j);
              pd2(i,j)=pd2cand(i,j);
              ll_y_both_curr(i,j)=ll_y_both_cand(i,j);
            }
          }
        }
      }
    }
    //update sigma
    rand=Rcpp::rnorm(1,sigma,propsigma);
    if(rand(0) > 0){
      sigmacand=rand(0);
      lly1candsum=0;
      lly2candsum=0;
      for(int i=0; i<M; i++) {
        for(int j=0; j<J; j++){
          lamd1cand(i,j)=lam01*exp(-D(i,j)*D(i,j)/(2*sigmacand*sigmacand));
          lamd2cand(i,j)=lam02*exp(-D(i,j)*D(i,j)/(2*sigmacand*sigmacand));
          pd1cand(i,j)=1-exp(-lamd1cand(i,j));
          pd2cand(i,j)=1-exp(-lamd2cand(i,j));
          if(X(j,2)==2){ //if double trap
            ll_y_both_cand(i,j)=R::dbinom(y_both(i,j),tf2(j),z(i)*pd2cand(i,j),TRUE);
            lly2candsum+=ll_y_both_cand(i,j);
            pd12cand(i,j)=2*pd1cand(i,j)-pd1cand(i,j)*pd1cand(i,j);
          }
          for(int k=0; k<K; k++){
            if(tf1(j,k)==2){
              ll_y_left_cand(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd12cand(i,j),TRUE);
              ll_y_right_cand(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd12cand(i,j),TRUE);
            }else if(tf1(j,k)==1){
              ll_y_left_cand(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd1cand(i,j),TRUE);
              ll_y_right_cand(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd1cand(i,j),TRUE);
            }else{
              ll_y_left_cand(i,j,k)=0;
              ll_y_right_cand(i,j,k)=0;
            }
            lly1candsum+=ll_y_left_cand(i,j,k);
            lly1candsum+=ll_y_right_cand(i,j,k);
          }
        }
      }
      rand2=Rcpp::runif(1);
      if(rand2(0)<exp((lly1candsum+lly2candsum)-(lly1sum+lly2sum))){
        sigma=sigmacand;
        for(int i=0; i<M; i++) {
          for(int j=0; j<J; j++){
            lamd1(i,j)=lamd1cand(i,j);
            lamd2(i,j)=lamd2cand(i,j);
            pd1(i,j)=pd1cand(i,j);
            pd12(i,j)=pd12cand(i,j);
            pd2(i,j)=pd2cand(i,j);
            ll_y_both_curr(i,j)=ll_y_both_cand(i,j);
            for(int k=0; k<K; k++){
              ll_y_left_curr(i,j,k)=ll_y_left_cand(i,j,k);
              ll_y_right_curr(i,j,k)=ll_y_right_cand(i,j,k);
            }
          }
        }
      }
    }
    //  Update left sides
    if(updates(2)){
      //   //Build mapL
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
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            for(int k=0; k<K; k++) {
              y_left_tmp(newID(i)-1,j,k)=y_left_obs(i,j,k);
            }
          }
        }
        swapped(0)=swapout;
        swapped(1)=swapin;
        //Calculate likelihoods pre and post switch. Only need to calc for left data
        llysum=0;
        llycandsum=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            for(int k=0; k<K; k++){
              if(tf1(j,k)==2){ //if double trap
                ll_y_left_cand(swapped(i)-1,j,k)=R::dbinom(y_left_tmp(swapped(i)-1,j,k),1,z(i)*pd12(swapped(i)-1,j),TRUE);
              }else if(tf1(j,k)==1){//single trap
                ll_y_left_cand(swapped(i)-1,j,k)=R::dbinom(y_left_tmp(swapped(i)-1,j,k),1,z(i)*pd1(swapped(i)-1,j),TRUE);
              }else{
                ll_y_left_cand(swapped(i)-1,j,k)=0;
              }
              llycandsum+=ll_y_left_cand(swapped(i)-1,j,k);
              llysum+=ll_y_left_curr(swapped(i)-1,j,k);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        if(rand(0)<(exp(llycandsum-llysum)*(backprob/jumpprob))){
          mapL(guy1-1,1)=swapin;
          mapL(guy2-1,1)=swapout;
          for (int i=0; i<M; i++) {//I think only need to loop over 2 guys
            ID_L(i)=newID(i);
            for(int j=0; j<J; j++) {
              for(int k=0; k<K; k++) {
                y_left_true(i,j,k)=y_left_tmp(i,j,k);
              }
            }
          }
          for(int i=0; i<2; i++) {
            for(int j=0; j<J; j++){
              for(int k=0; k<K; k++){
                ll_y_left_curr(swapped(i)-1,j,k)=ll_y_left_cand(swapped(i)-1,j,k);
              }
            }
          }
        }
      }
    }
    //Update Right sides
    if(updates(3)){
      //   //Build mapR
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
      
      ///////////////// //swap right sides//////////////
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
        for (int i=0; i<M; i++) {
          for(int j=0; j<J; j++) {
            for(int k=0; k<K; k++) {
              y_right_tmp(newID(i)-1,j,k)=y_right_obs(i,j,k);
            }
          }
        }
        swapped(0)=swapout;
        swapped(1)=swapin;
        //Calculate likelihoods pre and post switch. Only need to calc for left data
        llysum=0;
        llycandsum=0;
        for(int i=0; i<2; i++) {
          for(int j=0; j<J; j++){
            for(int k=0; k<K; k++){
              if(tf1(j,k)==2){ //if double trap
                ll_y_right_cand(swapped(i)-1,j,k)=R::dbinom(y_right_tmp(swapped(i)-1,j,k),1,z(i)*pd12(swapped(i)-1,j),TRUE);
              }else if(tf1(j,k)==1){//single trap
                ll_y_right_cand(swapped(i)-1,j,k)=R::dbinom(y_right_tmp(swapped(i)-1,j,k),1,z(i)*pd1(swapped(i)-1,j),TRUE);
              }else{
                ll_y_right_cand(swapped(i)-1,j,k)=0;
              }
              llycandsum+=ll_y_right_cand(swapped(i)-1,j,k);
              llysum+=ll_y_right_curr(swapped(i)-1,j,k);
            }
          }
        }
        //MH step
        rand=Rcpp::runif(1);
        if(rand(0)<(exp(llycandsum-llysum)*(backprob/jumpprob))){
          mapR(guy1-1,1)=swapin;
          mapR(guy2-1,1)=swapout;
          for (int i=0; i<M; i++) {//I think only need to loop over 2 guys
            ID_R(i)=newID(i);
            for(int j=0; j<J; j++) {
              for(int k=0; k<K; k++) {
                y_right_true(i,j,k)=y_right_tmp(i,j,k);
              }
            }
          }
          for(int i=0; i<2; i++) {
            for(int j=0; j<J; j++){
              for(int k=0; k<K; k++){
                ll_y_right_curr(swapped(i)-1,j,k)=ll_y_right_cand(swapped(i)-1,j,k);
              }
            }
          }
        }
      }
    }
    //Recalculate zeroguys
    for (int i=0; i<M; i++) {
      guycounts(i)=0;
      for(int j=0; j<J; j++) {
        guycounts(i)+=y_both(i,j);
        for(int k=0; k<K; k++) {
          guycounts(i)+=y_left_true(i,j,k)+y_right_true(i,j,k);
        }
      }
      zeroguys(i)=guycounts(i)==0;
    }
    //Update Psi
    // Calculate probability of no capture and update z
    for(int i=0; i<M; i++) {
      prob0(i)=1;
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          pbar2(i,j)=pow(1-pd2(i,j),tf2(j));
          prob0(i)*=pbar2(i,j);
        }
        for(int k=0; k<K; k++){
          if(tf1(j,k)==2){
            pbar1(i,j,k)=pow(1-pd12(i,j),2);
          }else if(tf1(j,k)==1){
            pbar1(i,j,k)=pow(1-pd1(i,j),2);
          }else{
            pbar1(i,j,k)=1;
          }
          prob0(i)*=pbar1(i,j,k);
        }
      }
      fc(i)=prob0(i)*psi(0)/(prob0(i)*psi(0) + 1-psi(0));
      swappable(i)=(knownvector(i)==0)&zeroguys(i);
      if(swappable(i)){
        NumericVector rand=Rcpp::rbinom(1,1,fc(i));
        z(i)=rand(0);
      }
    }
    
    //Calculate current likelihoods after z update
    for(int i=0; i<M; i++) {
      for(int j=0; j<J; j++){
        if(X(j,2)==2){ //if double trap
          ll_y_both_curr(i,j)=R::dbinom(y_both(i,j),tf2(j),z(i)*pd2(i,j),TRUE);
        }
        for(int k=0; k<K; k++){
          if(tf1(j,k)==2){
            ll_y_left_curr(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd12(i,j),TRUE);
            ll_y_right_curr(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd12(i,j),TRUE);
          }else if(tf1(j,k)==1){
            ll_y_left_curr(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd1(i,j),TRUE);
            ll_y_right_curr(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd1(i,j),TRUE);
          }else{
            ll_y_left_curr(i,j,k)=0;
            ll_y_right_curr(i,j,k)=0;
          }
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
        inbox=(ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
      }else{
        inbox=inoutCpp(ScandX,ScandY,vertices);
      }
      if(inbox(0)){
        llysum=0;
        llycandsum=0;
        for(int j=0; j<J; j++){
          dtmp(j)=pow(pow(ScandX(0)-X(j,0), 2.0)+pow(ScandY(0)-X(j,1),2.0),0.5);
          lamd1cand(i,j)=lam01*exp(-dtmp(j)*dtmp(j)/(2*sigma*sigma));
          lamd2cand(i,j)=lam02*exp(-dtmp(j)*dtmp(j)/(2*sigma*sigma));
          pd1cand(i,j)=1-exp(-lamd1cand(i,j));
          pd2cand(i,j)=1-exp(-lamd2cand(i,j));
          if(X(j,2)==2){ //if double trap
            ll_y_both_cand(i,j)=R::dbinom(y_both(i,j),tf2(j),z(i)*pd2cand(i,j),TRUE);
            llycandsum+=ll_y_both_cand(i,j);
            llysum+=ll_y_both_curr(i,j);
            pd12cand(i,j)=2*pd1cand(i,j)-pd1cand(i,j)*pd1cand(i,j);
          }
          for(int k=0; k<K; k++){
            if(tf1(j,k)==2){
              ll_y_left_cand(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd12cand(i,j),TRUE);
              ll_y_right_cand(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd12cand(i,j),TRUE);
            }else if(tf1(j,k)==1){//single trap
              ll_y_left_cand(i,j,k)=R::dbinom(y_left_true(i,j,k),1,z(i)*pd1cand(i,j),TRUE);
              ll_y_right_cand(i,j,k)=R::dbinom(y_right_true(i,j,k),1,z(i)*pd1cand(i,j),TRUE);
            }else{
              ll_y_left_cand(i,j,k)=0;
              ll_y_right_cand(i,j,k)=0;
            }
            llycandsum+=ll_y_left_cand(i,j,k);
            llycandsum+=ll_y_right_cand(i,j,k);
            llysum+=ll_y_left_curr(i,j,k);
            llysum+=ll_y_right_curr(i,j,k);
          }
        }
        rand=Rcpp::runif(1);
        if((rand(0)<exp(llycandsum-llysum))){
          s(i,0)=ScandX(0);
          s(i,1)=ScandY(0);
          for(int j=0; j<J; j++){
            D(i,j) = dtmp(j);
            lamd1(i,j) = lamd1cand(i,j);
            lamd2(i,j) = lamd2cand(i,j);
            pd1(i,j) = pd1cand(i,j);
            pd12(i,j) = pd12cand(i,j);
            pd2(i,j) = pd2cand(i,j);
            ll_y_both_curr(i,j) = ll_y_both_cand(i,j);
            for(int k=0; k<K; k++){
              ll_y_left_curr(i,j,k) = ll_y_left_cand(i,j,k);
              ll_y_right_curr(i,j,k) = ll_y_right_cand(i,j,k);
            }
          }
        }
      }
    }
    //Record output
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      if(storeLatent){
        sxout(iteridx,_)= s(_,0);
        syout(iteridx,_)= s(_,1);
        zout(iteridx,_)= z;
        ID_Lout(iteridx,_)=ID_L;
        ID_Rout(iteridx,_)=ID_R;
      }
      out(iteridx,0)=lam01;
      out(iteridx,1)=lam02;
      out(iteridx,2)=sigma;
      out(iteridx,3)=N;
      out(iteridx,4)=psi(0);
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
  to_return[6] = IDs;
  to_return[7] = y_left_tmp;
  return to_return;
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List HazACup(NumericMatrix y,NumericMatrix pd,NumericMatrix X,NumericMatrix D,NumericMatrix s,
             IntegerVector z,NumericMatrix ll_y_curr,
             double beta0, double beta1, double sigma, int M, int J,IntegerVector tf,double props,
             NumericVector xlim, NumericVector ylim) {
  NumericMatrix D_cand(M,J);
  NumericVector p0(J);//not providing to function will cause problem if no ACs are updated
  NumericVector haz(J);//same here. Shouldn't happen.
  NumericVector haz_cand(J);
  NumericVector p0_cand(J);
  NumericMatrix pd_cand(M,J);
  NumericMatrix ll_y_cand(M,J);
  double llysum;
  double llycandsum;
  NumericVector Scand(2);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  LogicalVector inbox(1);
  NumericVector rand;
  
  //Starting log likelihood obs mod
  llysum=0;
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++){
      llysum+=ll_y_curr(i,j); 
    }
  }
  for(int i=0; i<M; i++) {
    D_cand=clone(D);
    ScandX=Rcpp::rnorm(1,s(i,0),props);
    ScandY=Rcpp::rnorm(1,s(i,1),props);
    inbox= (ScandX<xlim[1]) & (ScandX>xlim[0]) & (ScandY<ylim[1]) & (ScandY>ylim[0]);
    if(inbox(0)){
      llycandsum=0;
      for(int j=0; j<J; j++){
        D_cand(i,j)=pow(pow(ScandX(0)-X(j,0), 2.0)+pow(ScandY(0)-X(j,1),2.0),0.5);
        haz_cand(j)=0;
        for(int i2=0; i2<M; i2++) {
          if(z(i2)==1){
            haz_cand(j)+=exp(-pow(D_cand(i2,j),2)/(2*pow(sigma,2)));
          }
        }
        p0_cand(j)=R::plogis(beta0+beta1*haz_cand(j),0,1,TRUE,FALSE);
        for(int i2=0; i2<M; i2++) {
          pd_cand(i2,j)=p0_cand(j)*exp(-pow(D_cand(i2,j),2)/(2*pow(sigma,2)));
          if(z(i2)==1){
            ll_y_cand(i2,j)=R::dbinom(y(i2,j),tf(j),pd_cand(i2,j),TRUE);
            llycandsum+=ll_y_cand(i2,j); 
          }else{
            ll_y_cand(i2,j)=0;
          }
        }
      }
      rand=Rcpp::runif(1);
      if((rand(0)<exp(llycandsum-llysum))){
        s(i,0)=ScandX(0);
        s(i,1)=ScandY(0);
        D(i,_)=D_cand(i,_);
        pd=pd_cand;
        p0=p0_cand;
        haz=haz_cand;
        ll_y_curr=ll_y_cand;
        llysum=llycandsum;
      }
    }
  }
  List to_return(6);
  to_return[0] = s;
  to_return[1] = D;
  to_return[2] = haz;
  to_return[3] = p0;
  to_return[4] = pd;
  to_return[5] = ll_y_curr;
  return to_return;
}


////////////////////////////////////////Integrated Likelihood///////////////////////////////////////////////
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double intlikRcpp(NumericVector parm, NumericMatrix ymat,IntegerMatrix X, int K, NumericMatrix G,
           NumericMatrix D,int n){
  RNGScope scope;
  double part1;
  double part2;
  double LLout;
  double nG=G.nrow();
  int J=X.nrow();
  NumericMatrix Pm(J,nG);
  NumericMatrix probcap(J,nG);
  NumericVector lik_cond(nG);
  NumericVector nv(n+1,1.0);
  double lik_cond_sum;
  NumericVector lik_marg(n+1);
  double p0=1/(1+exp(-parm(0)));
  double sigma=exp(parm(1));
  double n0=exp(parm(2));
  //calculate probcap
  for (int j=0; j<J; j++) {
    for (int g=0; g<nG; g++){
      probcap(j,g)=p0*exp(-1/(2*pow(sigma,2))*D(j,g)*D(j,g));
    }
  }
  nv(n)=n0;
  // calculate marginal likelihood
  for (int i=0; i<(n+1); i++){
    lik_cond_sum=0;
    for (int g=0; g<nG; g++){
      lik_cond(g)=0;
      for (int j=0; j<J; j++) {
        Pm(j,g)=R::dbinom(ymat(i,j),K,probcap(j,g),TRUE);
        lik_cond(g)+=Pm(j,g);
      }
      lik_cond_sum+=exp(lik_cond(g));
    }
    lik_marg(i)=lik_cond_sum*(1/nG);
  }
  part1=lgamma(n+n0 +1)-lgamma(n0+1);
  part2=0;
  for (int i=0; i<(n+1); i++){
    part2+=nv(i)*log(lik_marg(i));
  }
  LLout=-(part1+part2);
  double to_return;
  to_return = LLout;
  return to_return;
}

using namespace Rcpp;
// [[Rcpp::export]]
LogicalVector findPossible2D(IntegerVector z,IntegerMatrix G_true,
                           IntegerVector G_obs_true,int M,int ncat) {
  LogicalMatrix equals(M,ncat);
  LogicalVector allequals(M);
  for(int i=0; i<M; i++) {
    allequals(i)=TRUE;
    for(int l=0; l<ncat; l++) {
      equals(i,l)=G_true(i,l)==G_obs_true(l);
      allequals(i)=allequals(i)*equals(i,l);
    }
    allequals(i)=allequals(i)*(z(i)==1);
  }
  return allequals;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector arma_setdiff(arma::uvec& x, arma::uvec& y){
  x = arma::unique(x);
  y = arma::unique(y);
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    if (!q1.empty()) {
      x.shed_row(q1(0));
    }
  }
  Rcpp::NumericVector x2 = Rcpp::wrap(x);
  x2.attr("dim") = R_NilValue;
  return x2;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List indLL(IntegerVector z,arma::cube yI,arma::cube yP,double lambdaI,IntegerVector ID) {
  int M1 = size(yI)[0];
  int J = size(yI)[1];
  int K = size(yI)[2];
  arma::cube ll_yi=zeros(M1,J,K);
  double llsum=0;
  for(int i=0; i<M1; i++) {
    if(z(i)==1){
      for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
          ll_yi(i,j,k)=R::dpois(yI(i,j,k),yP(ID(i)-1,j,k)*lambdaI,TRUE);
          llsum+=ll_yi(i,j,k);
        }
      }
    }
  }
  List to_return(2);
  to_return[0] = ll_yi;
  to_return[1] = llsum;
  return to_return;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List siteup(IntegerVector z1,IntegerVector z2,arma::cube yI,arma::cube yP,arma::cube samp_mat,
            double lambdaI,IntegerVector packID,
            arma::cube ll_yi, arma::cube ll_yp,NumericMatrix lamd) {
  int M1 = size(yP)[0];
  int J = size(yP)[1];
  int K = size(yP)[2];
  int M2 = packID.length();
  LogicalVector match(M2);
  arma::cube ll_yi_cand = arma::zeros(M2,J,K);
  arma::cube ll_yp_cand = arma::zeros(M1,J,K);
  arma::cube yP_cand = arma::zeros(M1,J,K);
  NumericVector rand;
  
  for(int i=0; i<M1; i++) {
    if(z2(i)==1){
      int count=0;
      for(int l=0; l<M2; l++) {
        match(l)=packID(l)==i;
        if(match(l)==TRUE){
          count+=1;
        }
      }
      IntegerVector guys(count);
      int idx=0;
      for(int l=0; l<M2; l++) {
        if(match(l)==TRUE){
          guys(idx)=l;
          idx+=1;
        }
      }
      for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
          bool skip=FALSE;
          yP_cand(i,j,k)=yP(i,j,k)+samp_mat(i,j,k);
          if(yP_cand(i,j,k)<0){
            skip=TRUE;
          }
          if(skip==FALSE){
            int sumobs=0;
            for(int l=0; l<count; l++) {
              sumobs+=yI(guys(l),j,k);
            }
            if((sumobs>0)&(yP_cand(i,j,k)==0)){
              skip=TRUE;
            }
            if(skip==FALSE){//can propose this update
              double llyisum=0;
              double llyicandsum=0;
              for(int l=0; l<count; l++) {
                ll_yi_cand(guys(l),j,k)=R::dpois(yI(guys(l),j,k),z1(guys(l))*yP_cand(i,j,k)*lambdaI,TRUE);
                llyisum+=ll_yi(guys(l),j,k);
                llyicandsum+=ll_yi_cand(guys(l),j,k);
              }
              ll_yp_cand(i,j,k)=R::dpois(yP_cand(i,j,k),lamd(i,j),TRUE);
              rand=Rcpp::runif(1);
              if(rand(0) < exp((llyicandsum+ll_yp_cand(i,j,k))-
                 (llyisum+ll_yp(i,j,k)))){
                for(int l=0; l<count; l++) {
                  ll_yi(guys(l),j,k)=ll_yi_cand(guys(l),j,k);
                }
                yP(i,j,k)=yP_cand(i,j,k);
                ll_yp(i,j,k)=ll_yp_cand(i,j,k);
              }
            }
          }
        }
      }
    }else{//simulate from full conditional
      for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
          yP(i,j,k)=R::rpois(lamd(i,j));
          ll_yp(i,j,k)=R::dpois(yP(i,j,k),lamd(i,j),TRUE);
        }
      }
    }
  }
  List to_return(3);
  to_return[0] = ll_yi;
  to_return[1] = ll_yp;
  to_return[2] = yP;
  return to_return;
}


//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube calcllyevent(arma::cube y_true, arma::cube y_event, arma::cube Xcov3D,
                        arma::cube pi1,  arma::cube pi2,  arma::cube pi3, IntegerVector z,
                        int M, int J, int K) {
  NumericVector ptype(3);
  arma::cube ll_y_event=zeros<cube>(M,J,K);
  for(int i=0; i<M; i++) {
    if(z[i]==1){
      for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
          if((y_true(i,j,k)==1)&(Xcov3D(i,j,k)!=0)){
            ptype(0)=pi1(i,j,k);
            ptype(1)=pi2(i,j,k);
            ptype(2)=pi3(i,j,k);
            ll_y_event(i,j,k)=log(ptype(y_event(i,j,k)-1));
          }
        }
      }
    }
  }
  return ll_y_event;
}


//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix calcllyevent2D(IntegerMatrix y_true, IntegerMatrix y_event, IntegerMatrix Xcov3D,
                             NumericMatrix pi1,  NumericMatrix pi2,  NumericMatrix pi3,
                             int J, int K) {
  NumericVector ptype(3);
  NumericMatrix ll_y_event(J,K);
  for(int j=0; j<J; j++) {
    for(int k=0; k<K; k++) {
      if((y_true(j,k)==1)&(Xcov3D(j,k)!=0)){
        ptype(0)=pi1(j,k);
        ptype(1)=pi2(j,k);
        ptype(2)=pi3(j,k);
        ll_y_event(j,k)=log(ptype(y_event(j,k)-1));
      }
    }
  }
  return ll_y_event;
}

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List occDepUp(IntegerMatrix occidx, arma::cube pd, int upOcc,arma::cube y_true,
              arma::cube y_Occ_true,  arma::cube y_SCR,arma::cube y_event,
              arma::cube pi1,arma::cube pi2,arma::cube pi3,
              arma::cube ll_y,arma::cube ll_y_event,
              IntegerVector z,
              int M, int J, int K) {
  List out(5);
  NumericVector ptype(3);
  arma::cube y_Occ_true_cand=zeros<cube>(M,J,K);
  arma::cube y_cand=zeros<cube>(M,J,K);
  arma::cube y_event_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_event_cand=zeros<cube>(M,J,K);
  NumericVector probs(M);
  NumericVector normprobs(M);
  int thisj;
  int thisk;
  int lidx=occidx.nrow();
  IntegerVector cands=Rcpp::seq(0,M-1);
  IntegerVector choosei;
  int sumcount;
  NumericVector rand(1);
  // int l=0;
  for(int l=0; l<lidx; l++) { //for each occupancy observation
    thisj=occidx(l,0);
    thisk=occidx(l,1);
    //make distance based probs
    double sumprob=0;
    for(int i=0; i<M; i++){
      if(z(i)==1){
        probs(i)=pd(i,thisj,thisk);
        sumprob+=probs(i);
      }else{
        probs(i)=0;
      }
    }
    for(int i=0; i<M; i++){
      if(probs(i)!=0){
        normprobs(i)=probs(i)/sumprob;
      }else{
        normprobs(i)=0;
      }
    }
    //choose i index to propose update
    choosei=Rcpp::sample(cands,1,FALSE,normprobs);
    //update y_occ_true and y_event
    for(int i=0; i<M; i++){
      y_Occ_true_cand(i,thisj,thisk)=y_Occ_true(i,thisj,thisk);
      y_cand(i,thisj,thisk)=y_true(i,thisj,thisk);
    }
    if(y_Occ_true(choosei(0),thisj,thisk)==1){
      y_Occ_true_cand(choosei(0),thisj,thisk)=0;
      if(y_SCR(choosei(0),thisj,thisk)==0){
        y_cand(choosei(0),thisj,thisk)=0;//must be detected if site visited
      }
    }else{
      y_Occ_true_cand(choosei(0),thisj,thisk)=1;
      y_cand(choosei(0),thisj,thisk)=1;//must visit site to be detected
    }
    sumcount=0;
    for(int i=0; i<M; i++){
      sumcount+=y_Occ_true_cand(i,thisj,thisk);
    }
    if(sumcount>0){//can't be 0 if we saw someone there
      ptype(0)=pi1(choosei(0),thisj,thisk);
      ptype(1)=pi2(choosei(0),thisj,thisk);
      ptype(2)=pi3(choosei(0),thisj,thisk);
      //calculate event likelihood
      if((y_SCR(choosei(0),thisj,thisk)==1)&(y_Occ_true_cand(choosei(0),thisj,thisk)==0)){
        y_event_cand(choosei(0),thisj,thisk)=1;
        ll_y_event_cand(choosei(0),thisj,thisk)=log(ptype(0));
      }else if((y_SCR(choosei(0),thisj,thisk)==0)&(y_Occ_true_cand(choosei(0),thisj,thisk)==1)){
        y_event_cand(choosei(0),thisj,thisk)=2;
        ll_y_event_cand(choosei(0),thisj,thisk)=log(ptype(1));
      }else if((y_SCR(choosei(0),thisj,thisk)==1)&(y_Occ_true_cand(choosei(0),thisj,thisk)==1)){
        y_event_cand(choosei(0),thisj,thisk)=3;
        ll_y_event_cand(choosei(0),thisj,thisk)=log(ptype(2));
      }else{
        y_event_cand(choosei(0),thisj,thisk)=0;
        ll_y_event_cand(choosei(0),thisj,thisk)=0;
      }
      ll_y_cand(choosei(0),thisj,thisk)= R::dbinom(y_cand(choosei(0),thisj,thisk),1,pd(choosei(0),thisj,thisk),TRUE);
      double sumll=0;
      double sumllcand=0;
      sumll+=ll_y(choosei(0),thisj,thisk);
      sumll+=ll_y_event(choosei(0),thisj,thisk);
      sumllcand+=ll_y_cand(choosei(0),thisj,thisk);
      sumllcand+=ll_y_event_cand(choosei(0),thisj,thisk);
      rand=Rcpp::runif(1);
      if (rand(0) < exp(sumllcand - sumll)) {
        y_Occ_true(choosei(0),thisj,thisk)=y_Occ_true_cand(choosei(0),thisj,thisk);
        y_true(choosei(0),thisj,thisk)=y_cand(choosei(0),thisj,thisk);
        y_event(choosei(0),thisj,thisk)=y_event_cand(choosei(0),thisj,thisk);
        ll_y_event(choosei(0),thisj,thisk)=ll_y_event_cand(choosei(0),thisj,thisk);
        ll_y(choosei(0),thisj,thisk)=ll_y_cand(choosei(0),thisj,thisk);
      }
    }
  }
  out(0)=y_Occ_true;
  out(1)=y_true;
  out(2)=y_event;
  out(3)=ll_y;
  out(4)=ll_y_event;
  return out;
}

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List occDepUpMbradius(IntegerMatrix occidx, NumericMatrix s,NumericMatrix X,arma::cube pd, 
                      double upDist,arma::cube y_true,
                      arma::cube y_Occ_true,  arma::cube y_SCR,arma::cube y_event,
                      arma::cube pi1,arma::cube pi2,arma::cube pi3,
                      arma::cube ll_y,arma::cube ll_y_event,arma::cube state,
                      arma::cube linp_p0,arma::cube p0, NumericVector beta_p0,
                      IntegerVector sex, arma::cube Xcov3D1, arma::cube Xcov3D2,
                      arma::cube D,IntegerVector z,NumericVector sigma,
                      int M, int J, int K) {
  List out(5);
  NumericVector ptype(3);
  arma::cube y_Occ_true_cand=zeros<cube>(M,J,K);
  arma::cube y_cand=zeros<cube>(M,J,K);
  arma::cube y_event_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_event_cand=zeros<cube>(M,J,K);
  arma::cube state_cand=zeros<cube>(M,J,K);
  arma::cube linp_p0_cand=zeros<cube>(M,J,K);
  arma::cube p0_cand=zeros<cube>(M,J,K);
  arma::cube pd_cand=zeros<cube>(M,J,K);
  NumericVector dists(M);
  int thisj;
  int thisk;
  int lidx=occidx.nrow();
  IntegerVector cands=Rcpp::seq(0,M-1);
  IntegerVector choosei;
  int sumcount;
  NumericVector rand(1);
  // int l=0;
  for(int l=0; l<lidx; l++) { //for each occupancy observation
    thisj=occidx(l,0);
    thisk=occidx(l,1);
    for(int i=0; i<M; i++){
      if(z(i)==1){
        dists(i)=pow(pow(s(i,0)-X(thisj,0),2)+pow(s(i,1)-X(thisj,1),2),0.5);
        if(dists(i)<upDist){
          //update y_occ_true for all M to make sure you don' turn everyone off
          for(int i2=0; i2<M; i2++){
            y_Occ_true_cand(i2,thisj,thisk)=y_Occ_true(i2,thisj,thisk);
          }
          //fill in all k for y_cand to update Mb
          for(int k=0; k<K; k++){
            y_cand(i,thisj,k)=y_true(i,thisj,k);
          }
          if(y_Occ_true(i,thisj,thisk)==1){
            y_Occ_true_cand(i,thisj,thisk)=0;
            if(y_SCR(i,thisj,thisk)==0){
              y_cand(i,thisj,thisk)=0;//must be detected if site visited
            }
          }else{
            y_Occ_true_cand(i,thisj,thisk)=1;
            y_cand(i,thisj,thisk)=1;//must visit site to be detected
          }
          sumcount=0;
          for(int i2=0; i2<M; i2++){
            sumcount+=y_Occ_true_cand(i2,thisj,thisk);
          }
          if(sumcount>0){//can't be 0 if we saw someone there
            ptype(0)=pi1(i,thisj,thisk);
            ptype(1)=pi2(i,thisj,thisk);
            ptype(2)=pi3(i,thisj,thisk);
            //calculate event likelihood
            if((y_SCR(i,thisj,thisk)==1)&(y_Occ_true_cand(i,thisj,thisk)==0)){
              y_event_cand(i,thisj,thisk)=1;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(0));
            }else if((y_SCR(i,thisj,thisk)==0)&(y_Occ_true_cand(i,thisj,thisk)==1)){
              y_event_cand(i,thisj,thisk)=2;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(1));
            }else if((y_SCR(i,thisj,thisk)==1)&(y_Occ_true_cand(i,thisj,thisk)==1)){
              y_event_cand(i,thisj,thisk)=3;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(2));
            }else{
              y_event_cand(i,thisj,thisk)=0;
              ll_y_event_cand(i,thisj,thisk)=0;
            }
            //propose new behavioral state matrix
            state_cand(i,thisj,0)=0;
            for(int k=1; k<K; k++){
              if((y_cand(i,thisj,k-1)==1)|(state_cand(i,thisj,k-1)==1)){
                state_cand(i,thisj,k)=1;
              }else{
                state_cand(i,thisj,k)=0;
              }
            }
            //update likelihood at all k indices
            double sumll=0;
            double sumllcand=0;
            for(int k=0; k<K; k++){
              linp_p0_cand(i,thisj,k)=beta_p0(0)+beta_p0(1)*sex(i)+
                beta_p0(2)*Xcov3D1(i,thisj,k)+beta_p0(3)*Xcov3D2(i,thisj,k)+
                beta_p0(4)*state_cand(i,thisj,k);
              p0_cand(i,thisj,k)=1/(1+exp(-1*linp_p0_cand(i,thisj,k)));
              if(sex(i)==0){
                pd_cand(i,thisj,k)=p0_cand(i,thisj,k)*
                  exp(-pow(D(i,thisj,k),2)/(2*pow(sigma(0),2)));
              }else{
                pd_cand(i,thisj,k)=p0_cand(i,thisj,k)*
                  exp(-pow(D(i,thisj,k),2)/(2*pow(sigma(1),2)));
              }
              ll_y_cand(i,thisj,k)= R::dbinom(y_cand(i,thisj,k),1,pd_cand(i,thisj,k),TRUE);
              sumll+=ll_y(i,thisj,k);
              sumllcand+=ll_y_cand(i,thisj,k);
            }
            sumll+=ll_y_event(i,thisj,thisk);
            sumllcand+=ll_y_event_cand(i,thisj,thisk);
            
            rand=Rcpp::runif(1);
            if (rand(0) < exp(sumllcand - sumll)) {
              y_Occ_true(i,thisj,thisk)=y_Occ_true_cand(i,thisj,thisk);
              y_true(i,thisj,thisk)=y_cand(i,thisj,thisk);
              y_event(i,thisj,thisk)=y_event_cand(i,thisj,thisk);
              ll_y_event(i,thisj,thisk)=ll_y_event_cand(i,thisj,thisk);
              for(int k=0; k<K; k++){
                state(i,thisj,k)=state_cand(i,thisj,k);
                linp_p0(i,thisj,k)=linp_p0_cand(i,thisj,k);
                p0(i,thisj,k)=p0_cand(i,thisj,k);
                pd(i,thisj,k)=pd_cand(i,thisj,k);
                ll_y(i,thisj,k)=ll_y_cand(i,thisj,k);
              }
            }
          }
        }
      }
    }
  }
  out(0)=y_Occ_true;
  out(1)=y_true;
  out(2)=y_event;
  out(3)=ll_y;
  out(4)=y_cand;
  return out;
}

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List occDepUpMb4Dradius(IntegerMatrix occidx, NumericMatrix s,NumericMatrix X,arma::cube pd,
                        double upDist,arma::cube y_true,
                        arma::cube y_Occ_true,  arma::cube y_SCR,arma::cube y_event,
                        arma::cube pi1,arma::cube pi2,arma::cube pi3,
                        arma::cube ll_y,arma::cube ll_y_event,arma::cube state,
                        arma::cube linp_p0,arma::cube p0, NumericVector beta_p0,
                        IntegerVector sex, arma::cube Xcov3D1, arma::cube Xcov3D2,
                        arma::cube D,IntegerVector z,NumericVector sigma,
                        int M, int J, int K) {
  List out(9);
  NumericVector ptype(3);
  arma::cube y_Occ_true_cand=zeros<cube>(M,J,K);
  arma::cube y_cand=zeros<cube>(M,J,K);
  arma::cube y_event_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_event_cand=zeros<cube>(M,J,K);
  arma::cube state_cand=zeros<cube>(M,J,K);
  arma::cube linp_p0_cand=zeros<cube>(M,J,K);
  arma::cube p0_cand=zeros<cube>(M,J,K);
  arma::cube pd_cand=zeros<cube>(M,J,K);
  NumericVector dists(M);
  int thisj;
  int thisk;
  int lidx=occidx.nrow();
  IntegerVector cands=Rcpp::seq(0,M-1);
  IntegerVector choosei;
  int sumcount;
  NumericVector rand(1);
  // int l=0;
  for(int l=0; l<lidx; l++) { //for each occupancy observation
    thisj=occidx(l,0);
    thisk=occidx(l,1);
    for(int i=0; i<M; i++){
      if(z(i)==1){
        dists(i)=pow(pow(s(i,0)-X(thisj,0),2)+pow(s(i,1)-X(thisj,1),2),0.5);
        if(dists(i)<upDist){
          //update y_occ_true for all M to make sure you don't turn everyone off
          for(int i2=0; i2<M; i2++){
            y_Occ_true_cand(i2,thisj,thisk)=y_Occ_true(i2,thisj,thisk);
          }
          //fill in all k for y_cand to update Mb
          for(int k=0; k<K; k++){
            y_cand(i,thisj,k)=y_true(i,thisj,k);
          }
          if(y_Occ_true(i,thisj,thisk)==1){
            y_Occ_true_cand(i,thisj,thisk)=0;
            if(y_SCR(i,thisj,thisk)==0){
              y_cand(i,thisj,thisk)=0;//must be detected if site visited
            }
          }else{
            y_Occ_true_cand(i,thisj,thisk)=1;
            y_cand(i,thisj,thisk)=1;//must visit site to be detected
          }
          sumcount=0;
          for(int i2=0; i2<M; i2++){
            sumcount+=y_Occ_true_cand(i2,thisj,thisk);
          }
          if(sumcount>0){//can't be 0 if we saw someone there
            ptype(0)=pi1(i,thisj,thisk);
            ptype(1)=pi2(i,thisj,thisk);
            ptype(2)=pi3(i,thisj,thisk);
            //calculate event likelihood
            if((y_SCR(i,thisj,thisk)==1)&(y_Occ_true_cand(i,thisj,thisk)==0)){
              y_event_cand(i,thisj,thisk)=1;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(0));
            }else if((y_SCR(i,thisj,thisk)==0)&(y_Occ_true_cand(i,thisj,thisk)==1)){
              y_event_cand(i,thisj,thisk)=2;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(1));
            }else if((y_SCR(i,thisj,thisk)==1)&(y_Occ_true_cand(i,thisj,thisk)==1)){
              y_event_cand(i,thisj,thisk)=3;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(2));
            }else{
              y_event_cand(i,thisj,thisk)=0;
              ll_y_event_cand(i,thisj,thisk)=0;
            }
            //propose new behavioral state matrix
            state_cand(i,thisj,0)=0;
            for(int k=1; k<K; k++){
              if((y_cand(i,thisj,k-1)==1)|(state_cand(i,thisj,k-1)==1)){
                state_cand(i,thisj,k)=1;
              }else{
                state_cand(i,thisj,k)=0;
              }
            }
            //update likelihood at all k indices
            double sumll=0;
            double sumllcand=0;
            for(int k=0; k<K; k++){
              linp_p0_cand(i,thisj,k)=beta_p0(0)+beta_p0(1)*sex(i)+
                beta_p0(2)*Xcov3D1(i,thisj,k)+beta_p0(3)*Xcov3D2(i,thisj,k)+
                beta_p0(4)*state_cand(i,thisj,k);
              p0_cand(i,thisj,k)=1/(1+exp(-1*linp_p0_cand(i,thisj,k)));
              if(sex(i)==0){
                pd_cand(i,thisj,k)=p0_cand(i,thisj,k)*
                  exp(-pow(D(i,thisj,k),2)/(2*pow(sigma(0),2)));
              }else{
                pd_cand(i,thisj,k)=p0_cand(i,thisj,k)*
                  exp(-pow(D(i,thisj,k),2)/(2*pow(sigma(1),2)));
              }
              ll_y_cand(i,thisj,k)= R::dbinom(y_cand(i,thisj,k),1,pd_cand(i,thisj,k),TRUE);
              sumll+=ll_y(i,thisj,k);
              sumllcand+=ll_y_cand(i,thisj,k);
            }
            sumll+=ll_y_event(i,thisj,thisk);
            sumllcand+=ll_y_event_cand(i,thisj,thisk);
            
            rand=Rcpp::runif(1);
            if (rand(0) < exp(sumllcand - sumll)) {
              y_Occ_true(i,thisj,thisk)=y_Occ_true_cand(i,thisj,thisk);
              y_true(i,thisj,thisk)=y_cand(i,thisj,thisk);
              y_event(i,thisj,thisk)=y_event_cand(i,thisj,thisk);
              ll_y_event(i,thisj,thisk)=ll_y_event_cand(i,thisj,thisk);
              for(int k=0; k<K; k++){
                state(i,thisj,k)=state_cand(i,thisj,k);
                linp_p0(i,thisj,k)=linp_p0_cand(i,thisj,k);
                p0(i,thisj,k)=p0_cand(i,thisj,k);
                pd(i,thisj,k)=pd_cand(i,thisj,k);
                ll_y(i,thisj,k)=ll_y_cand(i,thisj,k);
              }
            }
          }
        }
      }
    }
  }
  out(0)=y_Occ_true;
  out(1)=y_true;
  out(2)=y_event;
  out(3)=ll_y;
  out(4)=ll_y_event;
  out(5)=state;
  out(6)=linp_p0;
  out(7)=p0;
  out(8)=pd;
  return out;
}

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List occDepUpRadius(IntegerMatrix occidx, NumericMatrix s,NumericMatrix X,
                    arma::cube pd, double upDist,arma::cube y_true,
                    arma::cube y_Occ_true,  arma::cube y_SCR,arma::cube y_event,
                    arma::cube pi1,arma::cube pi2,arma::cube pi3,
                    arma::cube ll_y,arma::cube ll_y_event,
                    IntegerVector z,
                    int M, int J, int K) {
  List out(5);
  NumericVector ptype(3);
  arma::cube y_Occ_true_cand=zeros<cube>(M,J,K);
  arma::cube y_cand=zeros<cube>(M,J,K);
  arma::cube y_event_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_cand=zeros<cube>(M,J,K);
  arma::cube ll_y_event_cand=zeros<cube>(M,J,K);
  NumericVector dists(M);
  int thisj;
  int thisk;
  int lidx=occidx.nrow();
  IntegerVector cands=Rcpp::seq(0,M-1);
  IntegerVector choosei;
  int sumcount;
  NumericVector rand(1);
  // int l=0;
  for(int l=0; l<lidx; l++) { //for each occupancy observation
    thisj=occidx(l,0);
    thisk=occidx(l,1);
    for(int i=0; i<M; i++){
      if(z(i)==1){
        dists(i)=pow(pow(s(i,0)-X(thisj,0),2)+pow(s(i,1)-X(thisj,1),2),0.5);
        if(dists(i)<upDist){
          //update y_occ_true and y_event
          for(int i2=0; i2<M; i2++){
            y_Occ_true_cand(i2,thisj,thisk)=y_Occ_true(i2,thisj,thisk);
            y_cand(i2,thisj,thisk)=y_true(i2,thisj,thisk);
          }
          if(y_Occ_true(i,thisj,thisk)==1){
            y_Occ_true_cand(i,thisj,thisk)=0;
            if(y_SCR(i,thisj,thisk)==0){
              y_cand(i,thisj,thisk)=0;//must be detected if site visited
            }
          }else{
            y_Occ_true_cand(i,thisj,thisk)=1;
            y_cand(i,thisj,thisk)=1;//must visit site to be detected
          }
          sumcount=0;
          for(int i2=0; i2<M; i2++){
            sumcount+=y_Occ_true_cand(i2,thisj,thisk);
          }
          if(sumcount>0){//can't be 0 if we saw someone there
            ptype(0)=pi1(i,thisj,thisk);
            ptype(1)=pi2(i,thisj,thisk);
            ptype(2)=pi3(i,thisj,thisk);
            //calculate event likelihood
            if((y_SCR(i,thisj,thisk)==1)&(y_Occ_true_cand(i,thisj,thisk)==0)){
              y_event_cand(i,thisj,thisk)=1;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(0));
            }else if((y_SCR(i,thisj,thisk)==0)&(y_Occ_true_cand(i,thisj,thisk)==1)){
              y_event_cand(i,thisj,thisk)=2;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(1));
            }else if((y_SCR(i,thisj,thisk)==1)&(y_Occ_true_cand(i,thisj,thisk)==1)){
              y_event_cand(i,thisj,thisk)=3;
              ll_y_event_cand(i,thisj,thisk)=log(ptype(2));
            }else{
              y_event_cand(i,thisj,thisk)=0;
              ll_y_event_cand(i,thisj,thisk)=0;
            }
            ll_y_cand(i,thisj,thisk)= R::dbinom(y_cand(i,thisj,thisk),1,pd(i,thisj,thisk),TRUE);
            double sumll=0;
            double sumllcand=0;
            sumll+=ll_y(i,thisj,thisk);
            sumll+=ll_y_event(i,thisj,thisk);
            sumllcand+=ll_y_cand(i,thisj,thisk);
            sumllcand+=ll_y_event_cand(i,thisj,thisk);
            rand=Rcpp::runif(1);
            if (rand(0) < exp(sumllcand - sumll)) {
              y_Occ_true(i,thisj,thisk)=y_Occ_true_cand(i,thisj,thisk);
              y_true(i,thisj,thisk)=y_cand(i,thisj,thisk);
              y_event(i,thisj,thisk)=y_event_cand(i,thisj,thisk);
              ll_y_event(i,thisj,thisk)=ll_y_event_cand(i,thisj,thisk);
              ll_y(i,thisj,thisk)=ll_y_cand(i,thisj,thisk);
            }
          }
        }
      }
    }
  }
  out(0)=y_Occ_true;
  out(1)=y_true;
  out(2)=y_event;
  out(3)=ll_y;
  out(4)=ll_y_event;
  return out;
}

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List typeCount(int ncat,int nrep,IntegerVector nlevels, List locustype,List ptype,
               IntegerMatrix G_obs_true, arma::cube G_obs) {
  List out(2);
  IntegerVector x_het(3);
  IntegerVector x_hom(2);
  int nobs=G_obs_true.nrow();
  int count;
  for(int l=0; l<ncat; l++) { 
    IntegerVector LT=locustype[l];
    IntegerMatrix PT=ptype[l];
    for(int k=0; k<(nlevels[l]); k++) { 
      if(LT(k)==1){//homozygote
        for(int k2=0; k2<(nlevels[l]); k2++) { 
          if(k==k2){//correct classification
            count=0;
            for(int m=0; m<nobs; m++) {
              for(int l2=0; l2<nrep; l2++) {
                if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                  count+=1;
                }
              }
            }
            x_hom(0)+=count;
          }else{//false allele
            count=0;
            for(int m=0; m<nobs; m++) {
              for(int l2=0; l2<nrep; l2++) {
                if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                  count+=1;
                }
              }
            }
            x_hom(1)+=count;
          }
        }
      }else{//heterozygote
        for(int k2=0; k2<(nlevels[l]); k2++) { 
          if(k==k2){//correct classification
            count=0;
            for(int m=0; m<nobs; m++) {
              for(int l2=0; l2<nrep; l2++) {
                if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                  count+=1;
                }
              }
            }
            x_het(0)+=count;
          }else{
            if(PT(k2,k)==2){//allelic dropout
              count=0;
              for(int m=0; m<nobs; m++) {
                for(int l2=0; l2<nrep; l2++) {
                  if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                    count+=1;
                  }
                }
              }
              x_het(1)+=count;
            }else{//false allele
              count=0;
              for(int m=0; m<nobs; m++) {
                for(int l2=0; l2<nrep; l2++) {
                  if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                    count+=1;
                  }
                }
              }
              x_het(2)+=count;
            }
          }
        }
      }
    }
  }
  out(0)=x_hom;
  out(1)=x_het;
  return out;
}
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List typeCountSamptype(int ncat,int nrep,IntegerVector nlevels, List locustype,List ptype,
                       IntegerMatrix G_obs_true, arma::cube G_obs,IntegerVector samptype,int samplevels) {
  List out(2);
  IntegerMatrix x_het(2,3);
  IntegerMatrix x_hom(2,2);
  int nobs=G_obs_true.nrow();
  IntegerVector count(3);
  for(int l=0; l<ncat; l++) { 
    IntegerVector LT=locustype[l];
    IntegerMatrix PT=ptype[l];
    for(int k=0; k<(nlevels[l]); k++) { 
      if(LT(k)==1){//homozygote
        for(int k2=0; k2<(nlevels[l]); k2++) { 
          if(k==k2){//correct classification
            for(int st=0; st<samplevels; st++) { 
              count(st)=0;
            }
            for(int m=0; m<nobs; m++) {
              for(int l2=0; l2<nrep; l2++) {
                if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                  count(samptype(m))+=1;
                }
              }
            }for(int st=0; st<samplevels; st++) { 
              x_hom(st,0)+=count(st);
            }
          }else{//false allele
            for(int st=0; st<samplevels; st++) { 
              count(st)=0;
            }
            for(int m=0; m<nobs; m++) {
              for(int l2=0; l2<nrep; l2++) {
                if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                  count(samptype(m))+=1;
                }
              }
            }for(int st=0; st<samplevels; st++) { 
              x_hom(st,1)+=count(st);
            }
          }
        }
      }else{//heterozygote
        for(int k2=0; k2<(nlevels[l]); k2++) { 
          if(k==k2){//correct classification
            for(int st=0; st<samplevels; st++) { 
              count(st)=0;
            }
            for(int m=0; m<nobs; m++) {
              for(int l2=0; l2<nrep; l2++) {
                if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                  count(samptype(m))+=1;
                }
              }
            }for(int st=0; st<samplevels; st++) { 
              x_het(st,0)+=count(st);
            }
          }else{
            if(PT(k2,k)==2){//allelic dropout
              for(int st=0; st<samplevels; st++) { 
                count(st)=0;
              }
              for(int m=0; m<nobs; m++) {
                for(int l2=0; l2<nrep; l2++) {
                  if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                    count(samptype(m))+=1;
                  }
                }
              }for(int st=0; st<samplevels; st++) { 
                x_het(st,1)+=count(st);
              }
            }else{//false allele
              for(int st=0; st<samplevels; st++) { 
                count(st)=0;
              }
              for(int m=0; m<nobs; m++) {
                for(int l2=0; l2<nrep; l2++) {
                  if((G_obs_true(m,l)==k)&(G_obs(m,l,l2)==k2)){
                    count(samptype(m))+=1;
                  }
                }
              }for(int st=0; st<samplevels; st++) { 
                x_het(st,2)+=count(st);
              }
            }
          }
        }
      }
    }
  }
  out(0)=x_hom;
  out(1)=x_het;
  return out;
}
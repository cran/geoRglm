#include <math.h> 
#include <stdio.h>
#include <memory.h>
#include <S.h>

#include "geoR.h"


#define FREE_ARG char*
# define Real double

#if(defined(SPLUS_VERSION))
# define Integer long
# define UNIF unif_rand(S_evaluator)
# define RNORM norm_rand(S_evaluator)
#else
# define Integer int
# define UNIF unif_rand()
# define RNORM norm_rand()
# define Salloc(n,t) (t*)S_alloc(n,sizeof(t))
#endif

#if(defined(SPLUS_VERSION) && SPLUS_VERSION > 5100) 
# define RANDIN seed_in((long *)NULL, S_evaluator)
# define RANDOUT seed_out((long *)NULL, S_evaluator)
#elif(SPLUS_VERSION > 5000)
# define RANDIN seedin((long *)NULL, S_evaluator)
# define RANDOUT seedout((long *)NULL, S_evaluator)
#else
# define RANDIN seed_in((long *)NULL)
# define RANDOUT seed_out((long *)NULL)
#endif



void binitprod(Integer *n, Real *xc, Real *yc, Real *sim, Integer *nbins, Real *lims, Real *maxdist, Integer *cbin, Real *vbin)
{
  Integer i, j, ind;
  Real register v;
  Real dist, dx, dy;
  
  ind =0; v=0;
  for (j=0; j < *n; j++)
  { 
      for (i=j+1; i<*n; i++) 
      {
	  dx = xc[i] - xc[j];
	  dy = yc[i] - yc[j];
	  dist = sqrt(dx*dx + dy*dy);
	  if(dist <= *maxdist)
	  {
	      v = sim[i]*sim[j];
	      ind = 0;
	      while (dist > lims[ind] && ind <= *nbins ) ind++ ;
	      if (dist <= lims[ind])
	      {
		  vbin[(ind-1)]+= v; 
		  cbin[(ind-1)]++;
	      }
	  }
      }
  }
  for (j=0; j < *nbins; j++) 
  {
      if (cbin[j]){
	vbin[j] = vbin[j]/cbin[j];
      }
  }
  
}

 
void cholesky(Real *inmat, Real *outmat, Integer n){
   /* returns L where L L'=inmatrix. This version is when input and output are just the lower triangles */
  
Integer i,j,k;
Real sum;

   for (i=0;i<n;i++) {
      for (j=i;j<n;j++) {
         for (sum=inmat[ n*i - i*(i+1)/2 + j],k=0;k<i;k++) sum -= outmat[ n*k - k*(k+1)/2 + i]*outmat[ n*k - k*(k+1)/2 + j];
         if (j == i) {
            if (sum < 0.0){
               printf("cholesky failed; inmat is not positive definite \n");  
               exit(1);
            }
            outmat[ n*i - i*(i+1)/2 + i]=sqrt(sum);
         }
         else outmat[ n*i - i*(i+1)/2 + j]=sum/outmat[ n*i - i*(i+1)/2 + i];
      }  
   }
}


Real trunc_u(Real mean, Real H){
  return (mean>H) ? H : mean;
}


void conddensity(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *linpred, Integer dim){
Integer l, k;

   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[l*dim+k]*z[k];
   }
   for (*logfcond=0,l=0; l<dim; l++){
     *logfcond+=obsdata[l]*(S[l]+linpred[l])-exp(S[l]+linpred[l]);
   }
}


void gradient(Real *S, Real *gradz, Real *BB, Real *z, Real *obsdata, Real *linpred, Real *Hz, Integer dim){
   Integer l, k;
   Real likeli;
   
   for (l=0; l<dim; l++) gradz[l]=-z[l];
   for (k=0; k<dim; k++){
      likeli=obsdata[k]- trunc_u(exp(S[k]+linpred[k]),Hz[k]);
      for (l=0; l<k+1; l++) gradz[l]+=BB[k*dim+l]*likeli;  
   }
}

Real calc_ss(Real *z, Integer dim){
Integer l;
Real sumsz;

   for (sumsz=0,l=0; l<dim; l++) sumsz+=pow(z[l],2);
   return sumsz;
}

void mcmcrun(Integer *n, Real *zz, Real *SS, Real *data, Real *meanS, Real *QQ, 
         Real *randnormal, Real *randunif, Real *Htrunc, Real *scale, Integer *nsim, Integer *subsample, Real *acc_rate){
         
  Integer i, ii, l, acc ;  
  Real *z = Salloc((*n),Real) ;
  Real *zprop = Salloc((*n),Real) ;
  Real *S = Salloc((*n),Real) ; 
  Real *Sprop = Salloc((*n),Real) ;
  Real *gradz = Salloc((*n),Real) ;
  Real *gradzprop = Salloc((*n),Real) ;
  Real *temp = Salloc((*n),Real) ;
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop  ; 
  
  for (l=0; l<(*n); l++){
      z[l] = zz[l]; 
      S[l] = 0;     
  } 
  conddensity(S,QQ,&logf,data,z,meanS,(*n));           
  gradient(S,gradz,QQ,z,data,meanS,Htrunc,(*n));
  logp_z = -calc_ss(z,(*n))/2;
  acc= 0;       
  for(i=0; i<(*nsim); i++){ 
     for(ii=0; ii<(*subsample); ii++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + randnormal[(i*(*subsample)+ii)*(*n)+l];
        conddensity(Sprop,QQ,&logfprop,data,zprop,meanS,(*n));
        gradient(Sprop,gradzprop,QQ,zprop,data,meanS,Htrunc,(*n));
        logp_zprop=-calc_ss(zprop,(*n))/2;  
        for (logq=0,logqprop=0,l=0; l<(*n); l++){
           logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
           logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
        }
        logq*=(-0.5/(*scale));
        logqprop*=(-0.5/(*scale));
        if (log(randunif[i*(*subsample)+ii])<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	   /*accept */
            logf=logfprop;
            logp_z=logp_zprop;
            temp=gradz;
            gradz=gradzprop;
            gradzprop=temp;
            temp=S;
            S=Sprop;
            Sprop=temp;
            temp=z;
            z=zprop;
            zprop=temp;
            acc++;
        }    
     } 
     for (l=0; l<(*n); l++) SS[i*(*n)+l] = S[l];  
  }
  *acc_rate = (Real) acc/((*nsim)*(*subsample));
  for (l=0; l<(*n); l++) zz[l] = z[l]; 
}

Real calc2_ss(Real *z, Real *Dmat, Integer dim){
Integer l,k;
Real sumsz,temp;
   for (sumsz=0,l=0; l<dim; l++){
       for (temp=0,k=0; k<l; k++)
            temp+=Dmat[l*dim+k]*z[k];
       sumsz+=(z[l]*temp*2+pow(z[l],2)*Dmat[l*dim+l]);
   }
   return sumsz;
}

void gradient2(Real *S, Real *gradz, Real *BB, Real *Dmat, Real *z, Real *obsdata, Real *linpred, Real *Hz, Integer dim){
   Integer l, k;
   Real likeli ;
   
   for (l=0; l<dim; l++) gradz[l]=0;
   for (k=0; k<dim; k++){
      likeli=obsdata[k]- trunc_u(exp(S[k]+linpred[k]),Hz[k]);
      for (l=0; l<dim; l++)
          if(l<k+1) gradz[l]+=BB[k*dim+l]*likeli - Dmat[l*dim+k]*z[k] ;
          else gradz[l]-=Dmat[k*dim+l]*z[k] ;
   }
}

void mcmcrun2(Integer *n, Real *zz, Real *SS, Real *data, Real *meanS, Real *QQ, Real *Dmat,
         Real *randnormal, Real *randunif, Real *Htrunc, Real *scale, Integer *nsim ,Integer *subsample, Real *acc_rate){
         
  Integer i, ii, l, acc ;  
  Real *z = Salloc((*n),Real) ; 
  Real *zprop = Salloc((*n),Real) ; 
  Real *S = Salloc((*n),Real) ;
  Real *Sprop = Salloc((*n),Real) ;
  Real *gradz = Salloc((*n),Real) ;
  Real *gradzprop = Salloc((*n),Real) ;
  Real *temp = Salloc((*n),Real) ;
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop  ;   
  for (l=0; l<(*n); l++){
      z[l] = zz[l];     
      S[l] = 0; 
  } 
  conddensity(S,QQ,&logf,data,z,meanS,(*n));          
  gradient2(S,gradz,QQ,Dmat,z,data,meanS,Htrunc,(*n));
  logp_z = -calc2_ss(z,Dmat,(*n))/2;
  acc= 0; 
  for(i=0; i<(*nsim); i++){ 
     for(ii=0; ii<(*subsample); ii++){ 
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + randnormal[(i*(*subsample)+ii)*(*n)+l];
        conddensity(Sprop,QQ,&logfprop,data,zprop,meanS,(*n));
        gradient2(Sprop,gradzprop,QQ,Dmat,zprop,data,meanS,Htrunc,(*n));
        logp_zprop=-calc2_ss(zprop,Dmat,(*n))/2;  
        for (logq=0,logqprop=0,l=0; l<(*n); l++){
           logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
           logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
        }
        logq*=(-0.5/(*scale));
        logqprop*=(-0.5/(*scale));      
        if (log(randunif[i*(*subsample)+ii])<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	   /*accept */
            logf=logfprop;
            logp_z=logp_zprop;
            temp=gradz;
            gradz=gradzprop;
            gradzprop=temp;
            temp=S;
            S=Sprop;
            Sprop=temp;
            temp=z;
            z=zprop;
            zprop=temp;
            acc++;
        }   
     }     
     for (l=0; l<(*n); l++) SS[i*(*n)+l] = S[l];     
  }
  *acc_rate = (Real) acc/((*nsim)*(*subsample));  
  for (l=0; l<(*n); l++) zz[l] = z[l]; 
}

void conddensity4(Real *S, Real *B, Real *logfcond, Real *obsdata, Real *z, Real *units, Integer dim){
Integer l, k;
           /* this version requires B to consist of just the entries in the lower triangular part of the matrix */          
   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=B[dim*k - k*(k+1)/2 + l]*z[k];
   }
   for (*logfcond=0,l=0; l<dim; l++){
     *logfcond+=(obsdata[l]*S[l]-units[l]*exp(S[l]));
   }
} 


Real calc4_ss(Real *z, Real *DDmat, Integer dim){
Integer l,k;
Real sumsz,temp;
     /* this version requires DDmat to consist of just the entries in the lower triangular part of the matrix */  
   for (sumsz=0,l=0; l<dim; l++){
       for (temp=0,k=0; k<l; k++)
            temp+=DDmat[dim*k - k*(k+1)/2 + l]*z[k];
       sumsz+=(z[l]*temp*2+pow(z[l],2)*DDmat[dim*l - l*(l+1)/2 + l]);
   }
   return sumsz;
}  


void gradient4(Real *S, Real *gradz, Real *B, Real *DDmat, Real *z, Real *obsdata, Real *units, Real *Hz, Integer dim, Real ss, Integer df){
   Integer l, k;
   Real likeli ;
     /* this version requires B and DDmat to consist of just the entries in the lower triangular part of the matrix */    
   for (l=0; l<dim; l++) gradz[l]=0;
   for (k=0; k<dim; k++){
      likeli=obsdata[k] - trunc_u(units[k]*exp(S[k]),Hz[k]);
      for (l=0; l<dim; l++){
          if(l<k+1) gradz[l]+=B[dim*l-l*(l+1)/2+k]*likeli - df*DDmat[dim*l-l*(l+1)/2+k]*z[k]/ss ; /* using symmetri of DDmat  */
          else gradz[l] -= df*DDmat[dim*k-k*(k+1)/2+l]*z[k]/ss ;   
      }  
   }
}

	  
Real logprior_phi(Real phi, Real e_mean, Integer phinum){
  switch(phinum){
  case 1:         /* uniform */  
    return 0;
    break;
  case 2:         /* exponential */  
    return exp(-phi/e_mean);
    break;
  case 3:         /* fixed */  
    printf("updating phi is not possible when phi is fixed \n");
    exit(1);
    break;
  case 4:         /* squared.reciprocal */  
    return 1/(phi*phi);
    break;
  case 5:         /* reciprocal */  
    return 1/phi;
    break;
  default: 
    return -1;
    break;
  }
  
}


void calc_Dmat(Real *B, Real *D, Real *outmat, Real *det_DivD_half, Integer dim, Integer no_linpar, Real *sqivD, Real *DivD,Real *sqDivD, Real *temp){   
/* calculation of the matrix Dmat */
/* this version is operating on the lower triangle of B and outmat (B is lower triangular, and outmat (Dmat) is symmetric) */

typedef Real *Doublearray; 
Integer l, k, q, p ;
Real sum, prod;

  /* computing sqivD = B^{-1}*D */
  for (q=0; q<no_linpar; q++){
      for (l=0; l<dim; l++){ 
          for (k=0,sum = 0; k<l; k++)
              sum += B[dim*k - k*(k+1)/2+l] *sqivD[k*no_linpar+q];
          sqivD[l*no_linpar+q] = (D[l*no_linpar+q]-sum)/B[dim*l - l*(l+1)/2+l];
      }
  }  
  /* computing DivD=sqivD^{T}*sqivD ( = D^{T}*(B*B^{T})^{-1}*D ) */
  for (q=0; q<no_linpar; q++){
      for (p=0; p<q+1; p++){
          for (l=0,DivD[no_linpar*p-p*(p+1)/2+q]=0; l<dim; l++)
              DivD[no_linpar*p-p*(p+1)/2+q]+=sqivD[l*no_linpar+q]*sqivD[l*no_linpar+p];  
      }
  } 
  if(no_linpar == 1){   
      for (l=0; l<dim; l++){
          for (k=0; k<l+1; k++){
              if(k==l) outmat[dim*l-l*(l+1)/2+l] = 1-sqivD[l*no_linpar]*sqivD[k*no_linpar]/DivD[0];
              else outmat[dim*k-k*(k+1)/2+l] = -sqivD[l*no_linpar]*sqivD[k*no_linpar]/DivD[0];
          }
      }  
  }
  else{
      cholesky(DivD,sqDivD,no_linpar); 
      /* temp =  sqDivD^{-1}*sqivD^{T}  */     
      for (q=0; q<no_linpar; q++){   
          for (l=0; l<dim; l++){  
              for (p=0,sum=0; p<q; p++)
                  sum += sqDivD[p*no_linpar-p*(p+1)/2+q]*temp[p*dim+l];
              temp[q*dim+l] = (sqivD[l*no_linpar+q]-sum)/sqDivD[q*no_linpar-q*(q+1)/2+q];  
          }
      } 
      /* computing outmat = I - temp^{T}*temp  ( =  I -  B^{-1}*D*(D^{T}*(B*B^{T})^{-1}*D)^{-1}*D^{T}*(B^{-1})^{T} ) */
      for (l=0; l<dim; l++){
          for (k=0; k<l+1; k++){
              for (q=0,sum=0; q<no_linpar; q++){
                  sum+=temp[q*dim+l]*temp[q*dim+k];     
              }
              if(k==l) outmat[dim*l-l*(l+1)/2+l] = 1 - sum;
              else outmat[dim*k-k*(k+1)/2+l]= - sum;
          }
      }  
  }
  /* calculating square root of the determinant of DivD */
  if(no_linpar == 1) *det_DivD_half = sqrt(DivD[0]);
  else{
      for(q=0,prod=1;q<no_linpar;q++) prod *= sqDivD[no_linpar*q-q*(q+1)/2+q];  
      *det_DivD_half = prod;
  } 
}


void initz(Real *Y, Real *B, Real *X, Integer dim){
  /* computing X = B^{-1}*Y */
Integer l, k;
Real sum;
  
  for (l=0; l<dim; l++){    
      for (k=0,sum=0; k<l; k++)
          sum += B[dim*k - k*(k+1)/2+l] *X[k];
      X[l] = (Y[l]-sum)/B[dim*l - l*(l+1)/2+l]; 
  }
}

       
void mcmcrun4(Integer *n, Real *data, Real *units, Real *DD, Integer *no_linpar, Integer *cornr, Real *kappa, Real *tausq,
	      Real *distcoords, Real *scale, Real *phiscale, Real *Htrunc, Integer *niter, Integer *subsample, 
	      Integer *burn_in, Real *ss_sigma, Integer *df, Integer *phinr,  Real *phi_discrete, Integer *nmphi, Real *e_mean, Real *SS,
	      Real *phisamples){
  
  typedef Real *Doublearray;  
#define PRN 1000           
  Integer i, j, jj, l, acc, acc_phi, dim2, itr, jstep, jphi ; 
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop, ss4, ss4prop, phi, phiprop, phistep, logp_phi, logp_phiprop ; 
  Real det_DDivDD[2001];
  Doublearray z, zprop, S, Sprop, gradz, gradzprop, temp;
  Doublearray Q, Qprop, AA, Dmat;
  Doublearray BB[2001], DDDMAT[2001];
  Doublearray sqivD, DivD, sqDivD, temp2; 
  
  z = Salloc((*n),Real) ;
  zprop = Salloc((*n),Real) ;
  S = Salloc((*n),Real) ;
  Sprop = Salloc((*n),Real) ;
  gradz = Salloc((*n),Real) ;
  gradzprop = Salloc((*n),Real) ;
  temp = Salloc((*n),Real) ;
  dim2 = ((*n)*(*n)+(*n))/2 ;
  Q = Salloc(dim2,Real) ;
  Qprop = Salloc(dim2,Real) ;
  AA = Salloc(dim2,Real) ;
  Dmat = Salloc(dim2,Real) ;
  sqivD = Salloc((*n)*(*no_linpar),Real) ; 
  DivD = Salloc(((*no_linpar)*(*no_linpar)+(*no_linpar))/2,Real) ;
  sqDivD = Salloc(((*no_linpar)*(*no_linpar)+(*no_linpar))/2,Real); 
  temp2 = Salloc((*no_linpar)*(*n),Real) ;      
  
  if(*nmphi>2001){
    printf(" number of values in phi.discrete must not exceed 2001 \n");
    exit(1);
  }
  for (l=0; l<(*n); l++) S[l]=SS[l];
  for (l=0; l<(*n); l++){
    z[l]=gradz[l]=gradzprop[l]=zprop[l]=Sprop[l]=0;
  }
  for (jj=0; jj<dim2; jj++){
    AA[jj]=Q[jj]=Qprop[jj]=Dmat[jj]=0;   
  }    
  phi = phisamples[0];
  if(*nmphi > 1){
    phistep=phi_discrete[1]-phi_discrete[0];
    jphi=floor((phi-phi_discrete[0])/phistep+0.5);
    if(jphi<0) jphi=0;   
    if(jphi>(*nmphi)-1) jphi=(*nmphi)-1;          
    for (j=0;j<*nmphi;j++){  
      phi=phi_discrete[0]+j*phistep;
      BB[j] = Salloc(dim2,Real) ;
      for (jj=0; jj<dim2; jj++) AA[jj]=corrfct(phi,(*kappa),distcoords[jj],(*cornr));  
      for (l=0; l<(*n); l++) AA[(*n)*l - l*(l+1)/2+l] += (*tausq); 
      cholesky(AA,BB[j],(*n)); 
      DDDMAT[j] = Salloc(dim2,Real) ;
      calc_Dmat(BB[j], DD, DDDMAT[j], &det_DDivDD[j], (*n), (*no_linpar), sqivD, DivD, sqDivD, temp2);
    }  
    phi=phi_discrete[0]+jphi*phistep; 
    Q=BB[jphi];
    Dmat=DDDMAT[jphi];
  }
  else{
    phistep=0;   /* redundant variable in this case */
    jphi=0;
    for (jj=0; jj<dim2; jj++) AA[jj]=corrfct(phi,*kappa,distcoords[jj],(*cornr));
    for (l=0; l<(*n); l++) AA[(*n)*l - l*(l+1)/2+l] += (*tausq); 
    cholesky(AA,Q,(*n)); 
    calc_Dmat(Q, DD, Dmat, &det_DDivDD[0], (*n), (*no_linpar), sqivD, DivD, sqDivD, temp2);  
  }
  initz(S,Q,z,(*n));
  conddensity4(S,Q,&logf,data,z,units,(*n));
  ss4 = calc4_ss(z,Dmat,(*n))+(*ss_sigma);  
  gradient4(S,gradz,Q,Dmat,z,data,units,Htrunc,(*n),ss4,(*df));  
  logp_z = -0.5*(*df)*log(ss4)+log(det_DDivDD[jphi]); 
  logp_phi = logprior_phi(phi, (*e_mean), (*phinr));
  itr = (*niter) + (*burn_in) ;
  acc= 0; acc_phi = 0 ; 
  RANDIN;
  
  for (i=0; i<itr; i++){  
    for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + RNORM*sqrt(*scale);
    conddensity4(Sprop,Q,&logfprop,data,zprop,units,(*n)); 
    ss4prop = calc4_ss(zprop,Dmat,(*n))+(*ss_sigma); 
    gradient4(Sprop,gradzprop,Q,Dmat,zprop,data,units,Htrunc,(*n),ss4prop,(*df)); 
    logp_zprop = -0.5*(*df)*log(ss4prop)+log(det_DDivDD[jphi]); 
    for (logq=0,logqprop=0,l=0; l<(*n); l++){
      logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
      logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
    }
    logq*=(-0.5/(*scale));
    logqprop*=(-0.5/(*scale)); 
    if (log(UNIF)<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
      /*accept */
      logf=logfprop;
      logp_z=logp_zprop;
      ss4 = ss4prop;
      temp=gradz;
      gradz=gradzprop;
      gradzprop=temp;
      temp=S;
      S=Sprop;
      Sprop=temp;
      temp=z;
      z=zprop;
      zprop=temp;
      acc++;
    }       
    if(*nmphi>1) jstep=floor(RNORM*sqrt(*phiscale)/phistep+0.5);
    else jstep=0;
    if ( jphi+jstep > -1 && jphi+jstep < (*nmphi) && !(jstep == 0)){  /* if within limits */  
      phiprop=phi+phistep*jstep;
      logp_phiprop = logprior_phi(phiprop, (*e_mean), (*phinr)); 
      conddensity4(Sprop,BB[jphi+jstep],&logfprop,data,z,units,(*n));                              
      ss4prop = calc4_ss(z,DDDMAT[jphi+jstep],(*n))+(*ss_sigma);                 
      logp_zprop = -0.5*(*df)*log(ss4prop)+log(det_DDivDD[jphi+jstep]);        
      if (log(UNIF)<logfprop-logf+logp_phiprop-logp_phi+logp_zprop-logp_z){    
	phi=phiprop;
	logf=logfprop;
	logp_phi=logp_phiprop;
	logp_z=logp_zprop;
	ss4 = ss4prop;
	jphi+=jstep;
	Q=BB[jphi];
	Dmat=DDDMAT[jphi];   
	temp=S;
	S=Sprop;
	Sprop=temp;
	acc_phi++;
	gradient4(S,gradz,Q,Dmat,z,data,units,Htrunc,(*n),ss4,(*df)); 
      } 
    }         
    if (((i+1-(*burn_in))%(*subsample))==0 && (i+1)>(*burn_in)){
      for (l=0; l<(*n); l++) SS[((i+1-(*burn_in))/(*subsample)-1)*(*n)+l] = S[l]; 
      phisamples[((i+1-(*burn_in))/(*subsample)-1)]=phi; 
    }
    if((i+1)==(*burn_in) && (*burn_in)>0){
      if(*nmphi > 1)
	printf("burn-in =  %d is finished; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",(*burn_in), (Real) acc/(*burn_in), 
	       (Real) acc_phi/(*burn_in));
      else
	printf("burn-in =  %d is finished; Acc.-rate = %1.2f \n",(*burn_in), (Real) acc/(*burn_in));
      acc =0 ; acc_phi =0 ;
    }
    if((i+(*burn_in)+1)%PRN==0 && (i+1)>(*burn_in)){  
      if(*nmphi > 1)
	printf("iter. numb. %d ; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",i+1,(Real) acc/PRN, (Real) acc_phi/PRN);
      else
	printf("iter. numb. %d ; Acc.-rate = %1.2f \n",i+1,(Real) acc/PRN);
      acc =0 ; acc_phi =0 ;
    }   
  }
  RANDOUT;  
}


void conddensity4boxcox(Real *S, Real *B, Real *logfcond, Real *obsdata, Real *z, Real *units, Integer dim, Real lambda){
Integer l, k, ctrl0;
           /* this version requires B to consist of just the entries in the lower triangular part of the matrix */    
   ctrl0 = 0;    
   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=B[dim*k - k*(k+1)/2 + l]*z[k];
      if(S[l] < -1/lambda){
         if(obsdata[l] > 0) ctrl0 = 1 ;
      }
   }
   if(ctrl0 == 1) *logfcond = 0;
   else {
      for (*logfcond=0,l=0; l<dim; l++){
         if(S[l] > -1/lambda) *logfcond += (obsdata[l]*log(lambda*S[l]+1)/lambda-units[l]*pow(lambda*S[l]+1,1/lambda));  
      }
   }
}   


void gradient4boxcox(Real *S, Real *gradz, Real *B, Real *DDmat, Real *z, Real *obsdata, Real *units, Real *Hz, Integer dim, Real ss, 
         Integer df, Real lambda){
   
   Integer l, k;
   Real likeli ;
     /* this version requires B and DDmat to consist of just the entries in the lower triangular part of the matrix */    
   for (l=0; l<dim; l++) gradz[l]=0;
   for (k=0; k<dim; k++){
      likeli=(obsdata[k]- trunc_u(units[k]*pow(lambda*S[k]+1,1/lambda),Hz[k]))*trunc_u(1/(lambda*S[k]+1),pow(Hz[k]/units[k],lambda));
      for (l=0; l<dim; l++){
          if(l<k+1) gradz[l]+=B[dim*l-l*(l+1)/2+k]*likeli - df*DDmat[dim*l-l*(l+1)/2+k]*z[k]/ss ; /* using symmetri of DDmat  */
          else gradz[l] -= df*DDmat[dim*k-k*(k+1)/2+l]*z[k]/ss ; 
      }  
   }
}

       
void mcmcrun4boxcox(Integer *n, Real *data, Real *units, Real *DD, Integer *no_linpar, Integer *cornr, Real *kappa, Real *tausq,
       Real *distcoords, Real *scale, Real *phiscale, Real *Htrunc, Integer *niter, Integer *subsample, 
       Integer *burn_in, Real *ss_sigma, Integer *df, Integer *phinr,  Real *phi_discrete, Integer *nmphi, Real *e_mean, Real *lambda, 
       Real *SS, Real *phisamples){
  
  typedef Real *Doublearray;  
  #define PRN 1000           
  Integer i, j, jj, l, acc, acc_phi, dim2, itr, jstep, jphi , ctrl1, ctrl2, ctrl3; 
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop, ss4, ss4prop, phi, phiprop, phistep, logp_phi, logp_phiprop ; 
  Real det_DDivDD[2001];
  Doublearray z, zprop, S, Sprop, gradz, gradzprop, temp;
  Doublearray Q, Qprop, AA, Dmat;
  Doublearray BB[2001], DDDMAT[2001];
  Doublearray sqivD, DivD, sqDivD, temp2; 
  
  z = Salloc((*n),Real) ;
  zprop = Salloc((*n),Real) ;
  S = Salloc((*n),Real) ;
  Sprop = Salloc((*n),Real) ;
  gradz = Salloc((*n),Real) ;
  gradzprop = Salloc((*n),Real) ;
  temp = Salloc((*n),Real) ;
  dim2 = ((*n)*(*n)+(*n))/2 ;
  Q = Salloc(dim2,Real) ;
  Qprop = Salloc(dim2,Real) ;
  AA = Salloc(dim2,Real) ;
  Dmat = Salloc(dim2,Real) ;
  sqivD = Salloc((*n)*(*no_linpar),Real) ; 
  DivD = Salloc(((*no_linpar)*(*no_linpar)+(*no_linpar))/2,Real) ;
  sqDivD = Salloc(((*no_linpar)*(*no_linpar)+(*no_linpar))/2,Real); 
  temp2 = Salloc((*no_linpar)*(*n),Real) ;      

  if(*nmphi>2001){
       printf(" number of values in phi.discrete must not exceed 2001 \n");
       exit(1);
  }
  for (l=0; l<(*n); l++) S[l]=SS[l];
  for (l=0; l<(*n); l++){
      z[l]=gradz[l]=gradzprop[l]=zprop[l]=Sprop[l]=0;
  }
  for (jj=0; jj<dim2; jj++){
      AA[jj]=Q[jj]=Qprop[jj]=Dmat[jj]=0;   
  }    
  phi = phisamples[0];
  if(*nmphi > 1){
      phistep=phi_discrete[1]-phi_discrete[0];
      jphi=floor((phi-phi_discrete[0])/phistep+0.5);
      if(jphi<0) jphi=0;   
      if(jphi>(*nmphi)-1) jphi=(*nmphi)-1;          
      for (j=0;j<*nmphi;j++){  
         phi=phi_discrete[0]+j*phistep;
         BB[j] = Salloc(dim2,Real) ;
         for (jj=0; jj<dim2; jj++) AA[jj]=corrfct(phi,(*kappa),distcoords[jj],(*cornr));
         for (l=0; l<(*n); l++) AA[(*n)*l - l*(l+1)/2+l] += (*tausq);    
         cholesky(AA,BB[j],(*n)); 
         DDDMAT[j] = Salloc(dim2,Real) ;
         calc_Dmat(BB[j], DD, DDDMAT[j], &det_DDivDD[j], (*n), (*no_linpar), sqivD, DivD, sqDivD, temp2);
      }  
      phi=phi_discrete[0]+jphi*phistep; 
      Q=BB[jphi];
      Dmat=DDDMAT[jphi];
  }
  else{
      phistep=0;   /* redundant variable in this case */
      jphi=0;
      for (jj=0; jj<dim2; jj++) AA[jj]=corrfct(phi,*kappa,distcoords[jj],(*cornr));
      for (l=0; l<(*n); l++) AA[(*n)*l - l*(l+1)/2+l] += (*tausq); 
      cholesky(AA,Q,(*n)); 
      calc_Dmat(Q, DD, Dmat, &det_DDivDD[0], (*n), (*no_linpar), sqivD, DivD, sqDivD, temp2);  
  }  
  initz(S,Q,z,(*n));
  conddensity4boxcox(S,Q,&logf,data,z,units,(*n),(*lambda));
  for(l=0;l<(*n);l++){
     if(S[l] < -1/(*lambda)){
        printf(" Bad starting value for MCMC \n");
        exit(1);
     }
  }
  ss4 = calc4_ss(z,Dmat,(*n))+(*ss_sigma);  
  gradient4boxcox(S,gradz,Q,Dmat,z,data,units,Htrunc,(*n),ss4,(*df),(*lambda));  
  logp_z = -0.5*(*df)*log(ss4)+log(det_DDivDD[jphi]); 
  logp_phi = logprior_phi(phi, (*e_mean), (*phinr));
  itr = (*niter) + (*burn_in) ;
  acc= 0; acc_phi = 0 ; ctrl2 = 0; ctrl3 = 0;
  RANDIN;
 
  for (i=0; i<itr; i++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + RNORM*sqrt(*scale);
        conddensity4boxcox(Sprop,Q,&logfprop,data,zprop,units,(*n),(*lambda));
        ctrl1 = 0;
        for(l=0;l<(*n);l++){
           if(Sprop[l] < -1/(*lambda)){
               if(data[l] > 0) ctrl1 = 1 ;
           }
        }
        if(ctrl1 == 0){
            ss4prop = calc4_ss(zprop,Dmat,(*n))+(*ss_sigma); 
            gradient4boxcox(Sprop,gradzprop,Q,Dmat,zprop,data,units,Htrunc,(*n),ss4prop,(*df),(*lambda)); 
            logp_zprop = -0.5*(*df)*log(ss4prop)+log(det_DDivDD[jphi]); 
            for (logq=0,logqprop=0,l=0; l<(*n); l++){
               logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
               logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
            }
            logq*=(-0.5/(*scale));
            logqprop*=(-0.5/(*scale)); 
            if (log(UNIF)<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	        /*accept */
                logf=logfprop;
                logp_z=logp_zprop;
                ss4 = ss4prop;
                temp=gradz;
                gradz=gradzprop;
                gradzprop=temp;
                temp=S;
                S=Sprop;
                Sprop=temp;
                temp=z;
                z=zprop;
                zprop=temp;
                acc++;
            }
        
        }
        else ctrl2++;      
        if(*nmphi>1) jstep=floor(RNORM*sqrt(*phiscale)/phistep+0.5);
        else jstep=0;
        if ( jphi+jstep > -1 && jphi+jstep < (*nmphi) && !(jstep == 0)){  /* if within limits */  
            phiprop=phi+phistep*jstep;
            logp_phiprop = logprior_phi(phiprop, (*e_mean), (*phinr)); 
            conddensity4boxcox(Sprop,BB[jphi+jstep],&logfprop,data,z,units,(*n),(*lambda));
            ctrl1 = 0;
            for(l=0;l<(*n);l++){
               if(Sprop[l] < -1/(*lambda)){
                   if(data[l] > 0) ctrl1 = 1 ;
               }
            }
            if(ctrl1 == 0){                         
               ss4prop = calc4_ss(z,DDDMAT[jphi+jstep],(*n))+(*ss_sigma);                 
               logp_zprop = -0.5*(*df)*log(ss4prop)+log(det_DDivDD[jphi+jstep]);        
               if (ctrl1 == 0 && log(UNIF)<logfprop-logf+logp_phiprop-logp_phi+logp_zprop-logp_z){    
                  phi=phiprop;
                  logf=logfprop;
                  logp_phi=logp_phiprop;
                  logp_z=logp_zprop;
                  ss4 = ss4prop;
                  jphi+=jstep;
                  Q=BB[jphi];
                  Dmat=DDDMAT[jphi];   
                  temp=S;
                  S=Sprop;
                  Sprop=temp;
                  acc_phi++;
                  gradient4boxcox(S,gradz,Q,Dmat,z,data,units,Htrunc,(*n),ss4,(*df),(*lambda)); 
               }
            }
            else ctrl3++;  
        }         
        if (((i+1-(*burn_in))%(*subsample))==0 && (i+1)>(*burn_in)){
              for (l=0; l<(*n); l++) SS[((i+1-(*burn_in))/(*subsample)-1)*(*n)+l] = S[l]; 
              phisamples[((i+1-(*burn_in))/(*subsample)-1)]=phi; 
        }
        if((i+1)==(*burn_in) && (*burn_in)>0){
	  if(*nmphi > 1)
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",(*burn_in), (Real) acc/(*burn_in), 
		   (Real) acc_phi/(*burn_in));
	  else
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f \n",(*burn_in), (Real) acc/(*burn_in));
	  acc =0 ; acc_phi =0 ;
        }
        if((i+(*burn_in)+1)%PRN==0 && (i+1)>(*burn_in)){  
	  if(*nmphi > 1)
	    printf("iter. numb. %d ; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",i+1,(Real) acc/PRN, (Real) acc_phi/PRN);
	  else
	    printf("iter. numb. %d ; Acc.-rate = %1.2f \n",i+1,(Real) acc/PRN);
	  acc =0 ; acc_phi =0 ;
        }   
  }
  if(ctrl2*10 > itr || ctrl3*10 > itr){
    printf("rejection of proposals cauced by density for proposal being zero: S: %d; phi: %d , out of %d iterations \n",ctrl2,ctrl3,itr);
  }
  RANDOUT;  
}


void conddensityboxcox(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *linpred, Real *units, Integer dim, Real lambda){
Integer l, k, ctrl0 ;

   ctrl0 = 0;    
   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[l*dim+k]*z[k];
      if(S[l]+linpred[l] < -1/lambda){
         if(obsdata[l] > 0) ctrl0 = 1 ;
      }
   }
   if(ctrl0 == 1) *logfcond = 0;
   else {
      for (*logfcond=0,l=0; l<dim; l++){
         if(S[l]+linpred[l] > -1/lambda) *logfcond+= (obsdata[l]*log(lambda*(S[l]+linpred[l])+1)/lambda-units[l]*pow(lambda*(S[l]+linpred[l])+1,1/lambda));
      }
   }
}


void gradientboxcox(Real *S, Real *gradz, Real *BB, Real *z, Real *obsdata, Real *linpred, Real *units, Real *Hz, Integer dim, Real lambda){
   Integer l, k;
   Real likeli;
   
   for (l=0; l<dim; l++) gradz[l]=-z[l];
   for (k=0; k<dim; k++){
      likeli=(obsdata[k]- trunc_u(units[k]*pow(lambda*(S[k]+linpred[k])+1,1/lambda),Hz[k]))*trunc_u(1/(lambda*(S[k]+linpred[k])+1),pow(Hz[k]/units[k],lambda));
      for (l=0; l<k+1; l++) gradz[l]+=BB[k*dim+l]*likeli;  
   }
}


void mcmcrunboxcox(Integer *n, Real *zz, Real *SS, Real *data, Real *units, Real *meanS, Real *QQ, 
         Real *randnormal, Real *randunif, Real *Htrunc, Real *scale, Integer *nsim, Integer *subsample, Real *lambda, Real *acc_rate){
         
  Integer i, ii, l, acc, ctrl1, ctrl2 ;  
  Real *z = Salloc((*n),Real) ;
  Real *zprop = Salloc((*n),Real) ;
  Real *S = Salloc((*n),Real) ; 
  Real *Sprop = Salloc((*n),Real) ;
  Real *gradz = Salloc((*n),Real) ;
  Real *gradzprop = Salloc((*n),Real) ;
  Real *temp = Salloc((*n),Real) ;
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop  ; 
  
  for (l=0; l<(*n); l++){
      z[l] = zz[l]; 
      S[l] = 0;     
  } 
  conddensityboxcox(S,QQ,&logf,data,z,meanS,units,(*n),(*lambda));
  for(l=0;l<(*n);l++){
       if(S[l]+meanS[l] < -1/(*lambda)){
           printf(" Bad starting value for MCMC \n");
           exit(1);
       }
  }           
  gradientboxcox(S,gradz,QQ,z,data,meanS,units,Htrunc,(*n),(*lambda));
  logp_z = -calc_ss(z,(*n))/2;
  acc= 0; ctrl2 = 0;     
  for(i=0; i<(*nsim); i++){ 
     for(ii=0; ii<(*subsample); ii++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + randnormal[(i*(*subsample)+ii)*(*n)+l];
        conddensityboxcox(Sprop,QQ,&logfprop,data,zprop,meanS,units,(*n),(*lambda));
        ctrl1 = 0;
        for(l=0;l<(*n);l++){
           if(Sprop[l]+meanS[l] < -1/(*lambda)){
               if(data[l] > 0) ctrl1 = 1 ;
           }
        }
        if(ctrl1 == 0){
           gradientboxcox(Sprop,gradzprop,QQ,zprop,data,meanS,units,Htrunc,(*n),(*lambda));
           logp_zprop=-calc_ss(zprop,(*n))/2;  
           for (logq=0,logqprop=0,l=0; l<(*n); l++){
              logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
              logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
           }
           logq*=(-0.5/(*scale));
           logqprop*=(-0.5/(*scale));
           if (log(randunif[i*(*subsample)+ii])<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	      /*accept */
               logf=logfprop;
               logp_z=logp_zprop;
               temp=gradz;
               gradz=gradzprop;
               gradzprop=temp;
               temp=S;
               S=Sprop;
               Sprop=temp;
               temp=z;
               z=zprop;
               zprop=temp;
               acc++; 
           } 
        }
        else ctrl2++;         
     } 
     for (l=0; l<(*n); l++) SS[i*(*n)+l] = S[l];  
  }
  if(ctrl2*10 > (*nsim)*(*subsample)){
     printf(" rejection of proposals for S cauced by density for proposal being zero: %d out of %d iterations \n",ctrl2,(*nsim)*(*subsample));
  }
  *acc_rate = (Real) acc/((*nsim)*(*subsample));
  for (l=0; l<(*n); l++) zz[l] = z[l]; 
}

void conddensity2boxcox(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *units, Integer dim, Real lambda){
Integer l, k, ctrl0 ;

   ctrl0 = 0;    
   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[l*dim+k]*z[k];
      if(S[l] < -1/lambda){
         if(obsdata[l] > 0) ctrl0 = 1 ;
      }
   }
   if(ctrl0 == 1) *logfcond = 0;
   else {
      for (*logfcond=0,l=0; l<dim; l++){
         if(S[l] > -1/lambda) *logfcond+= (obsdata[l]*log(lambda*S[l]+1)/lambda-units[l]*pow(lambda*S[l]+1,1/lambda));
      }
   }
}

void gradient2boxcox(Real *S, Real *gradz, Real *BB, Real *Dmat, Real *z, Real *obsdata, Real *units, Real *Hz, Integer dim, Real lambda){
   Integer l, k;
   Real likeli ;
   
   for (l=0; l<dim; l++) gradz[l]=0;
   for (k=0; k<dim; k++){
      likeli=(obsdata[k]- trunc_u(units[k]*pow(lambda*S[k]+1,1/lambda),Hz[k]))*trunc_u(1/(lambda*S[k]+1),pow(Hz[k]/units[k],lambda));
      for (l=0; l<dim; l++)
          if(l<k+1) gradz[l]+=BB[k*dim+l]*likeli - Dmat[l*dim+k]*z[k] ;
          else gradz[l]-=Dmat[k*dim+l]*z[k] ;
   }
}

void mcmcrun2boxcox(Integer *n, Real *zz, Real *SS, Real *data, Real *units, Real *QQ, Real *Dmat,
         Real *randnormal, Real *randunif, Real *Htrunc, Real *scale, Integer *nsim, Integer *subsample, Real *lambda, Real *acc_rate){
         
  Integer i, ii, l, acc, ctrl1, ctrl2 ;  
  Real *z = Salloc((*n),Real) ; 
  Real *zprop = Salloc((*n),Real) ; 
  Real *S = Salloc((*n),Real) ;
  Real *Sprop = Salloc((*n),Real) ;
  Real *gradz = Salloc((*n),Real) ;
  Real *gradzprop = Salloc((*n),Real) ;
  Real *temp = Salloc((*n),Real) ;
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop  ;   
  for (l=0; l<(*n); l++){
      z[l] = zz[l];     
      S[l] = 0; 
  } 
  conddensity2boxcox(S,QQ,&logf,data,z,units,(*n),(*lambda)); 
  for(l=0;l<(*n);l++){
       if(S[l] < -1/(*lambda)){
            printf(" Bad starting value for MCMC \n");
            exit(1);
       }
  }           
  gradient2boxcox(S,gradz,QQ,Dmat,z,data,units,Htrunc,(*n),(*lambda));
  logp_z = -calc2_ss(z,Dmat,(*n))/2;
  acc= 0; ctrl2 = 0;
  for(i=0; i<(*nsim); i++){ 
     for(ii=0; ii<(*subsample); ii++){ 
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + randnormal[(i*(*subsample)+ii)*(*n)+l];
        conddensity2boxcox(Sprop,QQ,&logfprop,data,zprop,units,(*n),(*lambda));
        ctrl1 = 0;
        for(l=0;l<(*n);l++){
           if(Sprop[l] < -1/(*lambda)){
               if(data[l] > 0) ctrl1 = 1 ;
           }
        }
        if(ctrl1 == 0){  
           gradient2boxcox(Sprop,gradzprop,QQ,Dmat,zprop,data,units,Htrunc,(*n),(*lambda));
           logp_zprop=-calc2_ss(zprop,Dmat,(*n))/2;  
           for (logq=0,logqprop=0,l=0; l<(*n); l++){
              logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
              logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
           }
           logq*=(-0.5/(*scale));
           logqprop*=(-0.5/(*scale));     
           if (log(randunif[i*(*subsample)+ii])<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	      /*accept */
               logf=logfprop;
               logp_z=logp_zprop;
               temp=gradz;
               gradz=gradzprop;
               gradzprop=temp;
               temp=S;
               S=Sprop;
               Sprop=temp;
               temp=z;
               z=zprop;
               zprop=temp;
               acc++;
           } 
        }
        else ctrl2++;   
     }     
     for (l=0; l<(*n); l++) SS[i*(*n)+l] = S[l];     
  }
  if(ctrl2*10 > (*nsim)*(*subsample)){
     printf(" rejection of proposals for S cauced by density for proposal being zero: %d out of %d iterations \n",ctrl2,(*nsim)*(*subsample));
  }
  *acc_rate = (Real) acc/((*nsim)*(*subsample));  
  for (l=0; l<(*n); l++) zz[l] = z[l]; 
}


void conddensitybinom(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *linpred, Real *units, Integer dim){
Integer l, k;

  for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[l*dim+k]*z[k];
   }  
   for (*logfcond=0,l=0; l<dim; l++){
      *logfcond+= (obsdata[l]*(S[l]+linpred[l])-units[l]*log(1+exp(S[l]+linpred[l])));
   }
}


void gradientbinom(Real *S, Real *gradz, Real *BB, Real *z, Real *obsdata, Real *linpred, Real *units, Integer dim){
   Integer l, k;
   Real likeli;
   
   for (l=0; l<dim; l++) gradz[l]=-z[l];
   for (k=0; k<dim; k++){
      likeli = obsdata[k] - units[k]*exp(S[k]+linpred[k])/(1+exp(S[k]+linpred[k]));    
      for (l=0; l<k+1; l++) gradz[l]+=BB[k*dim+l]*likeli;  
   }
}

void mcmcrunbinom(Integer *n, Real *zz, Real *SS, Real *data, 
                  Real *units, Real *meanS, Real *QQ, 
                  Real *randnormal, Real *randunif,  Real *scale, 
                  Integer *nsim, Integer *subsample, Real *acc_rate){
         
  Integer i, ii, l, acc ;  
  Real *z = Salloc((*n),Real) ;
  Real *zprop = Salloc((*n),Real) ;
  Real *S = Salloc((*n),Real) ; 
  Real *Sprop = Salloc((*n),Real) ;
  Real *gradz = Salloc((*n),Real) ;
  Real *gradzprop = Salloc((*n),Real) ;
  Real *temp = Salloc((*n),Real) ;
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop  ; 
  
  for (l=0; l<(*n); l++){
      z[l] = zz[l]; 
      S[l] = 0;     
  } 
  conddensitybinom(S,QQ,&logf,data,z,meanS,units,(*n));           
  gradientbinom(S,gradz,QQ,z,data,meanS,units,(*n));
  logp_z = -calc_ss(z,(*n))/2;
  acc= 0;       
  for(i=0; i<(*nsim); i++){ 
     for(ii=0; ii<(*subsample); ii++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + randnormal[(i*(*subsample)+ii)*(*n)+l];
        conddensitybinom(Sprop,QQ,&logfprop,data,zprop,meanS,units,(*n));     
        gradientbinom(Sprop,gradzprop,QQ,zprop,data,meanS,units,(*n));
        logp_zprop=-calc_ss(zprop,(*n))/2;  
        for (logq=0,logqprop=0,l=0; l<(*n); l++){
           logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
           logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
        }
        logq*=(-0.5/(*scale));
        logqprop*=(-0.5/(*scale));
        if (log(randunif[i*(*subsample)+ii])<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	   /*accept */
            logf=logfprop;
            logp_z=logp_zprop;
            temp=gradz;
            gradz=gradzprop;
            gradzprop=temp;
            temp=S;
            S=Sprop;
            Sprop=temp;
            temp=z;
            z=zprop;
            zprop=temp;
            acc++; 
        }      
     } 
     for (l=0; l<(*n); l++) SS[i*(*n)+l] = S[l];  
  }
  *acc_rate = (Real) acc/((*nsim)*(*subsample));
  for (l=0; l<(*n); l++) zz[l] = z[l]; 
}

void conddensity2binom(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *units, Integer dim){
Integer l, k ;

   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[l*dim+k]*z[k];
   }
   for (*logfcond=0,l=0; l<dim; l++){
      *logfcond+= (obsdata[l]*S[l]-units[l]*log(1+exp(S[l])));
   }
}

void gradient2binom(Real *S, Real *gradz, Real *BB, Real *Dmat, Real *z, Real *obsdata, Real *units, Integer dim){
   Integer l, k;
   Real likeli ;
   
   for (l=0; l<dim; l++) gradz[l]=0;
   for (k=0; k<dim; k++){
      likeli = obsdata[k] - units[k]*exp(S[k])/(1+exp(S[k]));
      for (l=0; l<dim; l++)
          if(l<k+1) gradz[l]+=BB[k*dim+l]*likeli - Dmat[l*dim+k]*z[k] ;
          else gradz[l]-=Dmat[k*dim+l]*z[k] ;
   }
}

void mcmcrun2binom(Integer *n, Real *zz, Real *SS, Real *data, Real *units, Real *QQ, Real *Dmat,
         Real *randnormal, Real *randunif, Real *scale, Integer *nsim, Integer *subsample, Real *acc_rate){
         
  Integer i, ii, l, acc ;  
  Real *z = Salloc((*n),Real) ; 
  Real *zprop = Salloc((*n),Real) ; 
  Real *S = Salloc((*n),Real) ;
  Real *Sprop = Salloc((*n),Real) ;
  Real *gradz = Salloc((*n),Real) ;
  Real *gradzprop = Salloc((*n),Real) ;
  Real *temp = Salloc((*n),Real) ;
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop  ;   
  for (l=0; l<(*n); l++){
      z[l] = zz[l];     
      S[l] = 0; 
  } 
  conddensity2binom(S,QQ,&logf,data,z,units,(*n));           
  gradient2binom(S,gradz,QQ,Dmat,z,data,units,(*n));
  logp_z = -calc2_ss(z,Dmat,(*n))/2;
  acc= 0; 
  for(i=0; i<(*nsim); i++){ 
     for(ii=0; ii<(*subsample); ii++){ 
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + randnormal[(i*(*subsample)+ii)*(*n)+l];
        conddensity2binom(Sprop,QQ,&logfprop,data,zprop,units,(*n));  
        gradient2binom(Sprop,gradzprop,QQ,Dmat,zprop,data,units,(*n));
        logp_zprop=-calc2_ss(zprop,Dmat,(*n))/2;  
        for (logq=0,logqprop=0,l=0; l<(*n); l++){
           logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
           logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
        }
        logq*=(-0.5/(*scale));
        logqprop*=(-0.5/(*scale));      
        if (log(randunif[i*(*subsample)+ii])<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
            /*accept */
            logf=logfprop;
            logp_z=logp_zprop;
            temp=gradz;
            gradz=gradzprop;
            gradzprop=temp;
            temp=S;
            S=Sprop;
            Sprop=temp;
            temp=z;
            z=zprop;
            zprop=temp;
            acc++;
        }   
     }     
     for (l=0; l<(*n); l++) SS[i*(*n)+l] = S[l];     
  }
  *acc_rate = (Real) acc/((*nsim)*(*subsample));  
  for (l=0; l<(*n); l++) zz[l] = z[l]; 
}


void conddensity4binom(Real *S, Real *B, Real *logfcond, Real *obsdata, Real *z, Real *units, Integer dim){
Integer l, k;
           /* this version requires B to consist of just the entries in the lower triangular part of the matrix */          
   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=B[dim*k - k*(k+1)/2 + l]*z[k];
   }
   for (*logfcond=0,l=0; l<dim; l++){
     *logfcond+=(obsdata[l]*(S[l])-units[l]*log(1+exp(S[l])));
   }
} 

void gradient4binom(Real *S, Real *gradz, Real *B, Real *DDmat, Real *z, Real *obsdata, Real *units, Integer dim, Real ss, Integer df){
   Integer l, k;
   Real likeli ;
     /* this version requires B and DDmat to consist of just the entries in the lower triangular part of the matrix */    
   for (l=0; l<dim; l++) gradz[l]=0;
   for (k=0; k<dim; k++){
      likeli = obsdata[k] - units[k]*exp(S[k])/(1+exp(S[k]));
      for (l=0; l<dim; l++){
          if(l<k+1) gradz[l]+=B[dim*l-l*(l+1)/2+k]*likeli - df*DDmat[dim*l-l*(l+1)/2+k]*z[k]/ss ; /* using symmetri of DDmat  */
          else gradz[l] -= df*DDmat[dim*k-k*(k+1)/2+l]*z[k]/ss ;   
      }  
   }
}
       
void mcmcrun4binom(Integer *n, Real *data, Real *units, Real *DD, Integer *no_linpar, Integer *cornr, Real *kappa, Real *tausq,
       Real *distcoords, Real *scale, Real *phiscale, Integer *niter, Integer *subsample, 
       Integer *burn_in, Real *ss_sigma, Integer *df, Integer *phinr,  Real *phi_discrete, Integer *nmphi, Real *e_mean, Real *SS, 
       Real *phisamples){
  
  typedef Real *Doublearray;  
  #define PRN 1000           
  Integer i, j, jj, l, acc, acc_phi, dim2, itr, jstep, jphi ; 
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop, ss4, ss4prop, phi, phiprop, phistep, logp_phi, logp_phiprop ; 
  Real det_DDivDD[2001];
  Doublearray z, zprop, S, Sprop, gradz, gradzprop, temp;
  Doublearray Q, Qprop, AA, Dmat;
  Doublearray BB[2001], DDDMAT[2001];
  Doublearray sqivD, DivD, sqDivD, temp2; 
  
  z = Salloc((*n),Real) ;
  zprop = Salloc((*n),Real) ;
  S = Salloc((*n),Real) ;
  Sprop = Salloc((*n),Real) ;
  gradz = Salloc((*n),Real) ;
  gradzprop = Salloc((*n),Real) ;
  temp = Salloc((*n),Real) ;
  dim2 = ((*n)*(*n)+(*n))/2 ;
  Q = Salloc(dim2,Real) ;
  Qprop = Salloc(dim2,Real) ;
  AA = Salloc(dim2,Real) ;
  Dmat = Salloc(dim2,Real) ;
  sqivD = Salloc((*n)*(*no_linpar),Real) ; 
  DivD = Salloc(((*no_linpar)*(*no_linpar)+(*no_linpar))/2,Real) ;
  sqDivD = Salloc(((*no_linpar)*(*no_linpar)+(*no_linpar))/2,Real); 
  temp2 = Salloc((*no_linpar)*(*n),Real) ;      

  if(*nmphi>2001){
       printf(" number of values in phi.discrete must not exceed 2001 \n");
       exit(1);
  }
  for (l=0; l<(*n); l++) S[l]=SS[l];
  for (l=0; l<(*n); l++){
      z[l]=gradz[l]=gradzprop[l]=zprop[l]=Sprop[l]=0;
  }
  for (jj=0; jj<dim2; jj++){
      AA[jj]=Q[jj]=Qprop[jj]=Dmat[jj]=0;   
  }    
  phi = phisamples[0];
  if(*nmphi > 1){
      phistep=phi_discrete[1]-phi_discrete[0];
      jphi=floor((phi-phi_discrete[0])/phistep+0.5);
      if(jphi<0) jphi=0;   
      if(jphi>(*nmphi)-1) jphi=(*nmphi)-1;          
      for (j=0;j<(*nmphi);j++){  
         phi=phi_discrete[0]+j*phistep;
         BB[j] = Salloc(dim2,Real) ;
         for (jj=0; jj<dim2; jj++) AA[jj]=corrfct(phi,(*kappa),distcoords[jj],(*cornr));  
         for (l=0; l<(*n); l++) AA[(*n)*l - l*(l+1)/2+l] += (*tausq); 
         cholesky(AA,BB[j],(*n)); 
         DDDMAT[j] = Salloc(dim2,Real) ;
         calc_Dmat(BB[j], DD, DDDMAT[j], &det_DDivDD[j], (*n), (*no_linpar), sqivD, DivD, sqDivD, temp2);
      }  
      phi=phi_discrete[0]+jphi*phistep; 
      Q=BB[jphi];
      Dmat=DDDMAT[jphi];
  }
  else{
      phistep=0;   /* redundant variable in this case */
      jphi=0;
      for (jj=0; jj<dim2; jj++) AA[jj]=corrfct(phi,*kappa,distcoords[jj],(*cornr));
      for (l=0; l<(*n); l++) AA[(*n)*l - l*(l+1)/2+l] += (*tausq); 
      cholesky(AA,Q,(*n)); 
      calc_Dmat(Q, DD, Dmat, &det_DDivDD[0], (*n), (*no_linpar), sqivD, DivD, sqDivD, temp2);  
  }  
  initz(S,Q,z,(*n));
  conddensity4binom(S,Q,&logf,data,z,units,(*n));
  ss4 = calc4_ss(z,Dmat,(*n))+(*ss_sigma);  
  gradient4binom(S,gradz,Q,Dmat,z,data,units,(*n),ss4,(*df));  
  logp_z = -0.5*(*df)*log(ss4)+log(det_DDivDD[jphi]); 
  logp_phi = logprior_phi(phi, (*e_mean), (*phinr));
  itr = (*niter) + (*burn_in) ;
  acc= 0; acc_phi = 0 ; 
  RANDIN;
 
  for (i=0; i<itr; i++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + RNORM*sqrt(*scale);
        conddensity4binom(Sprop,Q,&logfprop,data,zprop,units,(*n)); 
        ss4prop = calc4_ss(zprop,Dmat,(*n))+(*ss_sigma); 
        gradient4binom(Sprop,gradzprop,Q,Dmat,zprop,data,units,(*n),ss4prop,(*df)); 
        logp_zprop = -0.5*(*df)*log(ss4prop)+log(det_DDivDD[jphi]); 
        for (logq=0,logqprop=0,l=0; l<(*n); l++){
           logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
           logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
        }
        logq*=(-0.5/(*scale));
        logqprop*=(-0.5/(*scale)); 
        if (log(UNIF)<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	    /*accept */
            logf=logfprop;
            logp_z=logp_zprop;
            ss4 = ss4prop;
            temp=gradz;
            gradz=gradzprop;
            gradzprop=temp;
            temp=S;
            S=Sprop;
            Sprop=temp;
            temp=z;
            z=zprop;
            zprop=temp;
            acc++;
        }       
        if(*nmphi>1) jstep=floor(RNORM*sqrt(*phiscale)/phistep+0.5);
        else jstep=0;
        if ( jphi+jstep > -1 && jphi+jstep < (*nmphi) && !(jstep == 0)){  /* if within limits */  
            phiprop=phi+phistep*jstep;
            logp_phiprop = logprior_phi(phiprop, (*e_mean), (*phinr)); 
            conddensity4binom(Sprop,BB[jphi+jstep],&logfprop,data,z,units,(*n));                              
            ss4prop = calc4_ss(z,DDDMAT[jphi+jstep],(*n))+(*ss_sigma);                 
            logp_zprop = -0.5*(*df)*log(ss4prop)+log(det_DDivDD[jphi+jstep]);        
            if (log(UNIF)<logfprop-logf+logp_phiprop-logp_phi+logp_zprop-logp_z){    
               phi=phiprop;
               logf=logfprop;
               logp_phi=logp_phiprop;
               logp_z=logp_zprop;
               ss4 = ss4prop;
               jphi+=jstep;
               Q=BB[jphi];
               Dmat=DDDMAT[jphi];   
               temp=S;
               S=Sprop;
               Sprop=temp;
               acc_phi++;
               gradient4binom(S,gradz,Q,Dmat,z,data,units,(*n),ss4,(*df)); 
            } 
        }         
        if (((i+1-(*burn_in))%(*subsample))==0 && (i+1)>(*burn_in)){
              for (l=0; l<(*n); l++) SS[((i+1-(*burn_in))/(*subsample)-1)*(*n)+l] = S[l]; 
              phisamples[((i+1-(*burn_in))/(*subsample)-1)]=phi; 
        }
	if((i+1)==(*burn_in) && (*burn_in)>0){
	  if(*nmphi > 1)
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",(*burn_in), (Real) acc/(*burn_in), 
		   (Real) acc_phi/(*burn_in));
	  else
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f \n",(*burn_in), (Real) acc/(*burn_in));
	  acc =0 ; acc_phi =0 ;
	}
	if((i+(*burn_in)+1)%PRN==0 && (i+1)>(*burn_in)){  
	  if(*nmphi > 1)
	    printf("iter. numb. %d ; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",i+1,(Real) acc/PRN, (Real) acc_phi/PRN);
	  else
	    printf("iter. numb. %d ; Acc.-rate = %1.2f \n",i+1,(Real) acc/PRN);
	  acc =0 ; acc_phi =0 ;
	}   
  }
  RANDOUT;  
}

void conddensity5binom(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *linpred, Real *units, Integer dim){
Integer l, k;

  for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[dim*k - k*(k+1)/2 + l]*z[k];
   }  
   for (*logfcond=0,l=0; l<dim; l++){
      *logfcond+= (obsdata[l]*(S[l]+linpred[l])-units[l]*log(1+exp(S[l]+linpred[l])));
   }
}

void gradient5binom(Real *S, Real *gradz, Real *BB, Real *z, Real *obsdata, Real *linpred, Real *units, Integer dim, Real ss, Integer df){
   Integer l, k;
   Real likeli;
   
   for (l=0; l<dim; l++) gradz[l]=-z[l]*df/ss;
   for (k=0; k<dim; k++){
      likeli = obsdata[k] - units[k]*exp(S[k]+linpred[k])/(1+exp(S[k]+linpred[k]));    
      for (l=0; l<k+1; l++) gradz[l]+=BB[dim*l-l*(l+1)/2+k]*likeli;  
   }
}
       
void mcmcrun5binom(Integer *n, Real *data, Real *units, Real *meanS, Real *DDvbetaDD, Integer *cornr, Real *kappa, Real *tausq,
       Real *distcoords, Real *scale, Real *phiscale, Integer *niter, Integer *subsample, 
       Integer *burn_in, Real *ss_sigma, Integer *df, Integer *phinr,  Real *phi_discrete, Integer *nmphi, Real *e_mean, Real *SS, 
       Real *phisamples){
  
  typedef Real *Doublearray;  
  #define PRN 1000           
  Integer i, j, jj, l, k, acc, acc_phi, dim2, itr, jstep, jphi; 
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop, ss, ssprop, phi, phiprop, phistep, logp_phi, logp_phiprop ; 
  Doublearray z, zprop, S, Sprop, gradz, gradzprop, temp;
  Doublearray Q, Qprop, AA;
  Doublearray BB[2001];

  
  z = Salloc((*n),Real) ;
  zprop = Salloc((*n),Real) ;
  S = Salloc((*n),Real) ;
  Sprop = Salloc((*n),Real) ;
  gradz = Salloc((*n),Real) ;
  gradzprop = Salloc((*n),Real) ;
  temp = Salloc((*n),Real) ;
  dim2 = ((*n)*(*n)+(*n))/2 ;
  Q = Salloc(dim2,Real) ;
  Qprop = Salloc(dim2,Real) ;
  AA = Salloc(dim2,Real) ;

  if(*nmphi>2001){
       printf(" number of values in phi.discrete must not exceed 2001 \n");
       exit(1);
  }
  for (l=0; l<(*n); l++) S[l]=SS[l];
  for (l=0; l<(*n); l++){
      z[l]=gradz[l]=gradzprop[l]=zprop[l]=Sprop[l]=0;
  }
  for (jj=0; jj<dim2; jj++){
      AA[jj]=Q[jj]=Qprop[jj]=0;   
  }    
  phi = phisamples[0];
  if(*nmphi > 1){
      phistep=phi_discrete[1]-phi_discrete[0];
      jphi=floor((phi-phi_discrete[0])/phistep+0.5);
      if(jphi<0) jphi=0;   
      if(jphi>(*nmphi)-1) jphi=(*nmphi)-1;          
      for (j=0;j<(*nmphi);j++){  
         phi=phi_discrete[0]+j*phistep;
         BB[j] = Salloc(dim2,Real) ;
         for (l=0; l<(*n); l++){
	    for (k=0; k<l; k++){
               AA[(*n)*k - k*(k+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*k - k*(k+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+k];
            }
            AA[(*n)*l - l*(l+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*l - l*(l+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+l] + (*tausq);  
	 }

         cholesky(AA,BB[j],(*n));
      }
      phi=phi_discrete[0]+jphi*phistep; 
      Q=BB[jphi];
  }
  else{
      phistep=0;   /* redundant variable in this case */
      jphi=0;
      for (l=0; l<(*n); l++){
	 for (k=0; k<l; k++){
            AA[(*n)*k - k*(k+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*k - k*(k+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+k];
         }
         AA[(*n)*l - l*(l+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*l - l*(l+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+l] + (*tausq); 
      }
      cholesky(AA,Q,(*n)); 
  }
  initz(S,Q,z,(*n));
  conddensity5binom(S,Q,&logf,data,z,meanS,units,(*n));    
  ss = calc_ss(z,(*n))+(*ss_sigma);  
  gradient5binom(S,gradz,Q,z,data,meanS,units,(*n),ss,(*df));  
  logp_z = -0.5*(*df)*log(ss); 
  logp_phi = logprior_phi(phi, (*e_mean), (*phinr));
  itr = (*niter) + (*burn_in) ;
  acc= 0; acc_phi = 0 ; 
  RANDIN;
 
  for (i=0; i<itr; i++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + RNORM*sqrt(*scale); 
        conddensity5binom(Sprop,Q,&logfprop,data,zprop,meanS,units,(*n));    
        ssprop = calc_ss(zprop,(*n))+(*ss_sigma);  
        gradient5binom(Sprop,gradzprop,Q,zprop,data,meanS,units,(*n),ssprop,(*df));  
        logp_zprop = -0.5*(*df)*log(ssprop); 
        for (logq=0,logqprop=0,l=0; l<(*n); l++){
           logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
           logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
        }
        logq*=(-0.5/(*scale));
        logqprop*=(-0.5/(*scale)); 
        if (log(UNIF)<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	    /*accept */
            logf=logfprop;
            logp_z=logp_zprop;
            ss = ssprop;
            temp=gradz;
            gradz=gradzprop;
            gradzprop=temp;
            temp=S;
            S=Sprop;
            Sprop=temp;
            temp=z;
            z=zprop;
            zprop=temp;
            acc++;
        }       
        if(*nmphi>1) jstep=floor(RNORM*sqrt(*phiscale)/phistep+0.5);
        else jstep=0;
        if ( jphi+jstep > -1 && jphi+jstep < (*nmphi) && !(jstep == 0)){  /* if within limits */  
            phiprop=phi+phistep*jstep;
            logp_phiprop = logprior_phi(phiprop, (*e_mean), (*phinr));
            conddensity5binom(Sprop,BB[jphi+jstep],&logfprop,data,z,meanS,units,(*n));
            if (log(UNIF)<logfprop-logf+logp_phiprop-logp_phi){    
               phi=phiprop;
               logf=logfprop;
               logp_phi=logp_phiprop;
               jphi+=jstep;
               Q=BB[jphi];   
               temp=S;
               S=Sprop;
               Sprop=temp;
               acc_phi++;
               gradient5binom(S,gradz,Q,z,data,meanS,units,(*n),ssprop,(*df));  
            } 
        }         
        if (((i+1-(*burn_in))%(*subsample))==0 && (i+1)>(*burn_in)){
              for (l=0; l<(*n); l++) SS[((i+1-(*burn_in))/(*subsample)-1)*(*n)+l] = S[l]; 
              phisamples[((i+1-(*burn_in))/(*subsample)-1)]=phi; 
        }
	if((i+1)==(*burn_in) && (*burn_in)>0){
	  if(*nmphi > 1)
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",(*burn_in), (Real) acc/(*burn_in), 
		   (Real) acc_phi/(*burn_in));
	  else
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f \n",(*burn_in), (Real) acc/(*burn_in));
	  acc =0 ; acc_phi =0 ;
	}
	if((i+(*burn_in)+1)%PRN==0 && (i+1)>(*burn_in)){  
	  if(*nmphi > 1)
	    printf("iter. numb. %d ; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",i+1,(Real) acc/PRN, (Real) acc_phi/PRN);
	  else
	    printf("iter. numb. %d ; Acc.-rate = %1.2f \n",i+1,(Real) acc/PRN);
	  acc =0 ; acc_phi =0 ;
	}   
  }
  RANDOUT;  
}


void conddensity5boxcox(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *linpred, Real *units, Integer dim, Real lambda){
Integer l, k, ctrl0;
  
   ctrl0 = 0;    
   for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[dim*k - k*(k+1)/2 + l]*z[k]; 
      if(S[l] < -1/lambda){
         if(obsdata[l] > 0) ctrl0 = 1 ;
      }
   } 
   if(ctrl0 == 1) *logfcond = 0;
   else { 
     for (*logfcond=0,l=0; l<dim; l++){
       if((S[l]+linpred[l]) > -1/lambda) 
	 *logfcond+= (obsdata[l]*log(lambda*(S[l]+linpred[l])+1)/lambda-units[l]*pow(lambda*(S[l]+linpred[l])+1,1/lambda));  
     }
   }  
}


void gradient5boxcox(Real *S, Real *gradz, Real *BB, Real *z, Real *obsdata, Real *linpred, Real *units, Real *Hz, Integer dim, Real ss, 
   Integer df, Real lambda){

   Integer l, k;
   Real likeli;

   for (l=0; l<dim; l++) gradz[l]=-z[l]*df/ss;
   for (k=0; k<dim; k++){
     likeli=(obsdata[k]-trunc_u(units[k]*pow(lambda*(S[k]+linpred[k])+1,1/lambda),Hz[k]))*trunc_u(1/(lambda*(S[k]+linpred[k])+1),pow(Hz[k]/units[k],lambda));
     for (l=0; l<k+1; l++) gradz[l]+=BB[dim*l-l*(l+1)/2+k]*likeli ;
   }
}


void mcmcrun5boxcox(Integer *n, Real *data, Real *units, Real *meanS, Real *DDvbetaDD, Integer *cornr, Real *kappa, Real *tausq,
       Real *distcoords, Real *scale, Real *phiscale, Real *Htrunc, Integer *niter, Integer *subsample, 
       Integer *burn_in, Real *ss_sigma, Integer *df, Integer *phinr,  Real *phi_discrete, Integer *nmphi, Real *e_mean, Real *lambda,
       Real *SS, Real *phisamples){
  
  typedef Real *Doublearray;  
  #define PRN 1000           
  Integer i, j, jj, l, k, acc, acc_phi, dim2, itr, jstep, jphi, ctrl1, ctrl2, ctrl3; 
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop, ss, ssprop, phi, phiprop, phistep, logp_phi, logp_phiprop ; 
  Doublearray z, zprop, S, Sprop, gradz, gradzprop, temp;
  Doublearray Q, Qprop, AA;
  Doublearray BB[2001];
 
  z = Salloc((*n),Real) ;
  zprop = Salloc((*n),Real) ;
  S = Salloc((*n),Real) ;
  Sprop = Salloc((*n),Real) ;
  gradz = Salloc((*n),Real) ;
  gradzprop = Salloc((*n),Real) ;
  temp = Salloc((*n),Real) ;
  dim2 = ((*n)*(*n)+(*n))/2 ;
  Q = Salloc(dim2,Real) ;
  Qprop = Salloc(dim2,Real) ;
  AA = Salloc(dim2,Real) ;

  if(*nmphi>2001){
       printf(" number of values in phi.discrete must not exceed 2001 \n");
       exit(1);
  }
  for (l=0; l<(*n); l++) S[l]=SS[l];
  for (l=0; l<(*n); l++){
      z[l]=gradz[l]=gradzprop[l]=zprop[l]=Sprop[l]=0;
  }
  for (jj=0; jj<dim2; jj++){
      AA[jj]=Q[jj]=Qprop[jj]=0;   
  }    
  phi = phisamples[0];
  if(*nmphi > 1){
      phistep=phi_discrete[1]-phi_discrete[0];
      jphi=floor((phi-phi_discrete[0])/phistep+0.5);
      if(jphi<0) jphi=0;   
      if(jphi>(*nmphi)-1) jphi=(*nmphi)-1;          
      for (j=0;j<(*nmphi);j++){  
         phi=phi_discrete[0]+j*phistep;
         BB[j] = Salloc(dim2,Real) ;
         for (l=0; l<(*n); l++){
	    for (k=0; k<l; k++){
               AA[(*n)*k - k*(k+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*k - k*(k+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+k];
            }
            AA[(*n)*l - l*(l+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*l - l*(l+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+l] + (*tausq);
         }
         cholesky(AA,BB[j],(*n));
      }  
      phi=phi_discrete[0]+jphi*phistep; 
      Q=BB[jphi];
  }
  else{
      phistep=0;   /* redundant variable in this case */
      jphi=0;
      for (l=0; l<(*n); l++){
	 for (k=0; k<l; k++){
            AA[(*n)*k - k*(k+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*k - k*(k+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+k];
         }
         AA[(*n)*l - l*(l+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*l - l*(l+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+l] + (*tausq); 
      }
      cholesky(AA,Q,(*n)); 
  }  
  initz(S,Q,z,(*n));
  conddensity5boxcox(S,Q,&logf,data,z,meanS,units,(*n),(*lambda));    
  ss = calc_ss(z,(*n))+(*ss_sigma);  
  gradient5boxcox(S,gradz,Q,z,data,meanS,units,Htrunc,(*n),ss,(*df),(*lambda));  
  logp_z = -0.5*(*df)*log(ss); 
  logp_phi = logprior_phi(phi, (*e_mean), (*phinr));
  itr = (*niter) + (*burn_in) ;
  acc= 0; acc_phi = 0; ctrl2 = 0; ctrl3 = 0; 
  RANDIN;

  for (i=0; i<itr; i++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + RNORM*sqrt(*scale); 
        conddensity5boxcox(Sprop,Q,&logfprop,data,zprop,meanS,units,(*n),(*lambda));
        ctrl1 = 0;
        for(l=0;l<(*n);l++){
           if(Sprop[l]+meanS[l] < -1/(*lambda)){
               if(data[l] > 0) ctrl1 = 1 ;
           }
        }
        if(ctrl1 == 0){   
           ssprop = calc_ss(zprop,(*n))+(*ss_sigma);  
           gradient5boxcox(Sprop,gradzprop,Q,zprop,data,meanS,units,Htrunc,(*n),ssprop,(*df),(*lambda));  
           logp_zprop = -0.5*(*df)*log(ssprop); 
           for (logq=0,logqprop=0,l=0; l<(*n); l++){
              logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
              logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
           }
           logq*=(-0.5/(*scale));
           logqprop*=(-0.5/(*scale)); 
           if (log(UNIF)<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	       /*accept */
               logf=logfprop;
               logp_z=logp_zprop;
               ss = ssprop;
               temp=gradz;
               gradz=gradzprop;
               gradzprop=temp;
               temp=S;
               S=Sprop;
               Sprop=temp;
               temp=z;
               z=zprop;
               zprop=temp;
               acc++;
           }  
	}
	else ctrl2++;
        if(*nmphi>1) jstep=floor(RNORM*sqrt(*phiscale)/phistep+0.5);
        else jstep=0;
        if ( jphi+jstep > -1 && jphi+jstep < (*nmphi) && !(jstep == 0)){  /* if within limits */  
            phiprop=phi+phistep*jstep;
            logp_phiprop = logprior_phi(phiprop, (*e_mean), (*phinr));
            conddensity5boxcox(Sprop,BB[jphi+jstep],&logfprop,data,z,meanS,units,(*n),(*lambda));
	    ctrl1 = 0;
	    for(l=0;l<(*n);l++){
	      if(Sprop[l]+meanS[l] < -1/(*lambda)){
		if(data[l] > 0) ctrl1 = 1 ;
	      }
	    }
	    if(ctrl1 == 0){ 
	      if (log(UNIF)<logfprop-logf+logp_phiprop-logp_phi){    
		phi=phiprop;
		logf=logfprop;
		logp_phi=logp_phiprop;
		jphi+=jstep;
		Q=BB[jphi];   
		temp=S;
		S=Sprop;
		Sprop=temp;
		acc_phi++;
		gradient5boxcox(S,gradz,Q,z,data,meanS,units,Htrunc,(*n),ssprop,(*df),(*lambda));  
	      } 
	    }
	    else ctrl3++;
        }         
        if (((i+1-(*burn_in))%(*subsample))==0 && (i+1)>(*burn_in)){
              for (l=0; l<(*n); l++) SS[((i+1-(*burn_in))/(*subsample)-1)*(*n)+l] = S[l]; 
              phisamples[((i+1-(*burn_in))/(*subsample)-1)]=phi; 
        }
	if((i+1)==(*burn_in) && (*burn_in)>0){
	  if(*nmphi > 1)
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",(*burn_in), (Real) acc/(*burn_in), 
		   (Real) acc_phi/(*burn_in));
	  else
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f \n",(*burn_in), (Real) acc/(*burn_in));
	  acc =0 ; acc_phi =0 ;
	}
	if((i+(*burn_in)+1)%PRN==0 && (i+1)>(*burn_in)){  
	  if(*nmphi > 1)
	    printf("iter. numb. %d ; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",i+1,(Real) acc/PRN, (Real) acc_phi/PRN);
	  else
	    printf("iter. numb. %d ; Acc.-rate = %1.2f \n",i+1,(Real) acc/PRN);
	  acc =0 ; acc_phi =0 ;
	}   
  }
  if(ctrl2*10 > itr || ctrl3*10 > itr){
    printf("rejection of proposals cauced by density for proposal being zero: S: %d; phi: %d , out of %d iterations \n",ctrl2,ctrl3,itr);
  }
  RANDOUT;  
}


void conddensity5(Real *S, Real *BB, Real *logfcond, Real *obsdata, Real *z, Real *linpred, Real *units, Integer dim){
Integer l, k;

  for (l=0; l<dim; l++){
      for (S[l]=0,k=0; k<=l; k++) S[l]+=BB[dim*k - k*(k+1)/2 + l]*z[k];
   }  
   for (*logfcond=0,l=0; l<dim; l++){
      *logfcond+= (obsdata[l]*(S[l]+linpred[l])-units[l]*exp(S[l]+linpred[l]));
   }
}


void gradient5(Real *S, Real *gradz, Real *BB, Real *z, Real *obsdata, Real *linpred, Real *units, Real *Hz, Integer dim, Real ss, Integer df){
   Integer l, k;
   Real likeli;
   
   for (l=0; l<dim; l++) gradz[l]=-z[l]*df/ss;
   for (k=0; k<dim; k++){
      likeli = obsdata[k] -trunc_u(units[k]*exp(S[k]+linpred[k]),Hz[k]);    
      for (l=0; l<k+1; l++) gradz[l]+=BB[dim*l-l*(l+1)/2+k]*likeli;  
   }
}

 
void mcmcrun5(Integer *n, Real *data, Real *units, Real *meanS, Real *DDvbetaDD, Integer *cornr, Real *kappa, Real *tausq,
       Real *distcoords, Real *scale, Real *phiscale,Real *Htrunc, Integer *niter, Integer *subsample, 
       Integer *burn_in, Real *ss_sigma, Integer *df, Integer *phinr,  Real *phi_discrete, Integer *nmphi, Real *e_mean, Real *SS, 
       Real *phisamples){
  
  typedef Real *Doublearray;  
  #define PRN 1000           
  Integer i, j, jj, l, k, acc, acc_phi, dim2, itr, jstep, jphi ; 
  Real logfprop, logp_zprop,logf, logp_z, logq, logqprop, ss, ssprop, phi, phiprop, phistep, logp_phi, logp_phiprop ; 
  Doublearray z, zprop, S, Sprop, gradz, gradzprop, temp;
  Doublearray Q, Qprop, AA;
  Doublearray BB[2001];

  
  z = Salloc((*n),Real) ;
  zprop = Salloc((*n),Real) ;
  S = Salloc((*n),Real) ;
  Sprop = Salloc((*n),Real) ;
  gradz = Salloc((*n),Real) ;
  gradzprop = Salloc((*n),Real) ;
  temp = Salloc((*n),Real) ;
  dim2 = ((*n)*(*n)+(*n))/2 ;
  Q = Salloc(dim2,Real) ;
  Qprop = Salloc(dim2,Real) ;
  AA = Salloc(dim2,Real) ;

  if(*nmphi>2001){
       printf(" number of values in phi.discrete must not exceed 2001 \n");
       exit(1);
  }
  for (l=0; l<(*n); l++) S[l]=SS[l];
  for (l=0; l<(*n); l++){
      z[l]=gradz[l]=gradzprop[l]=zprop[l]=Sprop[l]=0;
  }
  for (jj=0; jj<dim2; jj++){
      AA[jj]=Q[jj]=Qprop[jj]=0;   
  }    
  phi = phisamples[0];
  if(*nmphi > 1){
      phistep=phi_discrete[1]-phi_discrete[0];
      jphi=floor((phi-phi_discrete[0])/phistep+0.5);
      if(jphi<0) jphi=0;   
      if(jphi>(*nmphi)-1) jphi=(*nmphi)-1;          
      for (j=0;j<(*nmphi);j++){  
         phi=phi_discrete[0]+j*phistep;
         BB[j] = Salloc(dim2,Real) ;
         for (l=0; l<(*n); l++){
	    for (k=0; k<l; k++){
               AA[(*n)*k - k*(k+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*k - k*(k+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+k];
            }
            AA[(*n)*l - l*(l+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*l - l*(l+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+l] + (*tausq); 
         }
         cholesky(AA,BB[j],(*n));
      }  
      phi=phi_discrete[0]+jphi*phistep; 
      Q=BB[jphi];
  }
  else{
      phistep=0;   /* redundant variable in this case */
      jphi=0;
      for (l=0; l<(*n); l++){
	 for (k=0; k<l; k++){
            AA[(*n)*k - k*(k+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*k - k*(k+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+k];
         }
         AA[(*n)*l - l*(l+1)/2+l]=corrfct(phi,(*kappa),distcoords[(*n)*l - l*(l+1)/2+l],(*cornr)) + DDvbetaDD[l*(*n)+l] + (*tausq); 
      }
      cholesky(AA,Q,(*n)); 
  }  
  initz(S,Q,z,(*n));
  conddensity5(S,Q,&logf,data,z,meanS,units,(*n));    
  ss = calc_ss(z,(*n))+(*ss_sigma);  
  gradient5(S,gradz,Q,z,data,meanS,units,Htrunc,(*n),ss,(*df));  
  logp_z = -0.5*(*df)*log(ss); 
  logp_phi = logprior_phi(phi, (*e_mean), (*phinr));
  itr = (*niter) + (*burn_in) ;
  acc= 0; acc_phi = 0 ; 
  RANDIN;
 
  for (i=0; i<itr; i++){  
        for (l=0; l<(*n); l++) zprop[l]=z[l]+0.5*gradz[l]*(*scale) + RNORM*sqrt(*scale); 
        conddensity5(Sprop,Q,&logfprop,data,zprop,meanS,units,(*n));    
        ssprop = calc_ss(zprop,(*n))+(*ss_sigma);  
        gradient5(Sprop,gradzprop,Q,zprop,data,meanS,units,Htrunc,(*n),ssprop,(*df));  
        logp_zprop = -0.5*(*df)*log(ssprop); 
        for (logq=0,logqprop=0,l=0; l<(*n); l++){
           logq+=pow(zprop[l]-(z[l]+0.5*gradz[l]*(*scale)),2);
           logqprop+=pow(z[l]-(zprop[l]+0.5*gradzprop[l]*(*scale)),2);
        }
        logq*=(-0.5/(*scale));
        logqprop*=(-0.5/(*scale)); 
        if (log(UNIF)<logfprop+logp_zprop+logqprop-logq-logf-logp_z){  
 	    /*accept */
            logf=logfprop;
            logp_z=logp_zprop;
            ss = ssprop;
            temp=gradz;
            gradz=gradzprop;
            gradzprop=temp;
            temp=S;
            S=Sprop;
            Sprop=temp;
            temp=z;
            z=zprop;
            zprop=temp;
            acc++;
        }       
        if(*nmphi>1) jstep=floor(RNORM*sqrt(*phiscale)/phistep+0.5);
        else jstep=0;
        if ( jphi+jstep > -1 && jphi+jstep < (*nmphi) && !(jstep == 0)){  /* if within limits */  
            phiprop=phi+phistep*jstep;
            logp_phiprop = logprior_phi(phiprop, (*e_mean), (*phinr));
            conddensity5(Sprop,BB[jphi+jstep],&logfprop,data,z,meanS,units,(*n));
            if (log(UNIF)<logfprop-logf+logp_phiprop-logp_phi){    
               phi=phiprop;
               logf=logfprop;
               logp_phi=logp_phiprop;
               jphi+=jstep;
               Q=BB[jphi];   
               temp=S;
               S=Sprop;
               Sprop=temp;
               acc_phi++;
               gradient5(S,gradz,Q,z,data,meanS,units,Htrunc,(*n),ssprop,(*df));  
            } 
        }         
        if (((i+1-(*burn_in))%(*subsample))==0 && (i+1)>(*burn_in)){
              for (l=0; l<(*n); l++) SS[((i+1-(*burn_in))/(*subsample)-1)*(*n)+l] = S[l]; 
              phisamples[((i+1-(*burn_in))/(*subsample)-1)]=phi; 
        }
	if((i+1)==(*burn_in) && (*burn_in)>0){
	  if(*nmphi > 1)
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",(*burn_in), (Real) acc/(*burn_in), 
		   (Real) acc_phi/(*burn_in));
	  else
	    printf("burn-in =  %d is finished; Acc.-rate = %1.2f \n",(*burn_in), (Real) acc/(*burn_in));
	  acc =0 ; acc_phi =0 ;
	}
	if((i+(*burn_in)+1)%PRN==0 && (i+1)>(*burn_in)){  
	  if(*nmphi > 1)
	    printf("iter. numb. %d ; Acc.-rate = %1.2f ; Acc-rate-phi = %1.2f \n",i+1,(Real) acc/PRN, (Real) acc_phi/PRN);
	  else
	    printf("iter. numb. %d ; Acc.-rate = %1.2f \n",i+1,(Real) acc/PRN);
	  acc =0 ; acc_phi =0 ;
	}   
  }
  RANDOUT;  
}




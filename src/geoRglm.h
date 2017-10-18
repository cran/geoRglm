
extern void binitprod(long *n, double *xc, double *yc, double *sim, long *nbins, double *lims, double *maxdist, long *cbin, double *vbin);

extern void mcmcrun4(long *n, double *data, double *units, double *DD, long *no_linpar, long *cornr, double *kappa, double *tausq,
	      double *coords1, double *coords2, double *scale, double *phiscale, double *Htrunc, long *niter, long *subsample, 
	      long *burn_in, long *messages, double *ss_sigma, long *df, double *phiprior,  double *phi_discrete, long *nmphi, double *SS,
	      double *phisamples, double *acc_rate, double *acc_rate_phi);

extern void mcmcrun4boxcox(long *n, double *data, double *units, double *DD, long *no_linpar, long *cornr, double *kappa, double *tausq,
		    double *coords1, double *coords2, double *scale, double *phiscale, double *Htrunc, long *niter, long *subsample, 
		    long *burn_in, long *messages, double *ss_sigma, long *df, double *phiprior,  double *phi_discrete, long *nmphi, 
		    double *lambda, double *SS, double *phisamples, double *acc_rate, double *acc_rate_phi);

extern void mcmcrun4binom(long *n, double *data, double *units, double *DD, long *no_linpar, long *cornr, double *kappa, double *tausq,
		   double *coords1, double *coords2, double *scale, double *phiscale, long *niter, long *subsample, 
		   long *burn_in, long *messages, double *ss_sigma, long *df, double *phiprior,  double *phi_discrete, long *nmphi, 
		   double *SS, double *phisamples, double *acc_rate, double *acc_rate_phi);

extern void mcmcrun5binom(long *n, double *data, double *units, double *meanS, double *DDvbetaDD, long *cornr, double *kappa, double *tausq,
		   double *coords1, double *coords2, double *scale, double *phiscale, long *niter, long *subsample, 
		   long *burn_in, long *messages, double *ss_sigma, long *df, double *phiprior,  double *phi_discrete, long *nmphi, 
		   double *SS, double *phisamples, double *acc_rate, double *acc_rate_phi);

extern void mcmcrun5boxcox(long *n, double *data, double *units, double *meanS, double *DDvbetaDD, long *cornr, double *kappa, double *tausq,
       double *coords1, double *coords2, double *scale, double *phiscale, double *Htrunc, long *niter, long *subsample, 
       long *burn_in, long *messages, double *ss_sigma, long *df, double *phiprior,  double *phi_discrete, long *nmphi, double *lambda,
       double *SS, double *phisamples, double *acc_rate, double *acc_rate_phi);

extern void mcmcrun5(long *n, double *data, double *units, double *meanS, double *DDvbetaDD, long *cornr, double *kappa, double *tausq,
	      double *coords1, double *coords2, double *scale, double *phiscale,double *Htrunc, long *niter, long *subsample, 
	      long *burn_in, long *messages, double *ss_sigma, long *df, double *phiprior,  double *phi_discrete, long *nmphi, double *SS, 
	      double *phisamples, double *acc_rate, double *acc_rate_phi);

extern void mcmc1poislog(long *n, double *zz, double *SS, double *data, double *meanS, double *QQ, double *QtivQ,
         double *randnormal, double *randunif, double *Htrunc, double *scale, long *nsim ,long *subsample, double *acc_rate);

extern void mcmc1poisboxcox(long *n, double *zz, double *SS, double *data, double *units, double *meanS, double *QQ, double *QtivQ,
         double *randnormal, double *randunif, double *Htrunc, double *scale, long *nsim, long *subsample, double *lambda, double *acc_rate);

extern void mcmc1binom(long *n, double *zz, double *SS, double *data, double *units, double *meanS, double *QQ, double *QtivQ,
         double *randnormal, double *randunif, double *scale, long *nsim, long *subsample, double *acc_rate);


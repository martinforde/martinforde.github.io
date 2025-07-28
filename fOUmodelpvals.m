%type following commands into command window in MATLAB to run
%clear all;
%global V X eta;
%SPXRealizedVariance1minIntervals03Jan22to15Jul24
%X=log(V/V(1));X=X(2:end);
%options = optimset('Display','iter','PlotFcns',@optimplotfval);options.MaxIter=100;options.TolFun=1e-12;options.TolX=1e-12;
%params=[.1 1.5 .01 1]
%fminsearch(@fOUMLEfunction,params,options)
%params=[ans,eta]
%fOUmodelpvals(params)

function fun=fOUmodelpvals(params) 

global V 
 
H=params(1);nu=params(2);lm=params(3);xi=params(4);

N=length(V);T=(N-1)/252;

tau=linspace(0,T,N); 

%kappa_0=lm^(-2*H)*nu^2*gamma(1+2*H);
R=real(.5*exp(-lm*tau).*(-lm^2)^(-2*H).*(-2*exp(lm*tau).*(-lm^2*tau).^(2*H)... 
         +((-lm)^(2*H) - lm^(2*H))*gamma(1 + 2*H) + lm^(2*H)*gamma_inc(1 + 2*H, -lm*tau) +  ...
         exp(2*lm*tau).*(-lm)^(2*H).*gamma_inc(1 + 2*H,lm*tau)));
Cov=toeplitz(R);
C=transpose(chol(Cov));
%Z_=inv(C)*(log(V/xi)+.5*kappa_0)/nu;
Z_=inv(C)*log(V/xi)/nu;

[H_, pvalSW, W] = swtest(Z_);
[h,pvalKS,ksstat,cv] =kstest(Z_);

[pvalKS,pvalSW]
normalitytest(transpose(Z_))


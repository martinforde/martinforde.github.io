%type following commands into command window in MATLAB to run
%clear all;
%global V X eta;
%SPXRealizedVariance1minIntervals03Jan22to15Jul24
%X=log(V/V(1));X=X(2:end);
%options = optimset('Display','iter','PlotFcns',@optimplotfval);options.MaxIter=100;options.TolFun=1e-12;options.TolX=1e-12;
%fminsearch(@fBMMLEFunction,.05,options)
%params=[ans,eta]
%RFSVModelpvals(params)

function fun=RFSVModelpvals(params) 
 

global V

H=params(1);nu=params(2);
N=length(V)-1;T=N/252;

NN=linspace(0,N-1,N);  
R=.5*(abs(NN+1).^(2*H)-2*abs(NN).^(2*H)+abs(NN-1).^(2*H));
Cov=toeplitz(R);
%C=transpose(toeplitz_cholesky_lower(N,Cov));
C=transpose(chol(Cov));
 
dBH_=diff(log(V/V(1))/nu*(N/T)^H);
Z_=inv(C)*dBH_;

%qqplot of Z_ values
qqplot(Z_)
[H_, pvalSW, W] = swtest(Z_);
[h,pvalKS,ksstat,cv] =kstest(Z_);

%Kolmogorov-Smirnov+ Shaprio-Wilks normality tests
disp('KS,SW p-vals=')
[pvalKS,pvalSW]

disp('Other normality tests:')
normalitytest(transpose(Z_)) 

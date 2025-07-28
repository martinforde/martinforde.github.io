
function fun=fOUMLEfunction(params,V)
global V  
 
H=params(1);nu=params(2);lm=params(3);xi=params(4);
H=max(min(H,.499),.001);nu=max(nu,0.0001);xi=max(xi,0.0001);
lm=max(lm,.00001);
  
N = length(V);
T = (N-1)/252;
dt = T/N;

X=log(V/xi);
 
tau=linspace(0,T,N);
%acf computed in closed-form using Garnier+Solna integral formula and Simplify command in Mathematica
R=real(.5*exp(-lm*tau).*(-lm^2)^(-2*H).*(-2*exp(lm*tau).*(-lm^2*tau).^(2*H)... 
         +((-lm)^(2*H) - lm^(2*H))*gamma(1 + 2*H) + lm^(2*H)*gamma_inc(1 + 2*H, -lm*tau) +  ...
         exp(2*lm*tau).*(-lm)^(2*H).*gamma_inc(1 + 2*H,lm*tau)));

Cov=toeplitz(R);

m=390;
%m needed if we wish to make market microstructure noise adjustment,
%appears to make little difference in practice
%LL=mtimes(transpose(X+.5*nu^2),mtimes(inv(nu^2*Cov+0/m),X+.5*nu^2));
LL=mtimes(transpose(X),mtimes(inv(nu^2*Cov+0/m),X));

%if det gives NaN try integer number higher than 1 for c1 until
%problems;  implementation uses that det(c1 A)= c1^N A where N
%is size of matrix A, and we then subtract off N log c_1 to negate this

c1=2;
R_=c1*Cov*nu^2;

fun=LL+log(det(R_))-N*log(c1);
[H nu xi lm]
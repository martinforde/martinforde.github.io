function fun=RbarSeries(tau,alpha,c,lambda)
 
SS=0;
 
 
N=30;
for n=0:N
    %SS=SS+c*tau^alpha*(-c*tau^alpha)^n/(alpha*(1+n)*gamma((1+n)*alpha));
    SS=SS+(-c)^(1 + n)*lambda^(-(1 + n)*alpha)*gamma_inc((1 + n)*alpha,lambda*tau)/gamma((1 + n)*alpha) - ...
        (-c)^(1 + n)*lambda^(-(1 + n)*alpha);
end
 
fun=SS;
 
end


function fun=gamma_inc(a,x)

fun=gamma(a)*(1-gammainc(x,a));

end

clear all

tic

NumSimulations = 50000;
NumTimeSteps = 200;
T = .0001;
T=.01;
V0=1;
H=.25;
Delta=1;

gm=.1;
p=-1;
 
Payoff=0;

dt=T/NumTimeSteps;
Sqrtdt=sqrt(dt);


M=20;
[uu,w]=lgwt(M,T,T+Delta);


for i=1:NumSimulations
 
    VIXSqrd=V0;
    
    dW=Sqrtdt*randn([NumTimeSteps 1]);
    
    for m=1:M
        logxi_t(m)=0.0;
    end 
    
    
    for j=1:NumTimeSteps
        
        WH(j) = 0; 
        t=j*dt;

        for k=1:j-1
            
            s=k*dt; 
            WH(j) = WH(j) + (t-s)^(H-.5)*dW(k);
            
        end   
        
        V(j)=exp(gm*WH(j)-.5*gm^2*t^(2*H)/(2*H));
        
        for m=1:M
            u=uu(m);
            logxi_t(m)=logxi_t(m)+gm*(u-t)^(H-.5)*dW(j);
        
        end 
           
    end
    
    t=T;
    
    VIXSqrd=0;
    
    for m=1:M
        u=uu(m);
        logxi_t(m)=logxi_t(m)-.25*gm^2*(u^(2*H)-(u-t)^(2*H))/(2*H);
        VIXSqrd=VIXSqrd+1/Delta*w(m)*exp(logxi_t(m));
        
    end 
            
    Payoff = Payoff + exp(p/T*(VIXSqrd-1));
    
end
 
T*log(Payoff/NumSimulations)
toc
%options = optimset('Display','iter','PlotFcns',@optimplotfval);options.MaxIter=5000;options.TolFun=1e-12;options.TolX=1e-12;
%fminsearch(@fun,params,options)
 
clear all;

M =250000;
N = 256;

%Parameters and market data from CBOE datashop https://datashop.cboe.com/option-eod-summary for 21 June 2024
alpha=0.5624027280856119;a=0.3209753476974414;
b=0.05093279984743844;c=0.00435798589530962;lmbda=5.136393131211633;theta=-0.09215241691687982; 

ImpVolMkt=[0.30833965626324,	0.295476337981542,	0.283463677875459,	0.271111343370089,	0.258451404309255,	0.245507644865208,	0.233694271876973,	0.221358004433688,	0.209097873714426,	0.197249412151294,	0.186334324110301,	0.174675684455109,	0.164666282229566,	0.15488794443707,	0.145823620475546,	0.137304198281544,	0.129231575271196,	0.121548584918863,	0.114417742487412,	0.108019587373161,	0.102797547251224,	0.0988783071539835,	0.0959001836718852,	0.0934990753558453,	0.0918556359282457,	0.0918871975138243,	0.0934062296063838];
ImpVolMkt2=[0.256853	0.247760	0.239149	0.230164	0.221471	0.213216	0.204941	0.196783	0.188797	0.180978	0.173596	0.166311	0.159104	0.152338	0.145637	0.139202	0.132974	0.127012	0.121383	0.116176	0.111556	0.108040	0.105025	0.102493	0.100738	0.099360	0.098370];

ImpVolBid=[0.307106846	0.294339683	0.282433411	0.270174503	0.257597541 0.244728216	0.232999977	0.22073676	0.208545866	0.196765632	0.185920088	0.17396049	0.164071874	0.154398142	0.145424921	0.136980885	0.128836002	0.121330108	0.114137677	0.107516176	0.102309591	0.098442181	0.095584017	0.093062924	0.091627293	0.091493972	0.092700753];
ImpVolAsk=[0.309547279	0.296591262	0.284475978	0.272033194	0.259292671	0.246276428	0.234380125	0.221972468	0.209644531	0.197729126	0.186745673	0.175382578	0.165255352	0.155374461	0.146220437	0.137626503	0.129626025	0.12176684	0.114697629	0.108522865	0.103285495	0.09931414	0.09621579	0.093932825	0.092082793	0.092275133	0.094089284];
ImpVolBid2=[0.255952	0.247347	0.238393	0.229819	0.220841	0.212648	0.204429	0.196553	0.188387	0.180613	0.173273	0.166026	0.158853	0.152117	0.145345	0.139031	0.132744	0.126805	0.121067	0.115819	0.111268	0.107749	0.104718	0.102287	0.100576	0.099158	0.097973];
ImpVolAsk2=[0.257741	0.248171	0.239896	0.230508	0.222095	0.213780	0.205450	0.197011	0.189205	0.181341	0.173917	0.166594	0.159355	0.152558	0.145928	0.139374	0.133203	0.127219	0.121700	0.116532	0.111844	0.108330	0.105333	0.102699	0.100899	0.099560	0.098764];

VIXImpVolMkt=[0.502489195 0.531472613 0.580483519 0.63710948 0.68409345 0.736859132 0.787044999 0.827805687 0.870005448 0.94647844 1.020997407 1.09475647 1.167988237 1.219598942];
VIXImpVolBid=[0.492706177	0.519256444	0.571592221	0.625867399	0.673715212	0.72339593	0.773522897	0.817418158	0.859238161	0.934672268	1.012307439	1.080370662	1.157540911	1.208042188];
VIXImpVolAsk=[0.511921243	0.543408429	0.589315274	0.648324234	0.694469382	0.750324131	0.800559657	0.838178933	0.88074487	0.958219175	1.029636149	1.10896741	1.178328718	1.231004519];
VIXFutureMkt=.1437955868051;

tic
lambda=lmbda;nu=1;Z0=0;
H=alpha-.5;g0=Z0;
Delta=1/12;

c1=1/gamma(alpha);

T=20/251;T1=T;
F=5489.83;F2=5509.62;
K=linspace(4500,5800,27)/F;K2=linspace(4500,5800,27)/F2;
K_VIX=[0.1200  0.1250    0.1300    0.1350    0.1400    0.1450    0.1500    0.1550    0.1600  .17 .18 .19 .2 .21];

%fix seed rng(1), seed should always be fixed for calibration, for same reason seed should be fixed when estimating Greeks using Monte Carlo 

V=zeros(M,N+1);
V2=zeros(M,N+1);
Z1=zeros(M,N+1);Z2=zeros(M,N+1);

dt = T/N; sqrtdt = sqrt(dt);
Z=randn([M N]);
dW=sqrtdt*Z;

X=zeros(M,1);
X2=X;
V(:,1)=a*(Z0-b)^2+c;
sqrt(V(1,1));
V2(:,1)=V(:,1);

tmp=4^(-H)*lambda^(-2*H);
tmp2=lambda^(-.5-H);

for j=2:N+1
       j;
       t=(j-1)*dt; 
       Z1(:,j)=Z0;
       Z2(:,j)=Z0;

       for k=1:j-1

            s=k*dt;
            sigma=tmp*gamma_inc(2*H, 2*(t-s)*lambda);
            s=(k-1)*dt;
            sigma=sigma-tmp*gamma_inc(2*H, 2*(t-s)*lambda);
            sigma=sqrt(sigma);

            Z1(:,j) =Z1(:,j) + c1*nu*sigma*sqrt(V(:,k)).*Z(:,k);  
            Z2(:,j) =Z2(:,j) - c1*nu*sigma*sqrt(V2(:,k)).*Z(:,k);  
            
       end 

       tmp3=theta*tmp2*(gamma(.5+H)-gamma_inc(.5+H, t*lambda));
       %tmp3=0;
       Z1(:,j) =Z1(:,j)+tmp3;
       Z2(:,j) =Z2(:,j)+tmp3;

       V(:,j)=a*(Z1(:,j)-b).^2+c;
       V2(:,j)=a*(Z2(:,j)-b).^2+c;
       X=X-.5*V(:,j-1)*dt+sqrt(V(:,j-1)).*dW(:,j-1);
       X2=X2-.5*V2(:,j-1)*dt-sqrt(V2(:,j-1)).*dW(:,j-1);

       S=exp(X);
       S2=exp(X2);

end

CallPrices=mean(.5*(max(S-K,0)+max(S2-K,0)));
CallPricesStd=std(.5*(max(S-K,0)+max(S2-K,0)));
CallPrices1=CallPrices;

for i=1:length(CallPrices)
    if K(i)>=1
       RSE(i)=norminv(.95)*CallPricesStd(i)/CallPrices1(i)/sqrt(M);
    else
       RSE(i)=norminv(.95)*CallPricesStd(i)/((CallPrices1(i)+K(i)-1))/sqrt(M);
    end
end

ImpVol=blsimpv(1,K,0,T1,CallPrices)
RSE

V=zeros(M,N+1);
V2=zeros(M,N+1);

T=17/251;T3=T;

dt = T/N; sqrtdt = sqrt(dt);
dW=sqrtdt*Z;

V(:,1)=a*(Z0-b)^2+c;
V2(:,1)=V(:,1);

for j=2:N+1
       t=(j-1)*dt; 
       Z1(:,j)=Z0;
       Z2(:,j)=Z0;

       for k=1:j-1

            s=k*dt;
            sigma=tmp*gamma_inc(2*H, 2*(t-s)*lambda);
            s=(k-1)*dt;
            sigma=sigma-tmp*gamma_inc(2*H, 2*(t-s)*lambda);
            sigma=sqrt(sigma);

            Z1(:,j) =Z1(:,j) + c1*nu*sigma*sqrt(V(:,k)).*Z(:,k);  
            Z2(:,j) =Z2(:,j) - c1*nu*sigma*sqrt(V2(:,k)).*Z(:,k);  
            
       end 

       tmp3=theta*tmp2*(gamma(.5+H)-gamma_inc(.5+H, t*lambda));
       Z1(:,j) =Z1(:,j)+tmp3;
       Z2(:,j) =Z2(:,j)+tmp3;

       V(:,j)=a*(Z1(:,j)-b).^2+c;
       V2(:,j)=a*(Z2(:,j)-b).^2+c;

end


NumPts=20;[u,w]=lgwt(NumPts,T,T+Delta);u=flip(u);

g=zeros(M,NumPts);
g2=zeros(M,NumPts);

g0=Z0;

t=T;
eta=nu;
alpha_star=2*alpha-1;
c_star=-a*eta^2*gamma(alpha_star)/gamma(alpha)^2;

for j=1:NumPts
    tau=u(j)-T;
    %see Eq 109 in chapter 6.2 of Sigurd Roemer article
    % ``Hybrid multifactor scheme for stochastic Volterra equations'' with completely monotone kernel for defN of Rbar
    Rb(j)=RbarSeries(Delta-tau,alpha_star,c_star,2*lambda);
    
    g(:,j)=g0+theta*tmp2*(gamma(.5+H)-gamma_inc(.5+H, (t+tau)*lambda));
    g2(:,j)=g0+theta*tmp2*(gamma(.5+H)-gamma_inc(.5+H, (t+tau)*lambda));

    for k=1:N
        
        s=(k-1)*dt;
        sigma=sqrt(4^(-H)*lambda^(-2*H)*gamma_inc(2*H, 2*(t-k*dt+tau)*lambda)-4^(-H)*lambda^(-2*H)*gamma_inc(2*H, 2*(t-s+tau)*lambda));
       
        g(:,j) =g(:,j) + c1*eta*sigma*sqrt(V(:,k)).*Z(:,k);  
        g2(:,j) =g2(:,j) - c1*eta*sigma*sqrt(V2(:,k)).*Z(:,k);  
    end
end

f=a*(g-b).^2+c;
f2=a*(g2-b).^2+c;

SS=zeros(M,1);SS2=zeros(M,1);

for j=1:NumPts
    SS=SS+(1-Rb(j))*f(:,j)*w(j);
    SS2=SS2+(1-Rb(j))*f2(:,j)*w(j);
end

VIXSqrd=SS/Delta;VIXSqrd2=SS2/Delta;
VIX=sqrt(VIXSqrd);VIX2=sqrt(VIXSqrd2);

VIXFuture=.5*(mean(VIX+VIX2));
%VIXCallPrices=.5*mean(max(VIX-K_VIX,0))+.5*mean(max(VIX2-K_VIX,0));
VIXCallPrices=mean(.5*max(VIX-K_VIX,0)+.5*max(VIX2-K_VIX,0));
VIXCallPricesStd=std(.5*max(VIX-K_VIX,0)+.5*max(VIX2-K_VIX,0));
%.5*std(max(VIX-K_VIX,0))+.5*mean(max(VIX2-K_VIX,0));
VIXImpVol=blsimpv(VIXFuture,K_VIX,0,T3,VIXCallPrices)

%RSE is relative standard error, used to estimate confidence interval
for i=1:length(VIXCallPrices)
    if K_VIX(i)>=VIXFuture
       RSE_VIX(i)=norminv(.95)*VIXCallPricesStd(i)/VIXCallPrices(i)/sqrt(M);
    else
       RSE_VIX(i)=norminv(.95)*VIXCallPricesStd(i)/((VIXCallPrices(i)+K_VIX(i)-VIXFuture))/sqrt(M);
    end
end

RSE_VIX
T=40/251;T2=T;
F2=5509.62;

V=zeros(M,N);
V2=zeros(M,N);
Z1=zeros(M,N);Z2=zeros(M,N);

dt = T/N; sqrtdt = sqrt(dt);
dW=sqrtdt*Z;

X=zeros(M,1);
X2=X;
V(:,1)=a*(Z0-b)^2+c;
V2(:,1)=V(:,1);

for j=2:N+1
       t=(j-1)*dt;
       Z1(:,j)=Z0;
       Z2(:,j)=Z0;

       for k=1:j-1

            s=k*dt;
            sigma=tmp*gamma_inc(2*H, 2*(t-s)*lambda);
            s=(k-1)*dt;
            sigma=sigma-tmp*gamma_inc(2*H, 2*(t-s)*lambda);
            sigma=sqrt(sigma);

            Z1(:,j) =Z1(:,j) + c1*nu*sigma*sqrt(V(:,k)).*Z(:,k);  
            Z2(:,j) =Z2(:,j) - c1*nu*sigma*sqrt(V2(:,k)).*Z(:,k);  
            
       end
       
       tmp3=theta*tmp2*(gamma(.5+H)-gamma_inc(.5+H, t*lambda));
       Z1(:,j) =Z1(:,j)+tmp3;
       Z2(:,j) =Z2(:,j)+tmp3;

       V(:,j)=a*(Z1(:,j)-b).^2+c;
       V2(:,j)=a*(Z2(:,j)-b).^2+c;
       X=X-.5*V(:,j-1)*dt+sqrt(V(:,j-1)).*dW(:,j-1);
       X2=X2-.5*V2(:,j-1)*dt-sqrt(V2(:,j-1)).*dW(:,j-1);

       S=exp(X);
       S2=exp(X2);

end

CallPrices=mean(.5*(max(S-K2,0)+max(S2-K2,0)));
CallPricesStd2=std(.5*(max(S-K2,0)+max(S2-K2,0)));
CallPrices2=CallPrices;

for i=1:length(CallPrices2)
    if K2(i)>=1
       RSE2(i)=norminv(.95)*CallPricesStd2(i)/CallPrices2(i)/sqrt(M);
    else
       RSE2(i)=norminv(.95)*CallPricesStd2(i)/((CallPrices2(i)+K2(i)-1))/sqrt(M);
    end
end

ImpVol2=blsimpv(1,K2,0,T2,CallPrices2)
RSE2

clf
subplot(3,2,1);
plot(K*F,ImpVol,'LineWidth',1.2,'Color',[.2 .4 .75]); 
hold on
plot(K*F,ImpVolMkt,'LineWidth',1.2,'Color',[.85 0 0]); 
plot(K*F,ImpVolBid,'+k','LineWidth',.75,'Color',[.4 .4 .4],'MarkerSize',3);
plot(K*F,ImpVolAsk,'+k','LineWidth',.75,'Color',[.4 .4 .4],'MarkerSize',3);
xline(F,'--','LineWidth',.5,'Color',[.2 .2 .2])
%plot(K*F,ImpVolBid,'-','LineWidth',.65,'Color',[.4 .4 .4]);
%plot(K*F,ImpVolAsk,'-','LineWidth',.65,'Color',[.4 .4 .4]);
grid on
ylim([min(min(ImpVol),min(ImpVolMkt)),max(max(ImpVol),max(ImpVolMkt))])
title("SPX smile T=20/251")
xlim([min(K*F),max(K*F)])
subplot(3,2,2);
plot(K2*F2,ImpVol2,'LineWidth',1.2,'Color',[.2 .4 .75]); 
hold on
plot(K2*F2,ImpVolMkt2,'LineWidth',1.2,'Color',[.85 0 0]); 
plot(K2*F2,ImpVolBid2,'+k','LineWidth',.75,'Color',[.4 .4 .4],'MarkerSize',3);
plot(K2*F2,ImpVolAsk2,'+k','LineWidth',.75,'Color',[.4 .4 .4],'MarkerSize',3);
xline(F2,'--','LineWidth',.5,'Color',[.2 .2 .2])
%plot(K*F2,ImpVolBid2,'-','LineWidth',.65,'Color',[.4 .4 .4]);
%plot(K*F2,ImpVolAsk2,'-','LineWidth',.65,'Color',[.4 .4 .4]);
title("SPX smile T=40/251")
xlim([min(K2*F2),max(K2*F2)])
ylim([min(min(ImpVol2),min(ImpVolMkt2)),max(max(ImpVol2),max(ImpVolMkt2))])
grid on
subplot(3,2,3);
plot(K_VIX,VIXImpVol,'LineWidth',1.2,'Color',[.2 .4 .75]); 
hold on
plot(K_VIX,VIXImpVolMkt,'LineWidth',1.2,'Color',[.85 0 0]); 
plot(K_VIX,VIXImpVolBid,'+k','LineWidth',.75,'Color',[.4 .4 .4],'MarkerSize',3);
plot(K_VIX,VIXImpVolAsk,'+k','LineWidth',.75,'Color',[.4 .4 .4],'MarkerSize',3);
xline(VIXFuture,'LineWidth',.85,'Color',[.2 .4 .75])
xline(VIXFutureMkt,'--','LineWidth',.5,'Color',[.2 .2 .2])
ylim([min(min(VIXImpVol),min(VIXImpVolMkt)),max(max(VIXImpVol),max(VIXImpVolMkt))])
%plot(K_VIX,VIXImpVolBid,'-','LineWidth',.65,'Color',[.4 .4 .4]);
%plot(K_VIX,VIXImpVolAsk,'-','LineWidth',.65,'Color',[.4 .4 .4]);
title("VIX smile T=17/251")
grid on
subplot(3,2,4);
plot(K*F,RSE,'LineWidth',1.2,'Color',[.2 .4 .75])
hold on
plot(K2*F2,RSE2,'LineWidth',1.2,'Color',[.4 .4 .4])
title("Relative Standard error for OTM SPX option prices")
legend('T=20/251','T=40/251')
grid on
subplot(3,2,5);
plot(K_VIX,RSE_VIX,'LineWidth',1.2,'Color',[.2 .4 .75])
title("Relative Standard error for OTM VIX option prices")
grid on
shg
toc
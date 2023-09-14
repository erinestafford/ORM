%clear all; close all; 
%Plots are separately adapted so that lines correspond to number of
%switches instead of numerical order of roots

N = 20;
pop=100.0;
beta = 0.0232/pop;
r = 0.1;%0.11;%0.1 0.1;
c = 0.2;%0.2;%0.23 0.2;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
sigma = 0.1;
I0 = 0.99*pop;
T = 305;
%% beta
out_beta = zeros(N,3);
beta_values = linspace(0.001/pop,.08/pop, N);
n_switches_beta=zeros(N,3);
count=1;

for b = beta_values
    [out_beta(count,:),n_switches_beta(count,:)]=sensitivity_all_roots(I0,b,gamma,r,c,sigma,pop,T,amin,amax);
    count=count+1;
end


%% I0
out_I0 = zeros(N,3);
nswitch_I0 = zeros(N,3);
I0_values = linspace(0.01*pop,.99*pop, N);
count=1;

for i = I0_values
    [temp,nswitch]=sensitivity_all_roots(i,beta,gamma,r,c,sigma,pop,T,amin,amax);
    out_I0(count,:)=temp;
    nswitch_I0(count,:)=nswitch;
    count=count+1;
end

%% gamma
gamma_values = linspace(0.001,0.02,30);
linspace(1/2,1/365,length(gamma_values));
out_gamma = zeros(length(gamma_values),3);
nswitch_gamma = zeros(length(gamma_values),3);
count=1;

for g = gamma_values
    [temp,nswitch]=sensitivity_all_roots(I0,beta,g,r,c,sigma,pop,T,amin,amax);
    out_gamma(count,:)=temp;
    nswitch_gamma(count,:)=nswitch;
    if isnan(out_gamma(count,1))
        out_gamma(count,1)=get_economic_output(I0,beta,g,r,c,sigma,pop,T, amin, amax);
        nswitch_gamma(count,1)=0;
    end
    count=count+1;
end

%% r/c
out_r = zeros(50,3);
nswitch_r = zeros(50,3);
r_values = linspace(0.001*c,c, 50);
count=1;

for rv = r_values
    [temp,nswitch]=sensitivity_all_roots(I0,beta,gamma,rv,c,sigma,pop,T,amin,amax);
    out_r(count,:)=temp;
    nswitch_r(count,:)=nswitch;
    count=count+1;
end


%% sigma
out_sigma = zeros(N,3);
nswitch_sigma = zeros(N,3);
s_values = linspace(0,1, N);
count=1;

for s = s_values
    [temp,nswitch]=sensitivity_all_roots(I0,beta,gamma,r,c,s,pop,T,amin,amax);
    out_sigma(count,:)=temp;
    nswitch_sigma(count,:)=nswitch;
    count=count+1;
end


%% T
out_T = zeros(N,3);
nswitch_T = zeros(N,3);
T_values = linspace(30,365, N);
count=1;

for termT = T_values
    [out_T(count,:),nswitch_T(count,:)]=sensitivity_all_roots(I0,beta,gamma,r,c,sigma,pop,termT,amin,amax);
    count=count+1;
end


%% amin
out_m = zeros(N,3);
nswitch_m = zeros(N,3);
m_values = linspace(0,0.5, N);
count=1;

for m = m_values
    [out_m(count,:),nswitch_m(count,:)]=sensitivity_all_roots(I0,beta,gamma,r,c,sigma,pop,T,m,amax);
    count=count+1;
end

%% amax
out_M = zeros(N,3);
nswitch_M = zeros(N,3);
M_values = linspace(0.5,1, N);
count=1;

for M = M_values
    [out_M(count,:),nswitch_M(count,:)]=sensitivity_all_roots(I0,beta,gamma,r,c,sigma,pop,T,amin,M);
    count=count+1;
end



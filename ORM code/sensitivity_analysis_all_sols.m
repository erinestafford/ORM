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

figure()
plot(beta_values,out_beta(:,1),'--r', 'linewidth',2)
hold on
plot(beta_values,out_beta(:,2),'-b', 'linewidth',2)
plot(beta_values,out_beta(:,3),':g', 'linewidth',2)
xlabel("\beta")
ylabel("Economic Output")

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
figure()
plot(I0_values,out_I0(:,1),'--r', 'linewidth',2)
hold on
plot(I0_values,out_I0(:,2),'-b', 'linewidth',2)
plot(I0_values,out_I0(:,3),'o-g', 'linewidth',2)
xlabel("I_0")
ylabel("Economic Output")
%%
% 
% out2_I0 = zeros(N);
% %nswitch_I0 = zeros(N,3);
% I0_values = linspace(0.01*pop,.999*pop, N);
% count=1;
% 
% for i = I0_values
%     [temp,temp2]=get_economic_output(i,beta,gamma,r,c,sigma,pop,T, amin, amax);
%     temp2
%     out2_I0(count)=temp;
%     count=count+1;
% end
% figure()
% plot(I0_values,out2_I0(:,1),'--r', 'linewidth',2)
% xlabel("I_0")
% ylabel("Economic Output")
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
figure()
plot(gamma_values,out_gamma(:,1),'--r', 'linewidth',2)
hold on
plot(gamma_values,out_gamma(:,2),'-b', 'linewidth',2)
plot(gamma_values,out_gamma(:,3),'o-g', 'linewidth',2)
xlabel("\gamma")
ylabel("Economic Output")
%nans correspont to always active?
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
figure()
plot(r_values./c,out_r(:,1),'--r', 'linewidth',2)
hold on
plot(r_values./c,out_r(:,2),'-b', 'linewidth',2)
plot(r_values./c,out_r(:,3),'o-g', 'linewidth',2)
xlabel("r/c")
ylabel("Economic Output")

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
figure()
plot(s_values,out_sigma(:,1),'--r', 'linewidth',2)
hold on
plot(s_values,out_sigma(:,2),'-b', 'linewidth',2)
plot(s_values,out_sigma(:,3),'o-g', 'linewidth',2)
xlabel("\sigma")
ylabel("Economic Output")

%% T
out_T = zeros(N,3);
nswitch_T = zeros(N,3);
T_values = linspace(30,365, N);
count=1;

for termT = T_values
    [out_T(count,:),nswitch_T(count,:)]=sensitivity_all_roots(I0,beta,gamma,r,c,sigma,pop,termT,amin,amax);
    count=count+1;
end
figure()
plot(T_values,out_T(:,1),'--r', 'linewidth',2)
hold on
plot(T_values,out_T(:,2),'-b', 'linewidth',2)
plot(T_values,out_T(:,3),'o-g', 'linewidth',2)
xlabel("Terminal time,T")
ylabel("Economic Output")

%% amin
out_m = zeros(N,3);
nswitch_m = zeros(N,3);
m_values = linspace(0,0.5, N);
count=1;

for m = m_values
    [out_m(count,:),nswitch_m(count,:)]=sensitivity_all_roots(I0,beta,gamma,r,c,sigma,pop,T,m,amax);
    count=count+1;
end
figure()
plot(m_values,out_m(:,1),'--r', 'linewidth',2)
hold on
plot(m_values,out_m(:,2),'-b', 'linewidth',2)
plot(m_values,out_m(:,3),'o-g', 'linewidth',2)
xlabel("Minimum mixing level, m")
ylabel("Economic Output")

%% amax
out_M = zeros(N,3);
nswitch_M = zeros(N,3);
M_values = linspace(0.5,1, N);
count=1;

for M = M_values
    [out_M(count,:),nswitch_M(count,:)]=sensitivity_all_roots(I0,beta,gamma,r,c,sigma,pop,T,amin,M);
    count=count+1;
end
figure()
plot(M_values,out_M(:,1),'--r', 'linewidth',2)
hold on
plot(M_values,out_M(:,2),'-b', 'linewidth',2)
plot(M_values,out_M(:,3),'o-g', 'linewidth',2)
xlabel("Maximum mixing level, M")
ylabel("Economic Output")


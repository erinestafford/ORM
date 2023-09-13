%clear all; close all;
%% beta
N = 100;
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
%%

out_beta = zeros(N,1);
beta_values = linspace(0.01/pop,.09/pop, N);
count=1;

for b = beta_values
    out_beta(count)=get_economic_output(I0,b,gamma,r,c,sigma,pop,T,amin,amax);
    count=count+1;
end
figure()
plot(beta_values,out_beta,'--r', 'linewidth',2)
xlabel("\beta")
ylabel("Economic Output")

%% I0
out_I0 = zeros(N,1);
I0_values = linspace(0.01*pop,0.99*pop, N);
count=1;

for i = I0_values
    out_I0(count)=get_economic_output(i,beta,gamma,r,c,sigma,pop,T,amin,amax);
    count=count+1;
end
figure()
plot(I0_values,out_I0,'--r', 'linewidth',2)
xlabel("I_0")
ylabel("Economic Output")

%% gamma
N=100;
out_gamma = zeros(N,1);
gamma_values = linspace(1/2,1/365, N);
count=1;

for g = gamma_values
    out_gamma(count)=get_economic_output(I0,beta,g,r,c,sigma,pop,T,amin,amax);
    count=count+1;
end
figure()
plot(gamma_values,out_gamma,'--r', 'linewidth',2)
xlabel("\gamma")
ylabel("Economic Output")

%% r/c
out_r = zeros(N,1);
r_values = linspace(0.01*c,2*c, N);
count=1;

for rv = r_values
    out_r(count)=get_economic_output(I0,beta,gamma,rv,c,sigma,pop,T,amin,amax);
    count=count+1;
end
figure()
plot(r_values./c,out_r,'--r', 'linewidth',2)
xlabel("r/c")
ylabel("Economic Output")

%% sigma
out_sigma = zeros(N,1);
s_values = linspace(0,1, N);
count=1;

for s = s_values
    out_sigma(count)=get_economic_output(I0,beta,gamma,r,c,s,pop,T,amin,amax);
    count=count+1;
end
figure()
plot(s_values,out_sigma,'--r', 'linewidth',2)
xlabel("\sigma")
ylabel("Economic Output")

%% T
N=100;
out_T = zeros(N,1);
T_values = linspace(1,365, N);
count=1;

for termT = T_values
    out_T(count)=get_economic_output(I0,beta,gamma,r,c,sigma,pop,termT,amin,amax);
    count=count+1;
end
figure()
plot(T_values,out_T,'--r', 'linewidth',2)
xlabel("Terminal time,T")
ylabel("Economic Output")

%% amin
out_m = zeros(N,1);
m_values = linspace(0,0.5, N);
count=1;

for m = m_values
    out_m(count)=get_economic_output(I0,beta,gamma,r,c,sigma,pop,T,m,amax);
    count=count+1;
end
figure()
plot(m_values,out_m,'--r', 'linewidth',2)
xlabel("Minimum mixing level, m")
ylabel("Economic Output")

%% amax
out_M = zeros(N,1);
M_values = linspace(0.5,1, N);
count=1;

for M = M_values
    out_M(count)=get_economic_output(I0,beta,gamma,r,c,sigma,pop,T,amin,M);
    count=count+1;
end
figure()
plot(M_values,out_M,'--r', 'linewidth',2)
xlabel("Maximum mixing level, M")
ylabel("Economic Output")


clear all; close all; clc;
%%
pop=1;
beta = 0.06;%0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1/110;
amax = 0.9;
amin = 0.2;
I0 = 0.99;
sigma = 0.5;
l0=-4;
t =0.001:0.01:115;

%%
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
options = odeset('Events',@lzero,'RelTol',1e-6,'AbsTol',1e-6);
[tout,y]=ode45(f,[0 t],[l0,I0],options);
% figure()
% plot(tout,y(:,2))
%%
bcfcn=@(yl,yr)[yl(2)-I0; yr(1)];
guess = @(t)ode45(f,[0 t],[-l0,I0],options).y(:,end);
solinit = bvpinit(t,guess);

sol = bvp4c(f, bcfcn, solinit);

%%
test_sol = deval(sol,t);
plot(t,test_sol(1,:),'linewidth',3)
figure()
plot(t,test_sol(2,:),'linewidth',3)
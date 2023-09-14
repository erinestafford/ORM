clear all; close all; clc;
%% Code to run the T vs lambda plot in the Supplement
pop=100;
beta = 0.0232/pop;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1/110;
amax = 0.9;
amin = 0.2;
I0 = 99;
sigma = 0.1;

l0 = -1.0;
dy = 0.01;
NSTPS = 10000;
t =0:0.01:800;
HT = 20;

%%
options = odeset('Events',@lzero,'MaxStep',.1,'RelTol',1.0e-6,'AbsTol',1.0e-6);
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);

y2 = l0;
l_vals = [];
tt=[];
c = 1;
while y2>-HT
    y = [y2, I0];
    tstart = 0.0;
    tend=0.1;
    for i=0:NSTPS
        tend = tend+0.1;
        [~,Y] = ode45(f,[tstart, tend],y,options);
        y = Y(end,:);
        tstart=tend;
        if Y(end,1)>0
            break;
        end
        if Y(end,1)<-HT
            break;
        end
    end
    if Y(end,1)<-HT
        break;
    end

    l_vals(c) = y2;
    tt(c) = tend;
    c = c+1;
    y2 =y2 - dy;
end

%%
plot(l_vals,tt, 'linewidth',3)
%%
hold on
plot(l_vals, 305*ones(length(tt),1),'linewidth',3)
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
sigma = 0.4;
ds = 0.01;
smax = 0.7;

l0 = -1.0;
dy = 0.01;
NSTPS = 10000;
t =0:0.01:800;
HT = 20;

%%
options = odeset('Events',@lzero,'MaxStep',.1,'RelTol',1.0e-6,'AbsTol',1.0e-6);
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
while sigma<= smax
y2 = l0;
l_vals = [];
tt=[];
count = 1;
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

    l_vals(count) = y2;
    tt(count) = tend;
    count = count+1;
    y2 =y2 - dy;
end
plot3(tt, sigma*ones(length(l_vals),1), l_vals, 'linewidth',3,'color','r')
hold on
sigma = sigma+ds;
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
end

clear all; close all; clc;
%%
pop=1;
beta = 0.06;%0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1/110;
amax = 1.0;
amin = 0.0;
I0 = 0.9;
sigma = 0.5;
term_I = 0.99;

l0 = -1.0;
dy = 0.01;
NSTPS = 10000;
t =0:0.01:800;
HT = 20;

%%
options2 = odeset('Events',@(t,y)term_cond(t,y, term_I),'MaxStep',.1,'RelTol',1.0e-8,'AbsTol',1.0e-8);
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);

y2 = l0;
l_vals = [];
tt=[];
pay_off = [];
count = 1;
while y2>-HT
    y = [y2, I0];
    tstart = 0.0;
    tend=0.1;
    for i=0:NSTPS
        tend = tend+0.1;
        [~,Y] = ode45(f,[tstart, tend],y,options2);
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

    [T,Y] = ode45(f,[0.0, tend],[y2, I0], options2);
    [a,~,~] = get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
    payoff=trapz(T,r*a.*(pop-sigma*Y(:,2)) - c*Y(:,2));

    l_vals(count) = y2;
    tt(count) = tend;
    pay_off(count) = payoff;
    count = count+1;
    y2 =y2 - dy;
end

%%
plot(pay_off,tt, 'linewidth',2)

%%
function [value, isterminal, direction] = term_cond(T, Y,term_I)
value= abs(Y(2)-term_I)<=1e-8;
isterminal = 1;   % Stop the integration
direction  = 0;
end

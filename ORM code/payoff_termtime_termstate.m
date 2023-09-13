%%
close all; clear all;
pop =1;
beta = 0.06;%0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
I0 = 0.7;
term_I = 0.7;
sigma = 0.5;
Tend = 1000;

%%
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
options = odeset('Events', @(t,y)term_cond(t, y, term_I),'RelTol',1e-8);
term_times = [];
payoffs = [];
count = 1;
for l0 = linspace(-50,10,5000)
    y0=[l0,I0];
    [T,Y] = ode45(f,[0, Tend],y0, options);
    if T(end)<1000
        [a,~,~] = get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
        payoff=trapz(T,r*a.*(pop-sigma*Y(:,2)) - c*Y(:,2));
        term_times(count)=T(end);
        payoffs(count) = payoff;
        count = count+1;
    end
end
figure()
plot(payoffs,term_times)

options2 = odeset('Events', @(t,y)term_cond2(t, y, 115),'RelTol',1e-8);
term_states = [];
payoffs = [];
count = 1;
for l0 = linspace(-50,10,5000)
    y0=[l0,I0];
    [T,Y] = ode45(f,[0, Tend],y0, options2);
    if T(end)<1000
        [a,~,~] = get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
        payoff=trapz(T,r*a.*(pop-sigma*Y(:,2)) - c*Y(:,2));
        term_states(count)=Y(end,2);
        payoffs(count) = payoff;
        count = count+1;
    end
end
figure()
plot(payoffs,term_states)


function [value, isterminal, direction] = term_cond(T, Y,term_I)
value= abs(Y(2)-term_I)<=1e-5;
isterminal = 1;   % Stop the integration
direction  = 0;
end

function [value, isterminal, direction] = term_cond2(T, Y,term_time)
value= abs(T-term_time)<=1e-5;
isterminal = 1;   % Stop the integration
direction  = 0;
end
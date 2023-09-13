close all; clear all;

%%
pop=1;
beta = 0.06;
r = 0.1;
c = 0.2;
gamma = 1/110;
amax = 0.9;
amin = 0.2;
I0 = 0.99;
%pop=1;
sigma = 0.5;
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
t =0:0.01:1000;
ivals = linspace(0,1,length(t));
figure()
plot(ivals, -r*(pop - sigma*ivals)./(beta*(pop-ivals).*ivals),'linewidth',2)
ylim([-40,0])
hold on

%%
[term_times,lvals] = find_term_times(200,130,beta,r,c,gamma,amax,amin,I0,pop,sigma);

%% generates vector plot
Opt1    = odeset('Events', @(t,y)lzero(t, y),'RelTol',1e-8);
lvals_vec = linspace(-1,-15,20);
for i =1:length(lvals_vec)
    [T,Y] = ode45(f,t,[lvals_vec(i),I0], Opt1);
    plot(Y(:,2),Y(:,1),'linewidth',2, 'color','r')
    ylim([-40,0])
end
%%
i_sing = pop - sqrt(pop*r*gamma*(1-sigma)/(c*beta));
l_sing = -r*(pop-sigma*i_sing)/(beta*(pop-i_sing)*i_sing);
plot(i_sing, l_sing, '*','linewidth',2, 'color','black', 'markersize',5)
%%
r1 = shooting_switch(115,-4.89,-4.88,beta,r,c,gamma,amax,amin,I0,pop,sigma);
%%
r2 = shooting_switch(115,lvals(10),lvals(8),beta,r,c,gamma,amax,amin,I0,pop,sigma);
%%
r3 = shooting_switch(115,lvals(4),-7.11,beta,r,c,gamma,amax,amin,I0,pop,sigma);

%%
ivals = linspace(0,1,length(t));
plot(ivals, -r*(pop - sigma*ivals)./(beta*(pop-ivals).*ivals),'linewidth',2)
ylim([-20,1])
hold on
[~,Y] = ode45(f,t,[r1,I0]);
Y(end,1)
plot(Y(:,2),Y(:,1),'linewidth',2)

[~,Y] = ode45(f,t,[r2,I0]);
Y(end,1)
plot(Y(:,2),Y(:,1),'linewidth',2)

[~,Y] = ode45(f,t,[r3,I0]);
Y(end,1)
plot(Y(:,2),Y(:,1),'linewidth',2)
legend('','\lambda =-4.8884','\lambda =-7.0806','\lambda = -7.1147')
%%
% [switches,indexes]= bisection3(term_times,switch_times,Final_time);
% objf2 = zeros(length(switches),1);
% count=1;
% for i = 1:length(switches)
%     t1 = switches(1,i);
%     t2 = switches(2,i);
%     a_vobj = zeros(length(trange),1);
%     I_vobj = zeros(length(trange),1);
%     
%     a_vobj(trange<=t1)=amax*ones(length(a_vobj(trange<=t1)),1);
%     a_vobj(trange>t1 & trange<=t2)=amin*ones(length(a_vobj(trange>t1 & trange<=t2)),1);
%     a_vobj(trange>t2)=amax*ones(length(a_vobj(trange>t2)),1);
%     
%     I_1= @(t) dM*I0./(amax*beta*I0 + (dM - amax*beta*I0)*exp(-dM*t));
%     I_2= @(t) dm*I_1(t1)./(amin*beta*I_1(t1) + (dm - amin*beta*I_1(t1))*exp(-dm*(t-t1)));
%     I_3= @(t) dM*I_2(t2)./(amax*beta*I_2(t2) + (dM - amax*beta*I_2(t2))*exp(-dM*(t-t2)));
% 
% 
%     I_vobj(trange<=t1) = I_1(trange(trange<=t1));
%     I_vobj(trange>t1 & trange<=t2) = I_2(trange(trange>t1&trange<=t2));
%     I_vobj(trange>t2) = I_3(trange(trange>t2));
%                               
%     objf2(count) = trapz(r*a_vobj.*(pop-sigma*I_vobj) - c*I_vobj);
%     count=count+1;
% end
% 
% %% Root Finding
% [switches,~]= bisection3(term_times,switch_times,Final_time)
% 
% function [switch_time,mid_time] = bisection(a,b,term_times,switch_times,time_val)
%    
%     a_time = term_times(round(a));
%     b_time = term_times(round(b));
%     xmid = (a+b)/2;
%     mid_time = term_times(round(xmid));
%     switch_time=switch_times(:,round(xmid));
%     
%     while abs(a-b)>1e-8
%         if (mid_time-time_val) == 0
%             break
%         elseif sign(a_time-time_val) == sign(b_time-time_val)
%             a = xmid;
%         else
%             b = xmid;
%         end
%     a_time = term_times(round(a));
%     b_time = term_times(round(b));
%     xmid = (a+b)/2;
%     mid_time = term_times(round(xmid));
%     switch_time=switch_times(:,round(xmid));
%     end
% end
% 
% function [switch_time,idxs] = bisection3(term_times,switch_times,time_val)
%    switch_time =zeros(2,3); 
%    idxs = zeros(1,3);
%    a = 1;
%    b = a+1;
%    a_time = term_times(round(a));
%    b_time = term_times(round(b));
%    while sign(a_time-time_val) == sign(b_time-time_val)
%        b = b+1;
%        b_time = term_times(round(b));
%    end
%    
%    [switch_time(:,1),idxs(1)] = bisection(a,b,term_times,switch_times,time_val);
%    a = b+1;
%    b = a+1;
%    a_time = term_times(round(a));
%    b_time = term_times(round(b));
%    while sign(a_time-time_val) == sign(b_time-time_val)
%        b = b+1;
%        b_time = term_times(round(b));
%    end
%    [switch_time(:,2),idxs(2)] = bisection(a,b+1,term_times,switch_times,time_val);
%    a = b+1;
%    b = a+1;
%    a_time = term_times(round(a));
%    b_time = term_times(round(b));
%    while sign(a_time-time_val) == sign(b_time-time_val)
%        b = b+1;
%        b_time = term_times(round(b));
%    end
%    [switch_time(:,3),idxs(3)] = bisection(a,b,term_times,switch_times,time_val);
% end

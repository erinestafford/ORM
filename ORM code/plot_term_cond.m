%%
close all; clear all;
global pop beta r c gamma amax amin I0 sigma Tend term_I
pop =1;
beta = 0.06;%0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
I0 = 0.7;
sigma = 0.5;
Tend = 1000;

f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
t =0:0.01:1000;
ivals = linspace(0,pop,length(t));
figure()
%ylim([-40,0])
%hold on;
plot(ivals, (c*ivals-amax*r*(pop-sigma*ivals))./(ivals.*(amax*beta*(pop*ivals)-gamma)),'linewidth',2,'color','red')
hold on;
%plot(ivals, (c*ivals+amin*r*(pop-sigma*ivals))./(ivals.*(amin*beta*(pop*ivals)-gamma)),'linewidth',2,'color','blue')
ylim([-40,30])
xlim([0,1])
plot(ivals, -r*(pop - sigma*ivals)./(beta*(pop-ivals).*ivals),'linewidth',2,'color','black')
plot(ivals, zeros(length(ivals),1),'LineStyle','--','color',[0.7 0.7 0.7])
%legend('Terminal condition, a(T)=M','Terminal condition, a(T)=m','Switching curve')

%%
options = odeset('Events', @(t,y)term_cond(t, y),'RelTol',1e-8);

f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
y0=[-20,I0];
[T,Y] = ode45(f,[0, Tend],y0, options);
plot(Y(:,2),Y(:,1),'linewidth',2, 'LineStyle',':', 'Color','k')
[a,~,~] = get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
payoff=trapz(T,r*a.*(pop-sigma*Y(:,2)) - c*Y(:,2));
text_x1 = Y(end,2);
text_y1 = Y(end,1)+1;
text_x2 = Y(1,2);
text_y2 = Y(1,1);
label = sprintf('%.1f',payoff);
label2 = sprintf('%.1f',T(end));
text(text_x1, text_y1, label, 'Color', 'k', 'FontSize',8);
text(text_x2, text_y2,label2, 'Color', 'k', 'FontSize',8);
term_I = Y(end,2);

options = odeset('Events', @(t,y)term_cond2(t, y),'RelTol',1e-8);

for l0 = linspace(-30,0,5000)
    y0=[l0,I0];
    [T,Y] = ode45(f,[0, Tend],y0, options);
    if T(end)<1000
        plot(Y(:,2),Y(:,1),'linewidth',1.5)
        [a,~,~] = get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
        payoff=trapz(T,r*a.*(pop-sigma*Y(:,2)) - c*Y(:,2));
        text_x1 = Y(end,2);
        text_y1 = Y(end,1);
        text_x2 = Y(1,2);
        text_y2 = Y(1,1);
        label = sprintf('%.1f',payoff);
        label2 = sprintf('%.1f',T(end));
        text(text_x1, text_y1, label, 'Color', 'k', 'FontSize',8);
        text(text_x2, text_y2,label2, 'Color', 'k', 'FontSize',8);
    end
end
% [t,I,curlambda,a] = opSIS_Icost(I0,beta,gamma,r,c,sigma,pop,Tend, amin, amax);
% plot(I,curlambda,'linewidth',2.5, 'Color','black')
% payoff=trapz(t,r*a.*(pop-sigma*I) - c*I);
% text_x1 = I(end);
% text_y1 = curlambda(end)+1;
% label = sprintf('%.1f',payoff);
% text(text_x1, text_y1, label, 'Color', 'k', 'FontSize',8);


function [value, isterminal, direction] = term_cond(T, Y)
global pop beta r c gamma amax amin I0 sigma Tend term_I
value= Y(1)<=(c*Y(2)-amax*r*(pop-sigma*Y(2)))./(Y(2).*(amax*beta*(pop*Y(2))-gamma));
isterminal = 1;   % Stop the integration
direction  = 0;
end

function [value, isterminal, direction] = term_cond2(T, Y)
global pop beta r c gamma amax amin I0 sigma Tend term_I
value= abs(Y(2)-term_I)<=1e-5;
isterminal = 1;   % Stop the integration
direction  = 0;
end
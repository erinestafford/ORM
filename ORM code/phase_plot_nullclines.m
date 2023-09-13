close all; clear all;
pop=100;
beta = 0.0232/pop;
r = 0.1;
c = 0.2;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
I0 = 0.99*pop;
sigma = 0.1;
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
t =0:1:1000;
ivals = linspace(0,pop,length(t));
plot(ivals, -r*(pop - sigma*ivals)./(beta*(pop-ivals).*ivals),'linewidth',2,'color','k')
ylim([-40,0])
xlim([0,100])
hold on
I_amax_null = pop - gamma/(amax*beta);
lamb_amax_null =-r*(pop - sigma*I_amax_null)/(beta*(pop-I_amax_null)*I_amax_null);
plot([1 1]*I_amax_null, [lamb_amax_null 0],'linewidth',2,'color','r')    

I_amin_null = pop - gamma/(amin*beta);
lamb_amin_null =-r*(pop - sigma*I_amin_null)/(beta*(pop-I_amin_null)*I_amin_null);
plot([1 1]*I_amin_null, [-20 lamb_amin_null],'linewidth',2,'color','r')    
amax_i_end=(c*pop*beta+r*gamma*sigma-2*pop*r*beta*amax+sqrt((c*pop*beta+r*gamma*sigma-2*pop*r*beta*amax)^2+4*pop*r*beta*(gamma - pop*beta*amax)*(-c + r*sigma*amax)))/(2*beta*(c-r*sigma*amax));
ivals_amax = linspace(amax_i_end, pop,length(t));
plot(ivals_amax, -(amax*r*sigma+c)./(gamma - amax*beta*(pop-2*ivals_amax)),'linewidth',2,'color','b')

amin_i_end=(c*pop*beta+r*gamma*sigma-2*pop*r*beta*amin+sqrt((c*pop*beta+r*gamma*sigma-2*pop*r*beta*amin)^2+4*pop*r*beta*(gamma - pop*beta*amin)*(-c + r*sigma*amin)))/(2*beta*(c-r*sigma*amin));
ivals_amin = linspace(0,amin_i_end,length(t));
plot(ivals_amin, -(amin*r*sigma+c)./(gamma - amin*beta*(pop-2*ivals_amin)),'linewidth',2,'color','b')

%%
Opt1    = odeset('Events', @(t,y)lzero(t, y),'RelTol',1e-8,'AbsTol',1e-9);
lvals=-15:1:-1;
c=1;
for i =1:length(lvals)
    [T,Y] = ode45(f,t,[lvals(i),I0], Opt1);
    plot(Y(:,2),Y(:,1),'linewidth',1, 'color','k','linewidth',1.5)
    c=c+1;
end
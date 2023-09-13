%%
% pop=100;
% beta = 0.06/pop;%0.0232/pop;%0.06;
% r = 0.1;%0.39;
% c = 0.2;%6.83;
% gamma = 1.0/110.0;
% amax = 0.9;
% amin = 0.2;
% I0 = 0.99*pop;
% sigma = 0.5;

% pop=500;
% beta = 0.0232/pop;
% r = 0.39;
% c = 6.83;
% gamma = 1/110.0;
% amax = 0.9;
% amin = 0.2;
% I0 = 0.99*pop;
% sigma = 0.5;

% pop=100.0;
% beta = 0.0232/pop;
% r = 0.39;
% c = 155.0/36.0;
% gamma = 1.0/36.0;
% amax = 0.75;
% amin = 0.1;
% I0 = 0.99*pop;
% sigma = 0.05;

% pop=100;
% beta = 0.06/pop;%0.1/pop;%0.06/pop;
% r = 0.1;
% c = 0.2;%155/36;%6.39;
% gamma = 1.0/110.0;%1.0/250;%1.0/36.0;
% amax = 0.9;%0.75;
% amin = 0.2;
% I0 = 0.99*pop;
% %https://bmcvetres.biomedcentral.com/articles/10.1186/s12917-016-0905-3
% sigma = 0.05;


pop=100;
beta = 0.0232/pop;%0.1/pop;%0.06/pop;
r = 0.1;
c = 0.2;%155/36;%6.39;
gamma = 1.0/110.0;%1.0/250;%1.0/36.0;
amax = 0.9;%0.75;
amin = 0.2;
I0 = 0.99*pop;
%https://bmcvetres.biomedcentral.com/articles/10.1186/s12917-016-0905-3
sigma = 0.1;

l0 = -1.0;
dy = 0.1;
NSTPS = 2000000;
HT = 500;
%TTARG=250.0;
TTARG=305.0;

f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
t =0:0.01:1000;
ivals = linspace(0,pop,length(t));
figure()
plot(ivals, -r*(pop - sigma*ivals)./(beta*(pop-ivals).*ivals),'linewidth',2,'color','k')
ylim([-40,0])
hold on;
%%
lvals = [];
tts = [];
y2 = l0;
count=1;
while y2>-HT
    tt = get_Ttime(y2,I0,beta,gamma,r,c,sigma,pop,amin,amax);
    if tt<(1e6/2)
        lvals(count) = y2;
        tts(count) = tt;
        count=count+1;
    else
        break
    end
    y2=y2-dy;
    
end
%%
options = odeset('Events', @(t,y)lzero(t, y),'RelTol',1e-8);
found_roots = [];
as=[];
cols = ['r','b','y'];
count=1;
for i = 2:length(tts)
    if ((tts(i)<TTARG)&&(TTARG<tts(i-1))) ||  ((tts(i)>TTARG)&&(TTARG>tts(i-1)))
        lval=fzero(@(x)rf_sa(x,I0,beta,gamma,r,c,sigma,pop,amin,amax,TTARG),[lvals(i),lvals(i-1)]);
        y0=[lval,I0];
        tstart = 0.0;
        tend=TTARG;
        [T,Y] = ode45(f,[tstart, tend],y0,options);
        [~,~,switch_ind]=get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
        T(switch_ind)
        plot(Y(:,2),Y(:,1),'linewidth',1.5, 'color',cols(count))
        count=count+1;
        found_roots=[found_roots;lval];
    end
        
end
%legend('switching curve',append('\lambda = ',num2str(found_roots(1))),append('\lambda = ',num2str(found_roots(2))),append('\lambda = ',num2str(found_roots(3))),'Location','northwest')
%%

y0=[found_roots(1),I0];
tstart = 0.0;
tend=TTARG;
[T,Y] = ode45(f,[tstart, tend],y0,options);
[~,~,switch_ind]=get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
T(switch_ind)
plot(Y(:,2),Y(:,1),'linewidth',1.5, 'color',cols(1))
        
%%
[~,Y] = ode45(f,[0, TTARG],[found_roots(1),I0],options);
[a1,num_switch1]=get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
ylim([-40,0])
%%

I_amax_null = pop - gamma/(amax*beta);
lamb_amax_null =-r*(pop - sigma*I_amax_null)/(beta*(pop-I_amax_null)*I_amax_null);
plot([1 1]*I_amax_null, [lamb_amax_null 0],'linewidth',1.5,'color',[.7 .7 .7])    
%%
I_amin_null = pop - gamma/(amin*beta);
lamb_amin_null =-r*(pop - sigma*I_amin_null)/(beta*(pop-I_amin_null)*I_amin_null);
plot([1 1]*I_amin_null, [-40 lamb_amin_null],'linewidth',1.5,'color',[.7 .7 .7])    
%%
amax_i_end=(c*pop*beta+r*gamma*sigma-2*pop*r*beta*amax+sqrt((c*pop*beta+r*gamma*sigma-2*pop*r*beta*amax)^2+4*pop*r*beta*(gamma - pop*beta*amax)*(-c + r*sigma*amax)))/(2*beta*(c-r*sigma*amax));
ivals_amax = linspace(amax_i_end, pop,length(t));
plot(ivals_amax, -(amax*r*sigma+c)./(gamma - amax*beta*(pop-2*ivals_amax)),'linewidth',1.5,'color',[.7 .7 .7])
%%
amin_i_end=(c*pop*beta+r*gamma*sigma-2*pop*r*beta*amin+sqrt((c*pop*beta+r*gamma*sigma-2*pop*r*beta*amin)^2+4*pop*r*beta*(gamma - pop*beta*amin)*(-c + r*sigma*amin)))/(2*beta*(c-r*sigma*amin));
ivals_amin = linspace(14.25,amin_i_end,length(t));
plot(ivals_amin, -(amin*r*sigma+c)./(gamma - amin*beta*(pop-2*ivals_amin)),'linewidth',1.5,'color',[.7 .7 .7])

%%
Opt1    = odeset('Events', @(t,y)lzero(t, y),'RelTol',1e-8);
lvals_vec = -30;%linspace(-1,-100,30);
t =0:0.01:6000;
for i =1:length(lvals_vec)
    [T,Y] = ode23s(f,t,[lvals_vec(i),I0], Opt1);
    plot(Y(:,2),Y(:,1),'linewidth',1, 'color','black')
    ylim([-100,0])
    xlim([0,100])
end
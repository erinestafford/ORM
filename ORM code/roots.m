%% Run to get the plot of the candidate solutions in the T-Lambda plane
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

pop=100.0;
beta = 0.0232/pop;
r = 0.1;
c = 0.23;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
I0 = 0.99*pop;
sigma = 0.1;

% pop=100;
% beta = 0.06/pop;%0.1/pop;%0.06/pop;
% r = 1.0;
% c = 3.0;%155/36;%6.39;
% gamma = 1.0/110.0;%1.0/250;%1.0/36.0;
% amax = 0.9;%0.75;
% amin = 0.2;
% I0 = 0.99*pop;
% %https://bmcvetres.biomedcentral.com/articles/10.1186/s12917-016-0905-3
% sigma = 0.05;

l0 = -1.0;
dy = 0.01;
NSTPS = 2000000;
HT = 500;
TTARG=305.0;

f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
t =0:0.01:1000;
ivals = linspace(0,pop,length(t));
figure()
plot(ivals, -r*(pop - sigma*ivals)./(beta*(pop-ivals).*ivals),'linewidth',2,'color','k')
ylim([-50,0])
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
        as = [as;[T(switch_ind)]];
        plot(Y(:,2),Y(:,1),'linewidth',1.5)
        count=count+1;
        found_roots=[found_roots;lval];
    end
        
end
%legend('switching curve',append('\lambda = ',num2str(found_roots(1))),append('\lambda = ',num2str(found_roots(2))),append('\lambda = ',num2str(found_roots(3))),'Location','northwest')


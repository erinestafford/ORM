close all; clear all;
pop=100;
beta = 0.0232/pop;%0.1/pop;%0.06/pop;
r = 0.39;
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7760962/
c = 155.0/30.0;%6.39;
%https://pdf.sciencedirectassets.com/279785/1-s2.0-S0022030295X71904/1-s2.0-S0022030295768643/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEJX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJGMEQCIE7Rvahno34ST8O2TstkfRqB%2F3aFYJbev3XeW70r%2BP0rAiA1wOwXfbMRQ4YzzvPI2CtTr%2BNTC4mpRbFJKxQbj7P1Ayq7BQje%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAUaDDA1OTAwMzU0Njg2NSIMu70TMdmlgo5c7WeWKo8FtsFr95EV1Vi06kVBiii9tyVsTUsGaxj9rN0b8f2sxrD5HszTs%2BEijkkCPA321impliC1hpEnEsN4wQFM%2BVUR7W6l1SwgI6iUrMCkBsjXtxMWGnWD%2FipdOjD%2BWadCBIBjlDVXAP6zD7iO%2FTO%2Bcg%2B%2FaKCmJHLC8AnXAD%2Bs5aDZxlUrZdIYAhRc%2F%2FYmT4SvurzrhDjGY12MK4iTmM%2F0jTzrrzHKrKF65d4CYyvri6V57f5ukVWKgsB1xBhtWou65PNHTF6WFq341nJpjWPRvNKdda4%2FAH42SY4NTtDuzzAHsoADcZE7miR2rq3EC8yJk%2F1CQJAcBVcug%2FaDKCI2Lmj%2Fo7gdN6peLkprPIck0HdayzXXzxvx2h30BKbjvb%2FeA5oe%2BkahFYoXYLVhqZI7v%2F9tFYeydLIm332AQ2UoRnw2eqyBZv95e%2FZQeiClxeY8JHC2FpREERWJkmfk4aGbWRn7NaisABVjVCvCfAMyFPxZ8tHfTYKK37iOl0%2BeGeDh50T%2BJw2vlp5oOj6vfboir2o0srMbDlCd9K6v%2FGKRcXZAzpWtu3EDMd9ygW9n%2FbicQceDs5cRsh%2FFSJpfqhT7TcxEVx6aOabK7mJoARN1cqPC9iJ5jhdn8dzg%2BPpZvUpdIXDtuUHj%2FXZR98WQzI5NZVFSk5o7Q30MdqTq9xDzs4hfqxItK5gGggWOet6J2MAuRN6vgIGc%2Fjo%2FipgeW6y7hz6na%2FKDP%2BxbBcWv%2Fq2Z4aX6HUksT%2BcKhVcdIkSUQgipOVtlcIz8vyFh8SKSsVrw1h3U4Y2W2xQ0tY9igNovimJXpi4giQe9dUcnrEQdVt%2F07vOKMokjM3fn1DiQHU5u3OmepAi0j2Fk0FlmI3mXpXQeNzCUgJ6kBjqyAQXhjKTTBIuW6JkSLxdtMOZNhhBxOq5daeV8U407PXdU890EDFcgXFHji9mco4yFbGDTU1Jlvel5m4oh%2BpbjWl%2FSTEUYY4ng050zkD7f9lMX2GxJbhOffHcsUGdS3TmDfFdUI%2FrGpNKeY0C2QrG2un3i7vK4m7PJt9zA28xvvZLJAZTm150cKrzEiCNz3UxdThnNk0Nn%2FKV4EA1MR3cyfZahiZJX8KXKJ03TL7r35NW7kbU%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230612T212630Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY5FHNTY4V%2F20230612%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=d541a27b71a19fc4cf09d6ba3c3648fbc33a16f3a505bb4857097c721f5067ef&hash=2744b5c125c1f1c509607135424370037b657a630f5b938e8b158a34c7dd4e55&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0022030295768643&tid=spdf-b2f176f8-4ab2-40c1-b408-c452271d030a&sid=8a867a729330d14c5d6b0f88163eb45e6b3dgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f1550025553560c565b09&rr=7d652c25feaf30ab&cc=us
gamma = 1.0/30.0;%1.0/36.0;%1.0/250;%1.0/36.0;
amax = 0.9;%0.75;
amin = 0.1;
I0 = 0.99*pop;
%https://bmcvetres.biomedcentral.com/articles/10.1186/s12917-016-0905-3
sigma = 0.05;
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
t =0:0.01:1000;
ivals = linspace(0,pop,length(t));
figure()
plot(ivals, -r*(pop - sigma*ivals)./(beta*(pop-ivals).*ivals),'linewidth',2,'color','k')
ylim([-400,0])
xlim([0,pop])
hold on
%nullclines
I_amax_null = pop - gamma/(amax*beta);
lamb_amax_null =-r*(pop - sigma*I_amax_null)/(beta*(pop-I_amax_null)*I_amax_null);
plot([1 1]*I_amax_null, [lamb_amax_null 0],'linewidth',2,'color','r')    
I_amin_null = pop - gamma/(amin*beta);
lamb_amin_null =-r*(pop - sigma*I_amin_null)/(beta*(pop-I_amin_null)*I_amin_null);
plot([1 1]*I_amin_null, [-400 lamb_amin_null],'linewidth',2,'color','r')    
amax_i_end=(c*pop*beta+r*gamma*sigma-2*pop*r*beta*amax+sqrt((c*pop*beta+r*gamma*sigma-2*pop*r*beta*amax)^2+4*pop*r*beta*(gamma - pop*beta*amax)*(-c + r*sigma*amax)))/(2*beta*(c-r*sigma*amax));
ivals_amax = linspace(amax_i_end, pop,length(t));
plot(ivals_amax, -(amax*r*sigma+c)./(gamma - amax*beta*(pop-2*ivals_amax)),'linewidth',2,'color','b')

amin_i_end=(c*pop*beta+r*gamma*sigma-2*pop*r*beta*amin+sqrt((c*pop*beta+r*gamma*sigma-2*pop*r*beta*amin)^2+4*pop*r*beta*(gamma - pop*beta*amin)*(-c + r*sigma*amin)))/(2*beta*(c-r*sigma*amin));
ivals_amin = linspace(0,amin_i_end,length(t));
plot(ivals_amin, -(amin*r*sigma+c)./(gamma - amin*beta*(pop-2*ivals_amin)),'linewidth',2,'color','b')

%intersections
-(amin*r*sigma+c)/(amin*beta*pop - gamma)
-(amax*r*sigma+c)/(amax*beta*pop - gamma)
%%
Opt1    = odeset('Events', @(t,y)lzero(t, y),'RelTol',1e-8,'AbsTol',1e-9);
lvals=-400:10:-1;
%term_times=[];
c=1;
for i =1:length(lvals)
    [T,Y] = ode45(f,t,[lvals(i),I0], Opt1);
    plot(Y(:,2),Y(:,1),'linewidth',1, 'color','k','linewidth',1.5)
    %term_times(c)=T(end);
    c=c+1;
end

%term_times=term_times';
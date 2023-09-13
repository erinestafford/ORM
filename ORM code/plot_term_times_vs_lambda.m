clear all; close all; clc;
%%
% pop=1;
% beta = 0.06;%0.0232/pop;%0.06;
% r = 0.1;%0.39;
% c = 0.2;%6.83;
% gamma = 1/110;
% amax = 0.9;
% amin = 0.2;
% I0 = 0.99;
% sigma = 0.5;

pop=100;
beta = 0.4/pop;%0.1/pop;%0.06/pop;
r = 0.39;
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7760962/
c = 155/36;%6.39;
%https://pdf.sciencedirectassets.com/279785/1-s2.0-S0022030295X71904/1-s2.0-S0022030295768643/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEJX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJGMEQCIE7Rvahno34ST8O2TstkfRqB%2F3aFYJbev3XeW70r%2BP0rAiA1wOwXfbMRQ4YzzvPI2CtTr%2BNTC4mpRbFJKxQbj7P1Ayq7BQje%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAUaDDA1OTAwMzU0Njg2NSIMu70TMdmlgo5c7WeWKo8FtsFr95EV1Vi06kVBiii9tyVsTUsGaxj9rN0b8f2sxrD5HszTs%2BEijkkCPA321impliC1hpEnEsN4wQFM%2BVUR7W6l1SwgI6iUrMCkBsjXtxMWGnWD%2FipdOjD%2BWadCBIBjlDVXAP6zD7iO%2FTO%2Bcg%2B%2FaKCmJHLC8AnXAD%2Bs5aDZxlUrZdIYAhRc%2F%2FYmT4SvurzrhDjGY12MK4iTmM%2F0jTzrrzHKrKF65d4CYyvri6V57f5ukVWKgsB1xBhtWou65PNHTF6WFq341nJpjWPRvNKdda4%2FAH42SY4NTtDuzzAHsoADcZE7miR2rq3EC8yJk%2F1CQJAcBVcug%2FaDKCI2Lmj%2Fo7gdN6peLkprPIck0HdayzXXzxvx2h30BKbjvb%2FeA5oe%2BkahFYoXYLVhqZI7v%2F9tFYeydLIm332AQ2UoRnw2eqyBZv95e%2FZQeiClxeY8JHC2FpREERWJkmfk4aGbWRn7NaisABVjVCvCfAMyFPxZ8tHfTYKK37iOl0%2BeGeDh50T%2BJw2vlp5oOj6vfboir2o0srMbDlCd9K6v%2FGKRcXZAzpWtu3EDMd9ygW9n%2FbicQceDs5cRsh%2FFSJpfqhT7TcxEVx6aOabK7mJoARN1cqPC9iJ5jhdn8dzg%2BPpZvUpdIXDtuUHj%2FXZR98WQzI5NZVFSk5o7Q30MdqTq9xDzs4hfqxItK5gGggWOet6J2MAuRN6vgIGc%2Fjo%2FipgeW6y7hz6na%2FKDP%2BxbBcWv%2Fq2Z4aX6HUksT%2BcKhVcdIkSUQgipOVtlcIz8vyFh8SKSsVrw1h3U4Y2W2xQ0tY9igNovimJXpi4giQe9dUcnrEQdVt%2F07vOKMokjM3fn1DiQHU5u3OmepAi0j2Fk0FlmI3mXpXQeNzCUgJ6kBjqyAQXhjKTTBIuW6JkSLxdtMOZNhhBxOq5daeV8U407PXdU890EDFcgXFHji9mco4yFbGDTU1Jlvel5m4oh%2BpbjWl%2FSTEUYY4ng050zkD7f9lMX2GxJbhOffHcsUGdS3TmDfFdUI%2FrGpNKeY0C2QrG2un3i7vK4m7PJt9zA28xvvZLJAZTm150cKrzEiCNz3UxdThnNk0Nn%2FKV4EA1MR3cyfZahiZJX8KXKJ03TL7r35NW7kbU%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230612T212630Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY5FHNTY4V%2F20230612%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=d541a27b71a19fc4cf09d6ba3c3648fbc33a16f3a505bb4857097c721f5067ef&hash=2744b5c125c1f1c509607135424370037b657a630f5b938e8b158a34c7dd4e55&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0022030295768643&tid=spdf-b2f176f8-4ab2-40c1-b408-c452271d030a&sid=8a867a729330d14c5d6b0f88163eb45e6b3dgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f1550025553560c565b09&rr=7d652c25feaf30ab&cc=us
gamma = 1.0/36.0;%1.0/250;%1.0/36.0;
amax = 1.0;%0.75;
amin = 0.2;
I0 = 0.99*pop;
%https://bmcvetres.biomedcentral.com/articles/10.1186/s12917-016-0905-3
sigma = 0.05;

t =0:0.01:800;
%%
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
options = odeset('Events',@lzero,'MaxStep',.0001);
% [t,Y,te,ye,ie] = ode45(f,[0:0.01:200],[-10,I0],options);
% figure()
% plot(t,Y(:,2))
% Y(end,1)
% t(end)
%%
c = 1;
lvals=-90:.01:-40;
for l=lvals
    [t,Y] = ode45(f,t,[l,I0],options);
    tt(c)=t(end);
    c = c+1;
end
%%
plot(lvals,tt)
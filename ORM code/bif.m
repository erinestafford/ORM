%%
% pop=1;
% beta = 0.06;%0.0232/pop;%0.06;
% r = 0.1;%0.39;
% c = 0.2;%6.83;
% gamma = 1/110;
% amax = 0.9;
% amin = 0.2;
% I0 = 0.99;

pop=100;
beta = 0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1/110;
amax = 0.9;
amin = 0.2;
I0 = 0.99;
ds = 0.001;
smin=0.45;
smax=0.7;

l0 = -1.0;
dy = 0.01;
HT = 100;
TTARG=305.0;

dl = 0.1;

f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);

sigma = smin;
sigmas = [];
roots=[];
j = 1;
sigmas(j)=sigma;
while sigma<smax
    y2s = [];
    tts=[];
    c=1;
    y2 = l0;
    while y2>-HT
        y2=y2-dl;
        tt = ttime_sigma(y2,sigma);
        if tt<(1e6/2)
            lvals(c) = y2;
            tts(c) = tt;
            c=c+1;
        else
            break
        end
    end
    c=1;
    roots_temp = [0,0,0];
    for i = 2:length(tts)
    if ((tts(i)<TTARG)&&(TTARG<tts(i-1))) ||  ((tts(i)>TTARG)&&(TTARG>tts(i-1)))
        lval=bisect_ttarg_sigma(@rf_sigma,lvals(i),lvals(i-1),1e-8,TTARG,sigma);
        roots_temp(c)=lval;
        c=c+1;
    end
    end
    
    
    roots=[roots;roots_temp];
    sigmas(j)=sigma;
    j = j+1;
  
    
    sigma=sigma+ds;
    
end



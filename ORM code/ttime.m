function term = ttime(l0)
NSTPS=2000000;HT=20;amin=0.2;amax=0.9;r=0.1;sigma=0.4;beta=0.06;pop=1;c=0.2;gamma=1/110.0;
 options = odeset('Events',@lzero,'MaxStep',.1,'RelTol',1.0e-6,'AbsTol',1.0e-6);
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
tstart = 0.0;
term=1e6;
y0=[l0,0.99];
tend=tstart;
for i=1:NSTPS
    tend = tend+0.1;
    [~,Y] = ode45(f,[tstart, tend],y0,options);
    y0 = Y(end,:);
    tstart=tend;
    if Y(end,1)>0
        term = tend;
        break;
    end
    if Y(end,1)<-HT
        term = 1e6;
    end
end
end